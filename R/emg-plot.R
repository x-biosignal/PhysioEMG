# Declare ggplot2 non-standard-evaluation column names (avoids R CMD check's
# "no visible binding for global variable" notes).
utils::globalVariables(c(
  "time", "value", "series", "channel", "freq", "synergy", "muscle", "weight",
  "activation", "n_synergies", "vaf", "sample", "ymax", "lo", "hi", "cci",
  "label", "selected"))

.emg_channel_index <- function(x, channel) {
  n <- ncol(SummarizedExperiment::assay(x, defaultAssay(x)))
  if (is.numeric(channel)) {
    idx <- as.integer(channel)
  } else {
    labels <- tryCatch(as.character(SummarizedExperiment::colData(x)$label),
                       error = function(e) NULL)
    idx <- match(as.character(channel), labels)
  }
  if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > n) {
    stop("channel out of range: ", channel, call. = FALSE)
  }
  idx
}

#' Plot EMG Signal and Envelope
#'
#' Plots the raw signal and its amplitude envelope for one channel, optionally
#' shading the resting noise floor from [emgQualityCheck()].
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param channel Channel index or label to plot (default: 1).
#' @param method Envelope method passed to [emgEnvelope()] (default: "rms").
#' @param window_ms Envelope window in ms (default: 50).
#' @param qc If TRUE (default), shade the resting noise floor and annotate SNR.
#' @param assay_name Input assay name (default: first assay).
#' @return A ggplot object.
#' @seealso [emgEnvelope()], [emgQualityCheck()]
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @export
#' @examples
#' pe <- make_emg(n_time = 2000, n_channels = 1)
#' plotEmgEnvelope(pe)
plotEmgEnvelope <- function(x, channel = 1, method = "rms", window_ms = 50,
                            qc = TRUE, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  ch <- .emg_channel_index(x, channel)
  sr <- samplingRate(x)
  sig <- SummarizedExperiment::assay(x, assay_name)[, ch]
  env <- SummarizedExperiment::assay(
    emgEnvelope(x, method = method, window_ms = window_ms,
                assay_name = assay_name), "envelope")[, ch]
  tt <- (seq_along(sig) - 1) / sr

  df <- rbind(
    data.frame(time = tt, value = sig, series = "raw"),
    data.frame(time = tt, value = env, series = "envelope"))

  p <- ggplot2::ggplot(df, ggplot2::aes(time, value, colour = series)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::scale_colour_manual(values = c(raw = "grey70",
                                            envelope = "#1b6ca8")) +
    ggplot2::labs(x = "Time (s)", y = "Amplitude", colour = NULL,
                  title = "EMG signal and amplitude envelope")

  if (qc) {
    qcr <- emgQualityCheck(x, assay_name = assay_name)
    nf <- qcr$baseline_noise[ch]
    if (is.finite(nf)) {
      p <- p + ggplot2::annotate("rect", xmin = -Inf, xmax = Inf,
                                 ymin = 0, ymax = nf, alpha = 0.12,
                                 fill = "#d1495b")
    }
    p <- p + ggplot2::labs(subtitle = sprintf(
      "SNR = %.1f dB%s", qcr$snr_db[ch],
      if (isTRUE(qcr$ecg_contamination[ch])) " | ECG-contaminated" else ""))
  }
  p + PhysioCore::theme_physio()
}

#' Plot EMG Activation Onsets and Offsets
#'
#' Plots the rectified signal with detected activation intervals shaded and
#' onset/offset times marked.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param channel Channel index or label to plot (default: 1).
#' @param method Onset method passed to [emgOnsetDetect()] (default:
#'   "hodges_bui").
#' @param assay_name Input assay name (default: first assay).
#' @param ... Further arguments for [emgOnsetDetect()].
#' @return A ggplot object.
#' @seealso [emgOnsetDetect()]
#' @references Hodges, P.W. & Bui, B.H. (1996). "A comparison of computer-based
#'   methods for the determination of onset of muscle contraction using
#'   electromyography." Electroencephalography and Clinical Neurophysiology,
#'   101(6), 511-519. doi:10.1016/S0921-884X(96)95190-5
#' @export
#' @examples
#' pe <- make_emg_contraction()
#' plotEmgOnset(pe)
plotEmgOnset <- function(x, channel = 1, method = "hodges_bui",
                         assay_name = NULL, ...) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  ch <- .emg_channel_index(x, channel)
  sr <- samplingRate(x)
  sig <- abs(SummarizedExperiment::assay(x, assay_name)[, ch])
  tt <- (seq_along(sig) - 1) / sr

  det <- emgOnsetDetect(x, method = method, assay_name = assay_name, ...)
  on <- det$onsets[det$onsets$channel == ch, , drop = FALSE]
  off <- det$offsets[det$offsets$channel == ch, , drop = FALSE]

  p <- ggplot2::ggplot(data.frame(time = tt, value = sig),
                       ggplot2::aes(time, value)) +
    ggplot2::geom_line(colour = "grey40", linewidth = 0.3)

  n_ev <- min(nrow(on), nrow(off))
  if (n_ev > 0) {
    bands <- data.frame(lo = on$time_sec[seq_len(n_ev)],
                        hi = off$time_sec[seq_len(n_ev)])
    p <- p +
      ggplot2::geom_rect(data = bands, inherit.aes = FALSE,
                         ggplot2::aes(xmin = lo, xmax = hi,
                                      ymin = -Inf, ymax = Inf),
                         alpha = 0.15, fill = "#66a182") +
      ggplot2::geom_vline(xintercept = on$time_sec[seq_len(n_ev)],
                          colour = "#2e8540", linewidth = 0.4) +
      ggplot2::geom_vline(xintercept = off$time_sec[seq_len(n_ev)],
                          colour = "#d1495b", linewidth = 0.4,
                          linetype = "dashed")
  }
  p + ggplot2::labs(x = "Time (s)", y = "Rectified amplitude",
                    title = "EMG activation onsets and offsets") +
    PhysioCore::theme_physio()
}

#' Plot EMG Spectral-Fatigue Regression
#'
#' Scatters the windowed median or mean frequency against time with the fitted
#' fatigue-regression line and slope annotation.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param feature "mdf" (median frequency) or "mnf" (mean frequency).
#' @param window_sec,overlap Windowing passed to [emgFatigue()].
#' @param assay_name Input assay name (default: first assay).
#' @return A ggplot object (faceted by channel).
#' @seealso [emgFatigue()], [emgFatigueSlope()]
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
#' @examples
#' pe <- make_emg_fatigue()
#' plotEmgFatigue(pe)
plotEmgFatigue <- function(x, feature = c("mdf", "mnf"), window_sec = 1.0,
                           overlap = 0.5, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  feature <- match.arg(feature)
  col <- if (feature == "mdf") "median_freq" else "mean_freq"

  fat <- emgFatigue(x, window_sec = window_sec, overlap = overlap,
                    assay_name = assay_name)
  slp <- emgFatigueSlope(x, feature = feature, window_sec = window_sec,
                         overlap = overlap, assay_name = assay_name)
  df <- data.frame(time = fat$time_sec, freq = fat[[col]],
                   channel = factor(fat$channel))

  ann <- data.frame(
    channel = factor(slp$channel),
    label = sprintf("slope = %.1f Hz/min\nR2 = %.2f, p = %.3g",
                    slp$slope_hz_per_min, slp$r_squared, slp$p_value))

  ggplot2::ggplot(df, ggplot2::aes(time, freq)) +
    ggplot2::geom_point(colour = "#1b6ca8", alpha = 0.6, size = 1) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = "#d1495b", linewidth = 0.6) +
    ggplot2::geom_text(data = ann, inherit.aes = FALSE,
                       ggplot2::aes(x = Inf, y = Inf, label = label),
                       hjust = 1.05, vjust = 1.3, size = 3) +
    ggplot2::facet_wrap(~ channel, labeller = ggplot2::label_both) +
    ggplot2::labs(x = "Time (s)",
                  y = sprintf("%s frequency (Hz)", toupper(feature)),
                  title = "EMG spectral-fatigue regression") +
    PhysioCore::theme_physio()
}

#' Plot Muscle Synergies
#'
#' Renders a three-panel figure for a muscle-synergy decomposition: the synergy
#' weight matrix W (heatmap), the activation time courses H, and the VAF curve.
#'
#' @param x A [muscleSynergy()] result, or a [muscleSynergyOrder()] result (whose
#'   VAF curve and selected order are used).
#' @return A patchwork/ggplot object combining the three panels.
#' @seealso [muscleSynergy()], [muscleSynergyOrder()]
#' @references Ting, L.H. & Chvatal, S.A. (2010). "Decomposing muscle activity in
#'   motor tasks: methods and interpretation." In Motor Control (Oxford).
#' @export
#' @examples
#' set.seed(1)
#' tt <- seq_len(400)
#' H <- sapply(c(120, 200, 280), function(c0) exp(-((tt - c0) / 30)^2))
#' W <- rbind(c(1, 1, 0, 0, 0, 0), c(0, 0, 1, 1, 0, 0), c(0, 0, 0, 0, 1, 1))
#' data <- H %*% W + matrix(abs(rnorm(400 * 6, sd = 0.02)), 400, 6)
#' pe <- PhysioExperiment(assays = list(env = data), samplingRate = 1000)
#' ord <- muscleSynergyOrder(pe, max_synergies = 5, seed = 1)
#' plotSynergy(ord)
plotSynergy <- function(x) {
  is_order <- is.list(x) && !is.null(x$vaf_curve)
  fit <- if (is_order) x$fit else x
  stopifnot(is.list(fit), all(c("W", "H") %in% names(fit)))
  vaf_curve <- if (is_order) x$vaf_curve else
    data.frame(n_synergies = nrow(fit$W), vaf = fit$vaf)
  selected <- if (is_order) x$selected else nrow(fit$W)

  W <- fit$W; H <- fit$H
  k <- nrow(W); M <- ncol(W)

  wdf <- data.frame(
    synergy = factor(rep(seq_len(k), times = M)),
    muscle = factor(rep(seq_len(M), each = k)),
    weight = as.vector(W))
  pW <- ggplot2::ggplot(wdf, ggplot2::aes(muscle, synergy, fill = weight)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "#1b6ca8") +
    ggplot2::labs(x = "Muscle", y = "Synergy", title = "Weights (W)") +
    PhysioCore::theme_physio()

  hdf <- data.frame(
    sample = rep(seq_len(nrow(H)), times = k),
    synergy = factor(rep(seq_len(k), each = nrow(H))),
    activation = as.vector(H))
  pH <- ggplot2::ggplot(hdf, ggplot2::aes(sample, activation, colour = synergy)) +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::labs(x = "Time (samples)", y = "Activation",
                  title = "Activations (H)") +
    PhysioCore::theme_physio()

  vdf <- data.frame(n_synergies = vaf_curve$n_synergies, vaf = vaf_curve$vaf,
                    selected = vaf_curve$n_synergies == selected)
  pV <- ggplot2::ggplot(vdf, ggplot2::aes(n_synergies, vaf)) +
    ggplot2::geom_line(colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = selected), size = 2) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey50",
                                            `TRUE` = "#d1495b"),
                                 guide = "none") +
    ggplot2::labs(x = "Number of synergies", y = "VAF", title = "VAF curve") +
    PhysioCore::theme_physio()

  patchwork::wrap_plots(pW, pH, pV, ncol = 1)
}

#' Plot Agonist-Antagonist Co-Activation
#'
#' Overlays the agonist and antagonist envelopes and shades their common
#' (co-contraction) area, annotated with the co-contraction index.
#'
#' @param x A PhysioExperiment object of EMG envelopes.
#' @param agonist,antagonist Channel index or label of the two muscles.
#' @param method Co-contraction index for the annotation (see
#'   [emgCoContraction()]).
#' @param normalize Peak-normalize each channel before plotting (default: TRUE).
#' @param assay_name Input assay name (default: first assay).
#' @return A ggplot object.
#' @seealso [emgCoContraction()]
#' @references Falconer, K. & Winter, D.A. (1985). "Quantitative assessment of
#'   co-contraction at the ankle joint in walking." Electromyography and Clinical
#'   Neurophysiology, 25(2-3), 135-149.
#' @export
#' @examples
#' t <- seq_len(1000)
#' a <- exp(-((t - 400) / 120)^2); b <- exp(-((t - 600) / 120)^2)
#' pe <- PhysioExperiment(assays = list(env = cbind(a, b)),
#'   colData = S4Vectors::DataFrame(label = c("TA", "GAS")), samplingRate = 1000)
#' plotCoactivation(pe, "TA", "GAS")
plotCoactivation <- function(x, agonist, antagonist,
                             method = c("falconer_winter", "rudolph", "frost"),
                             normalize = TRUE, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  ai <- .emg_channel_index(x, agonist)
  bi <- .emg_channel_index(x, antagonist)
  sr <- samplingRate(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  a <- abs(data[, ai]); b <- abs(data[, bi])
  if (normalize) {
    if (max(a) > 0) a <- a / max(a)
    if (max(b) > 0) b <- b / max(b)
  }
  tt <- (seq_along(a) - 1) / sr

  cc <- emgCoContraction(x, agonist = agonist, antagonist = antagonist,
                         method = method, normalize = normalize,
                         assay_name = assay_name)

  lab_a <- tryCatch(as.character(SummarizedExperiment::colData(x)$label[ai]),
                    error = function(e) NULL)
  lab_b <- tryCatch(as.character(SummarizedExperiment::colData(x)$label[bi]),
                    error = function(e) NULL)
  if (is.null(lab_a) || is.na(lab_a)) lab_a <- paste0("ch", ai)
  if (is.null(lab_b) || is.na(lab_b)) lab_b <- paste0("ch", bi)

  common <- data.frame(time = tt, ymax = pmin(a, b))
  lines <- rbind(
    data.frame(time = tt, value = a, series = lab_a),
    data.frame(time = tt, value = b, series = lab_b))

  ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = common, inherit.aes = FALSE,
                         ggplot2::aes(x = time, ymin = 0, ymax = ymax),
                         fill = "#edae49", alpha = 0.5) +
    ggplot2::geom_line(data = lines,
                       ggplot2::aes(time, value, colour = series),
                       linewidth = 0.5) +
    ggplot2::labs(x = "Time (s)", y = "Normalized amplitude", colour = NULL,
                  title = "Agonist-antagonist co-activation",
                  subtitle = sprintf("%s CCI = %.1f", method, cc$summary)) +
    PhysioCore::theme_physio()
}
