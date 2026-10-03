#' Detect EMG Muscle Activation Onset and Offset
#'
#' Identifies when muscle activation begins and ends. The Hodges-Bui method
#' rectifies the signal, smooths it with a centered moving average
#' (\code{smooth_ms}, as in the original method's low-pass conditioning), and
#' thresholds the smoothed trace at a multiple of baseline standard
#' deviations. The Teager-Kaiser method applies the Teager-Kaiser energy
#' operator, then the same rectify-and-smooth conditioning, before
#' thresholding. Without smoothing (\code{smooth_ms = 0}) the rectified
#' band-limited signal dips below threshold at every carrier zero-crossing, so
#' contiguous suprathreshold runs rarely reach \code{min_duration_ms} and
#' genuine bursts go undetected. The Bonato method uses a statistical
#' double-threshold on the (optionally whitened) signal, declaring activation
#' when at least \code{r} of \code{m} successive samples exceed a chi-square
#' threshold set from a target false-alarm probability. The AGLR method applies
#' Staude & Wolf's approximated generalized-likelihood-ratio change detector for
#' a variance increase and refines the onset to the estimated change point.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param method Detection method: "hodges_bui" (baseline SD threshold),
#'   "teager_kaiser" (Teager-Kaiser energy operator), "bonato" (Bonato
#'   double-threshold) or "aglr" (approximated generalized likelihood ratio).
#' @param threshold_sd Number of baseline SDs above mean for the threshold in
#'   the "hodges_bui" and "teager_kaiser" methods (default: 3).
#' @param baseline_sec Duration of the baseline period in seconds from the
#'   signal start, used to estimate baseline statistics (default: 0.2).
#' @param min_duration_ms Minimum activation duration in ms to accept
#'   (default: 50).
#' @param smooth_ms Width in ms of the centered moving-average window applied
#'   to the rectified signal ("hodges_bui") or rectified Teager-Kaiser trace
#'   ("teager_kaiser") before baseline estimation and thresholding
#'   (default: 25, in the 10--50 ms range examined by Hodges & Bui 1996).
#'   Set to 0 to disable smoothing (pre-0.2.1 behavior). Ignored by the
#'   "bonato" and "aglr" methods, which operate on unsmoothed samples by
#'   design.
#' @param r,m Bonato double-threshold parameters: activation requires at least
#'   \code{r} of \code{m} successive samples above the detection threshold
#'   (defaults: 5 of 10). For "aglr", \code{m} is the sliding analysis-window
#'   length in samples.
#' @param whiten If TRUE (default), whiten the signal with an AR model estimated
#'   from the baseline before applying the "bonato"/"aglr" statistics, so the
#'   baseline approximates white Gaussian noise.
#' @param false_alarm Target false-alarm probability that sets the detection
#'   threshold for the "bonato" (per m-sample window) and "aglr" (per sample)
#'   methods (default: 0.05).
#' @param assay_name Input assay name (default: first assay).
#' @return A list with two data.frames:
#'   \describe{
#'     \item{onsets}{A data.frame with columns \code{channel} (integer channel
#'       index), \code{sample} (sample index of onset), and \code{time_sec}
#'       (onset time in seconds). For "bonato"/"aglr" a \code{stat} column holds
#'       the peak detector statistic of the activation.}
#'     \item{offsets}{A data.frame with the matching \code{channel},
#'       \code{sample}, \code{time_sec} (and \code{stat}) for each offset.}
#'   }
#'   If no activations are detected, both data.frames have zero rows.
#' @seealso [emgEnvelope()] for computing amplitude envelopes,
#'   [emgAmplitudeNormalize()] for amplitude normalization,
#'   [emgFatigue()] for fatigue analysis
#' @references Hodges, P.W. & Bui, B.H. (1996). "A comparison of
#'   computer-based methods for the determination of onset of muscle contraction
#'   using electromyography." Electroencephalography and Clinical
#'   Neurophysiology, 101(6), 511-519. doi:10.1016/S0921-884X(96)95190-5
#' @references Bonato, P., D'Alessio, T. & Knaflitz, M. (1998). "A statistical
#'   method for the measurement of muscle activation intervals from surface
#'   myoelectric signal during gait." IEEE Transactions on Biomedical
#'   Engineering, 45(3), 287-299. doi:10.1109/10.661154
#' @references Staude, G. & Wolf, W. (1999). "Objective motor response onset
#'   detection in surface myoelectric signals." Medical Engineering & Physics,
#'   21(6-7), 449-467. doi:10.1016/S1350-4533(99)00067-3
#' @references Solnik, S., Rider, P., Steinweg, K., DeVita, P. &
#'   Hortobagyi, T. (2010). "Teager-Kaiser energy operator signal conditioning
#'   improves EMG onset detection." European Journal of Applied Physiology,
#'   110(3), 489-498. doi:10.1007/s00421-010-1521-8
#' @export
#' @examples
#' pe <- make_emg_contraction()
#' emgOnsetDetect(pe, method = "bonato")$onsets
emgOnsetDetect <- function(x, method = c("hodges_bui", "teager_kaiser",
                                         "bonato", "aglr"),
                            threshold_sd = 3, baseline_sec = 0.2,
                            min_duration_ms = 50, smooth_ms = 25,
                            r = 5L, m = 10L,
                            whiten = TRUE, false_alarm = 0.05,
                            assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  stopifnot(false_alarm > 0, false_alarm < 1, m >= 1L, r >= 1L, r <= m)
  stopifnot(is.numeric(smooth_ms), length(smooth_ms) == 1L, smooth_ms >= 0)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  baseline_samples <- max(1L, as.integer(round(baseline_sec * sr)))
  baseline_samples <- min(baseline_samples, n_time)
  baseline_idx <- seq_len(baseline_samples)
  min_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))
  smooth_win <- if (smooth_ms > 0) {
    max(1L, as.integer(round(smooth_ms / 1000 * sr)))
  } else 1L
  if (smooth_ms > 0 && smooth_win <= 1L &&
      method %in% c("hodges_bui", "teager_kaiser")) {
    warning("smooth_ms = ", smooth_ms, " ms spans < 2 samples at ", sr,
            " Hz; smoothing is effectively disabled and band-limited bursts ",
            "may go undetected", call. = FALSE)
  }
  has_stat <- method %in% c("bonato", "aglr")

  onsets_list <- list()
  offsets_list <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    stat_vec <- NULL
    refine <- NULL

    if (method == "hodges_bui") {
      rect_sig <- .emg_moving_avg(abs(sig), smooth_win)
      bl <- rect_sig[baseline_idx]
      threshold <- mean(bl, na.rm = TRUE) + threshold_sd * sd(bl, na.rm = TRUE)
      active <- rect_sig > threshold

    } else if (method == "teager_kaiser") {
      tkeo <- numeric(n_time)
      tkeo[1] <- sig[1]^2
      tkeo[n_time] <- sig[n_time]^2
      for (i in 2:(n_time - 1)) {
        tkeo[i] <- sig[i]^2 - sig[i - 1] * sig[i + 1]
      }
      tkeo <- .emg_moving_avg(abs(tkeo), smooth_win)
      bl <- tkeo[baseline_idx]
      threshold <- mean(bl, na.rm = TRUE) + threshold_sd * sd(bl, na.rm = TRUE)
      active <- tkeo > threshold

    } else if (method == "bonato") {
      det <- .bonato_detect(sig, baseline_idx, r, m, whiten, false_alarm)
      active <- det$active
      stat_vec <- det$stat

    } else if (method == "aglr") {
      det <- .aglr_detect(sig, baseline_idx, m, whiten, false_alarm)
      active <- det$active
      stat_vec <- det$stat
      refine <- det$refine
    }

    transitions <- diff(as.integer(active))
    onset_samples <- which(transitions == 1) + 1
    offset_samples <- which(transitions == -1)
    if (active[1]) onset_samples <- c(1L, onset_samples)
    if (active[n_time]) offset_samples <- c(offset_samples, n_time)

    n_events <- min(length(onset_samples), length(offset_samples))
    if (n_events > 0) {
      for (e in seq_len(n_events)) {
        on_s <- onset_samples[e]
        off_s <- offset_samples[e]
        if (!is.null(refine)) on_s <- refine(on_s)   # AGLR change-point onset
        if (off_s - on_s >= min_samples) {
          on_row <- data.frame(channel = ch, sample = on_s,
                               time_sec = (on_s - 1) / sr,
                               stringsAsFactors = FALSE)
          off_row <- data.frame(channel = ch, sample = off_s,
                                time_sec = (off_s - 1) / sr,
                                stringsAsFactors = FALSE)
          if (has_stat) {
            peak <- max(stat_vec[on_s:off_s], na.rm = TRUE)
            on_row$stat <- peak
            off_row$stat <- peak
          }
          onsets_list[[length(onsets_list) + 1]] <- on_row
          offsets_list[[length(offsets_list) + 1]] <- off_row
        }
      }
    }
  }

  empty_df <- if (has_stat) {
    data.frame(channel = integer(0), sample = integer(0),
               time_sec = numeric(0), stat = numeric(0))
  } else {
    data.frame(channel = integer(0), sample = integer(0),
               time_sec = numeric(0))
  }
  onsets_df <- if (length(onsets_list)) do.call(rbind, onsets_list) else empty_df
  offsets_df <- if (length(offsets_list)) do.call(rbind, offsets_list) else empty_df

  list(onsets = onsets_df, offsets = offsets_df)
}

# Centered moving average with edge-shrinking windows, O(n) via prefix sums
# (same edge convention as the RMS envelope in emgEnvelope()). NA samples are
# excluded from each window mean (na.rm semantics, matching the baseline
# statistics); a window with no finite sample yields NA.
.emg_moving_avg <- function(v, win) {
  if (win <= 1L) return(v)
  n <- length(v)
  half <- win %/% 2L
  ok <- !is.na(v)
  cs <- c(0, cumsum(ifelse(ok, v, 0)))
  cn <- c(0, cumsum(ok))
  idx <- seq_len(n)
  lo <- pmax(1L, idx - half)
  hi <- pmin(n, idx + half)
  cnt <- cn[hi + 1L] - cn[lo]
  out <- (cs[hi + 1L] - cs[lo]) / cnt
  out[cnt == 0L] <- NA_real_
  out
}

# Single-sample exceedance probability p0 such that the r-of-m window
# false-alarm probability P(Binom(m, p0) >= r) equals `false_alarm`.
.solve_p0 <- function(false_alarm, r, m) {
  f <- function(p) stats::pbinom(r - 1L, m, p, lower.tail = FALSE) - false_alarm
  if (f(1e-6) >= 0) return(1e-6)
  if (f(1 - 1e-6) <= 0) return(1 - 1e-6)
  stats::uniroot(f, c(1e-6, 1 - 1e-6))$root
}

# AR (prediction-error) whitening using coefficients estimated from `baseline`.
.emg_whiten <- function(sig, baseline_idx, order = 4L) {
  bl <- sig[baseline_idx]
  centered <- sig - mean(bl, na.rm = TRUE)
  fit <- tryCatch(
    stats::ar.burg(bl - mean(bl, na.rm = TRUE), aic = FALSE,
                   order.max = min(order, length(bl) - 1L), demean = FALSE),
    error = function(e) NULL)
  if (is.null(fit) || length(fit$ar) == 0L) return(centered)
  e <- as.numeric(stats::filter(centered, c(1, -fit$ar), sides = 1))
  e[is.na(e)] <- 0
  e
}

# Bonato double-threshold detector: r-of-m samples of the whitened, baseline-
# standardized chi-square statistic above a false-alarm-calibrated threshold.
.bonato_detect <- function(sig, baseline_idx, r, m, whiten, false_alarm) {
  e <- if (whiten) .emg_whiten(sig, baseline_idx) else
    sig - mean(sig[baseline_idx], na.rm = TRUE)
  sigma0 <- stats::sd(e[baseline_idx], na.rm = TRUE)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- stats::sd(e, na.rm = TRUE)
  if (!is.finite(sigma0) || sigma0 <= 0) sigma0 <- 1
  z2 <- (e / sigma0)^2                       # ~ chi-square(1) under baseline
  xi <- stats::qchisq(1 - .solve_p0(false_alarm, r, m), df = 1)
  exceed <- as.integer(z2 > xi)
  n <- length(exceed)
  cs <- c(0, cumsum(exceed))
  hi <- pmin(seq_len(n) + m - 1L, n)         # forward r-of-m window
  rc <- cs[hi + 1L] - cs[seq_len(n)]
  list(active = rc >= r, stat = z2)
}

# AGLR variance-change detector (Staude & Wolf, approximated GLR) with
# change-point onset refinement.
.aglr_detect <- function(sig, baseline_idx, m, whiten, false_alarm) {
  e <- if (whiten) .emg_whiten(sig, baseline_idx) else
    sig - mean(sig[baseline_idx], na.rm = TRUE)
  var0 <- stats::var(e[baseline_idx], na.rm = TRUE)
  if (!is.finite(var0) || var0 <= 0) var0 <- stats::var(e, na.rm = TRUE)
  if (!is.finite(var0) || var0 <= 0) var0 <- 1
  n <- length(e)
  cse2 <- c(0, cumsum(e^2))
  glr <- function(sum_sq, len) {
    th1 <- sum_sq / len
    if (th1 > var0) (len / 2) * (th1 / var0 - 1 - log(th1 / var0)) else 0
  }
  g <- numeric(n)
  for (k in seq_len(n)) {
    a <- max(1L, k - m + 1L)
    g[k] <- glr(cse2[k + 1L] - cse2[a], k - a + 1L)
  }
  h <- stats::qchisq(1 - false_alarm, df = 1) / 2
  refine <- function(k_on) {                 # argmax two-segment change point
    a <- max(1L, k_on - m + 1L)
    best <- a; best_llr <- -Inf
    for (j in (a - 1L):(k_on - 1L)) {
      llr <- glr(cse2[k_on + 1L] - cse2[j + 1L], k_on - j)
      if (llr > best_llr) { best_llr <- llr; best <- j + 1L }
    }
    best
  }
  list(active = g > h, stat = g, refine = refine)
}
