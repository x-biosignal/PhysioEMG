#' Detect EMG Muscle Activation Onset and Offset
#'
#' Identifies when muscle activation begins and ends using threshold-based
#' or energy-based methods. The Hodges-Bui method thresholds the rectified
#' signal at a multiple of baseline standard deviations. The Teager-Kaiser
#' method applies the Teager-Kaiser energy operator before thresholding,
#' which can be more sensitive to sudden onsets.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param method Detection method: "hodges_bui" (baseline SD threshold) or
#'   "teager_kaiser" (Teager-Kaiser energy operator).
#' @param threshold_sd Number of baseline SDs above mean for threshold
#'   (default: 3).
#' @param baseline_sec Duration of baseline period in seconds from signal
#'   start (default: 0.2).
#' @param min_duration_ms Minimum activation duration in ms to accept
#'   (default: 50).
#' @param assay_name Input assay name (default: first assay).
#' @return A list with two data.frames:
#'   \describe{
#'     \item{onsets}{A data.frame with columns \code{channel} (integer channel
#'       index), \code{sample} (sample index of onset), and \code{time_sec}
#'       (onset time in seconds).}
#'     \item{offsets}{A data.frame with columns \code{channel}, \code{sample},
#'       and \code{time_sec} for each activation offset. Rows correspond to
#'       matching onsets.}
#'   }
#'   If no activations are detected, both data.frames have zero rows.
#' @seealso [emgEnvelope()] for computing amplitude envelopes,
#'   [emgAmplitudeNormalize()] for amplitude normalization,
#'   [emgFatigue()] for fatigue analysis
#' @references Hodges, P.W. & Bui, B.H. (1996). "A comparison of
#'   computer-based methods for the determination of onset of muscle contraction
#'   using electromyography." Electroencephalography and Clinical
#'   Neurophysiology, 101(6), 511-519. doi:10.1016/S0921-884X(96)95190-5
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
emgOnsetDetect <- function(x, method = c("hodges_bui", "teager_kaiser"),
                            threshold_sd = 3, baseline_sec = 0.2,
                            min_duration_ms = 50,
                            assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  baseline_samples <- max(1L, as.integer(round(baseline_sec * sr)))
  min_samples <- max(1L, as.integer(round(min_duration_ms / 1000 * sr)))

  onsets_list <- list()
  offsets_list <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    if (method == "hodges_bui") {
      rect_sig <- abs(sig)
      bl <- rect_sig[seq_len(baseline_samples)]
      bl_mean <- mean(bl, na.rm = TRUE)
      bl_sd <- sd(bl, na.rm = TRUE)
      threshold <- bl_mean + threshold_sd * bl_sd
      active <- rect_sig > threshold

    } else if (method == "teager_kaiser") {
      tkeo <- numeric(n_time)
      tkeo[1] <- sig[1]^2
      tkeo[n_time] <- sig[n_time]^2
      for (i in 2:(n_time - 1)) {
        tkeo[i] <- sig[i]^2 - sig[i - 1] * sig[i + 1]
      }
      tkeo <- abs(tkeo)

      bl <- tkeo[seq_len(baseline_samples)]
      bl_mean <- mean(bl, na.rm = TRUE)
      bl_sd <- sd(bl, na.rm = TRUE)
      threshold <- bl_mean + threshold_sd * bl_sd
      active <- tkeo > threshold
    }

    transitions <- diff(as.integer(active))
    onset_samples <- which(transitions == 1) + 1
    offset_samples <- which(transitions == -1)

    if (active[1]) onset_samples <- c(1L, onset_samples)
    if (active[n_time]) offset_samples <- c(offset_samples, n_time)

    n_events <- min(length(onset_samples), length(offset_samples))
    if (n_events > 0) {
      for (e in seq_len(n_events)) {
        duration <- offset_samples[e] - onset_samples[e]
        if (duration >= min_samples) {
          onsets_list[[length(onsets_list) + 1]] <- data.frame(
            channel = ch, sample = onset_samples[e],
            time_sec = (onset_samples[e] - 1) / sr, stringsAsFactors = FALSE)
          offsets_list[[length(offsets_list) + 1]] <- data.frame(
            channel = ch, sample = offset_samples[e],
            time_sec = (offset_samples[e] - 1) / sr, stringsAsFactors = FALSE)
        }
      }
    }
  }

  onsets_df <- if (length(onsets_list) > 0) do.call(rbind, onsets_list) else
    data.frame(channel = integer(0), sample = integer(0), time_sec = numeric(0))
  offsets_df <- if (length(offsets_list) > 0) do.call(rbind, offsets_list) else
    data.frame(channel = integer(0), sample = integer(0), time_sec = numeric(0))

  list(onsets = onsets_df, offsets = offsets_df)
}
