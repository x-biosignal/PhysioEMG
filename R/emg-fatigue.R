#' EMG Fatigue Analysis
#'
#' Tracks median and mean frequency over time to assess muscle fatigue.
#' Decreasing median frequency indicates fatigue due to reduced motor unit
#' conduction velocity. The signal is divided into overlapping windows and
#' the power spectral density is computed via FFT for each window.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param window_sec Analysis window in seconds (default: 1.0).
#' @param overlap Overlap fraction between windows (default: 0.5).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per channel per window, containing columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{window}{Integer window number (1-indexed).}
#'     \item{time_sec}{Start time of the window in seconds.}
#'     \item{median_freq}{Median frequency (Hz) at which 50 percent of the
#'       spectral power is below.}
#'     \item{mean_freq}{Power-weighted mean frequency (Hz).}
#'     \item{rms_amplitude}{Root mean square amplitude of the window.}
#'   }
#' @seealso [emgFatigueIndex()] for a summary fatigue metric,
#'   [emgSpectralMoments()] for spectral moment analysis,
#'   [emgEnvelope()] for amplitude envelope extraction
#' @references De Luca, C.J. (1984). "Myoelectrical manifestations of localized
#'   muscular fatigue in humans." Critical Reviews in Biomedical Engineering,
#'   11(4), 251-279.
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
emgFatigue <- function(x, window_sec = 1.0, overlap = 0.5,
                        assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  win_samples <- as.integer(round(window_sec * sr))
  step_samples <- as.integer(round(win_samples * (1 - overlap)))

  results <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    win_idx <- 0L

    start <- 1L
    while (start + win_samples - 1 <= n_time) {
      win_idx <- win_idx + 1L
      end <- start + win_samples - 1

      segment <- sig[start:end]

      # Power spectrum
      n <- length(segment)
      psd <- (Mod(fft(segment))^2) / n
      freqs <- seq(0, sr / 2, length.out = n %/% 2 + 1)
      psd_half <- psd[seq_len(length(freqs))]

      # Remove DC component
      psd_half[1] <- 0

      total_power <- sum(psd_half)

      if (total_power > 0) {
        # Median frequency: frequency at which 50% of power is below
        cum_power <- cumsum(psd_half)
        mdf_idx <- which(cum_power >= total_power / 2)[1]
        mdf <- freqs[mdf_idx]

        # Mean frequency: power-weighted average
        mnf <- sum(freqs * psd_half) / total_power
      } else {
        mdf <- NA_real_
        mnf <- NA_real_
      }

      rms <- sqrt(mean(segment^2, na.rm = TRUE))

      results[[length(results) + 1]] <- data.frame(
        channel = ch, window = win_idx,
        time_sec = (start - 1) / sr,
        median_freq = mdf, mean_freq = mnf,
        rms_amplitude = rms,
        stringsAsFactors = FALSE)

      start <- start + step_samples
    }
  }

  do.call(rbind, results)
}

#' EMG Fatigue Index
#'
#' Computes a fatigue index as the ratio of final to initial median frequency.
#' Values less than 1 indicate fatigue (frequency decrease). The signal is
#' internally analyzed with [emgFatigue()] using 0.5-second windows at 50
#' percent overlap, then the initial and final portions are compared.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param initial_pct Percentage of signal used for initial estimate (default: 0.2).
#' @param final_pct Percentage of signal used for final estimate (default: 0.2).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per channel, containing columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{fatigue_index}{Ratio of final to initial median frequency. Values
#'       less than 1 indicate fatigue.}
#'     \item{initial_mdf}{Mean median frequency (Hz) in the initial portion.}
#'     \item{final_mdf}{Mean median frequency (Hz) in the final portion.}
#'   }
#' @seealso [emgFatigue()] for detailed windowed fatigue tracking,
#'   [emgSpectralMoments()] for spectral moment analysis,
#'   [emgEnvelope()] for amplitude envelope extraction
#' @references De Luca, C.J. (1984). "Myoelectrical manifestations of localized
#'   muscular fatigue in humans." Critical Reviews in Biomedical Engineering,
#'   11(4), 251-279.
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
emgFatigueIndex <- function(x, initial_pct = 0.2, final_pct = 0.2,
                             assay_name = NULL) {
  fatigue_data <- emgFatigue(x, window_sec = 0.5, overlap = 0.5,
                              assay_name = assay_name)

  channels <- unique(fatigue_data$channel)
  results <- list()

  for (ch in channels) {
    ch_data <- fatigue_data[fatigue_data$channel == ch, ]
    n_windows <- nrow(ch_data)

    n_initial <- max(1L, as.integer(round(n_windows * initial_pct)))
    n_final <- max(1L, as.integer(round(n_windows * final_pct)))

    initial_mdf <- mean(ch_data$median_freq[seq_len(n_initial)], na.rm = TRUE)
    final_mdf <- mean(ch_data$median_freq[(n_windows - n_final + 1):n_windows], na.rm = TRUE)

    results[[length(results) + 1]] <- data.frame(
      channel = ch,
      fatigue_index = final_mdf / initial_mdf,
      initial_mdf = initial_mdf,
      final_mdf = final_mdf,
      stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}

#' EMG Spectral Moments
#'
#' Computes spectral moments (M0, M1, M2) over sliding windows. M0 is total
#' power, M1 is the first spectral moment (related to mean frequency), and M2
#' is the second moment (related to bandwidth). These can be combined to derive
#' the mean frequency (M1/M0) and spectral bandwidth.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param window_sec Analysis window in seconds (default: 1.0).
#' @param overlap Overlap fraction (default: 0.5).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per channel per window, containing columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{window}{Integer window number (1-indexed).}
#'     \item{m0}{Zeroth spectral moment (total power).}
#'     \item{m1}{First spectral moment (frequency-weighted power).}
#'     \item{m2}{Second spectral moment (frequency-squared-weighted power).}
#'   }
#' @seealso [emgFatigue()] for median and mean frequency tracking,
#'   [emgFatigueIndex()] for summary fatigue metric,
#'   [emgEnvelope()] for time-domain amplitude analysis
#' @references De Luca, C.J. (1984). "Myoelectrical manifestations of localized
#'   muscular fatigue in humans." Critical Reviews in Biomedical Engineering,
#'   11(4), 251-279.
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
emgSpectralMoments <- function(x, window_sec = 1.0, overlap = 0.5,
                                assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  win_samples <- as.integer(round(window_sec * sr))
  step_samples <- as.integer(round(win_samples * (1 - overlap)))

  results <- list()

  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    win_idx <- 0L

    start <- 1L
    while (start + win_samples - 1 <= n_time) {
      win_idx <- win_idx + 1L
      end <- start + win_samples - 1

      segment <- sig[start:end]
      n <- length(segment)
      psd <- (Mod(fft(segment))^2) / n
      freqs <- seq(0, sr / 2, length.out = n %/% 2 + 1)
      psd_half <- psd[seq_len(length(freqs))]
      psd_half[1] <- 0  # Remove DC

      df <- if (length(freqs) > 1) freqs[2] - freqs[1] else 1
      m0 <- sum(psd_half) * df
      m1 <- sum(freqs * psd_half) * df
      m2 <- sum(freqs^2 * psd_half) * df

      results[[length(results) + 1]] <- data.frame(
        channel = ch, window = win_idx,
        m0 = m0, m1 = m1, m2 = m2,
        stringsAsFactors = FALSE)

      start <- start + step_samples
    }
  }

  do.call(rbind, results)
}
