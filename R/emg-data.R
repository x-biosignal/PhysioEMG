#' Create a Basic EMG PhysioExperiment
#'
#' Generates a synthetic multi-channel EMG \code{PhysioExperiment} object with
#' simulated muscle activity. Each channel contains baseline noise with a burst
#' of higher-amplitude activity in the middle 40 percent of the signal
#' (from 30 to 70 percent), mimicking a typical voluntary contraction.
#'
#' @param n_time Number of time points (default: 2000).
#' @param n_channels Number of EMG channels (default: 4).
#' @param sr Sampling rate in Hz (default: 1000).
#' @return A \code{PhysioExperiment} object with a single \code{"raw"} assay
#'   containing simulated EMG data (time x channels matrix), channel metadata
#'   in \code{colData}, and the specified sampling rate.
#' @seealso [make_emg_contraction()] for EMG with a defined contraction window,
#'   [make_emg_fatigue()] for EMG with fatigue progression,
#'   [emgEnvelope()] for extracting amplitude envelopes
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
#' @examples
#' pe <- make_emg()
#' pe
#' dim(SummarizedExperiment::assay(pe, "raw"))
make_emg <- function(n_time = 2000, n_channels = 4, sr = 1000) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  for (ch in seq_len(n_channels)) {
    signal <- rnorm(n_time, sd = 0.01)
    burst_start <- as.integer(n_time * 0.3)
    burst_end <- as.integer(n_time * 0.7)
    signal[burst_start:burst_end] <- rnorm(burst_end - burst_start + 1,
                                           sd = 0.5 + 0.1 * ch)
    data[, ch] <- signal
  }

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("EMG", seq_len(n_channels)),
      type = rep("EMG", n_channels)
    ),
    samplingRate = sr
  )
}

#' Create EMG PhysioExperiment with Known Contraction
#'
#' Generates a synthetic single-channel EMG \code{PhysioExperiment} with a
#' clearly defined contraction period. The contraction region has substantially
#' higher amplitude than the baseline, making it suitable for testing onset
#' detection algorithms.
#'
#' @param n_time Number of time points (default: 5000).
#' @param sr Sampling rate in Hz (default: 1000).
#' @param contraction_start Proportion of signal where contraction begins
#'   (default: 0.3).
#' @param contraction_end Proportion of signal where contraction ends
#'   (default: 0.7).
#' @param baseline_sd Standard deviation of baseline noise (default: 0.01).
#' @param contraction_sd Standard deviation of contraction activity
#'   (default: 0.5).
#' @return A \code{PhysioExperiment} object with a single \code{"raw"} assay
#'   containing a one-channel EMG signal with a known contraction window. The
#'   contraction region is defined by \code{contraction_start} and
#'   \code{contraction_end} as proportions of total signal length.
#' @seealso [emgOnsetDetect()] for detecting the contraction onset,
#'   [emgEnvelope()] for extracting the amplitude envelope,
#'   [make_emg()] for basic multi-channel EMG,
#'   [make_emg_fatigue()] for EMG with fatigue progression
#' @references Hodges, P.W. & Bui, B.H. (1996). "A comparison of
#'   computer-based methods for the determination of onset of muscle contraction
#'   using electromyography." Electroencephalography and Clinical
#'   Neurophysiology, 101(6), 511-519. doi:10.1016/S0921-884X(96)95190-5
#' @export
#' @examples
#' pe <- make_emg_contraction()
#' pe
#' onset <- emgOnsetDetect(pe)
#' onset$onsets
make_emg_contraction <- function(n_time = 5000, sr = 1000,
                                  contraction_start = 0.3,
                                  contraction_end = 0.7,
                                  baseline_sd = 0.01,
                                  contraction_sd = 0.5) {
  signal <- rnorm(n_time, sd = baseline_sd)
  start_idx <- as.integer(n_time * contraction_start)
  end_idx <- as.integer(n_time * contraction_end)
  signal[start_idx:end_idx] <- rnorm(end_idx - start_idx + 1, sd = contraction_sd)
  data <- matrix(signal, ncol = 1)

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "EMG1", type = "EMG"),
    samplingRate = sr
  )
}

#' Create EMG PhysioExperiment with Fatigue Progression
#'
#' Generates a synthetic single-channel EMG \code{PhysioExperiment} that
#' simulates a fatiguing isometric contraction. The median frequency of the
#' signal decreases progressively over time (from 80 Hz to approximately
#' 44 Hz across 10 segments), mimicking the spectral compression
#' characteristic of muscle fatigue.
#'
#' @param n_time Number of time points (default: 10000).
#' @param sr Sampling rate in Hz (default: 1000).
#' @return A \code{PhysioExperiment} object with a single \code{"raw"} assay
#'   containing a one-channel EMG signal exhibiting progressive median
#'   frequency decrease across 10 equal-length segments.
#' @seealso [emgFatigue()] for tracking median frequency over time,
#'   [emgFatigueIndex()] for computing a summary fatigue metric,
#'   [emgSpectralMoments()] for spectral moment analysis,
#'   [make_emg()] for basic multi-channel EMG
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
#' @examples
#' pe <- make_emg_fatigue()
#' pe
#' fatigue <- emgFatigue(pe)
#' head(fatigue)
make_emg_fatigue <- function(n_time = 10000, sr = 1000) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  n_segments <- 10
  seg_len <- n_time %/% n_segments
  signal <- numeric(n_time)

  for (i in seq_len(n_segments)) {
    start <- (i - 1) * seg_len + 1
    end <- min(i * seg_len, n_time)
    center_freq <- 80 - 4 * (i - 1)
    idx <- start:end
    signal[idx] <- sin(2 * pi * center_freq * t[idx]) * 0.3 +
                   rnorm(length(idx), sd = 0.2)
  }

  data <- matrix(signal, ncol = 1)

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(label = "EMG1", type = "EMG"),
    samplingRate = sr
  )
}
