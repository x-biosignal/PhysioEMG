#' Create test EMG PhysioExperiment
#' @param n_time Number of time points
#' @param n_channels Number of EMG channels
#' @param sr Sampling rate in Hz
#' @return PhysioExperiment with simulated EMG data
make_emg <- function(n_time = 2000, n_channels = 4, sr = 1000) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  data <- matrix(NA_real_, nrow = n_time, ncol = n_channels)
  for (ch in seq_len(n_channels)) {
    # Simulated EMG: white noise burst (active) + low noise (rest)
    signal <- rnorm(n_time, sd = 0.01)  # baseline noise
    # Active burst from 30-70% of signal
    burst_start <- as.integer(n_time * 0.3)
    burst_end <- as.integer(n_time * 0.7)
    signal[burst_start:burst_end] <- rnorm(burst_end - burst_start + 1, sd = 0.5 + 0.1 * ch)
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

#' Create test EMG data with known fatigue pattern
#' @param n_time Number of time points
#' @param sr Sampling rate
#' @return PhysioExperiment with frequency-shifting EMG
make_emg_fatigue <- function(n_time = 10000, sr = 1000) {
  t <- seq(0, (n_time - 1) / sr, length.out = n_time)

  # Simulated fatiguing contraction: median frequency decreases over time
  n_segments <- 10
  seg_len <- n_time %/% n_segments
  signal <- numeric(n_time)

  for (i in seq_len(n_segments)) {
    start <- (i - 1) * seg_len + 1
    end <- min(i * seg_len, n_time)
    center_freq <- 80 - 4 * (i - 1)  # decreasing from 80 to 44 Hz
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

#' Create multi-channel EMG for synergy analysis
#' @param n_time Number of time points
#' @param n_channels Number of muscles
#' @param n_synergies True number of underlying synergies
#' @param sr Sampling rate
#' @return PhysioExperiment with synergy-structured EMG
make_emg_synergy <- function(n_time = 1000, n_channels = 8, n_synergies = 3, sr = 1000) {
  # Generate synergy activation patterns (time x synergies)
  H <- matrix(0, nrow = n_time, ncol = n_synergies)
  for (s in seq_len(n_synergies)) {
    center <- n_time * s / (n_synergies + 1)
    H[, s] <- exp(-((seq_len(n_time) - center)^2) / (2 * (n_time / 8)^2))
  }

  # Synergy weights (synergies x channels)
  W <- matrix(abs(rnorm(n_synergies * n_channels, sd = 0.5)),
              nrow = n_synergies, ncol = n_channels)

  # Data = H %*% W + noise
  data <- H %*% W + matrix(rnorm(n_time * n_channels, sd = 0.05),
                            nrow = n_time, ncol = n_channels)
  data[data < 0] <- 0  # EMG is non-negative after rectification

  PhysioExperiment(
    assays = list(raw = data),
    colData = S4Vectors::DataFrame(
      label = paste0("Muscle", seq_len(n_channels)),
      type = rep("EMG", n_channels)
    ),
    samplingRate = sr
  )
}
