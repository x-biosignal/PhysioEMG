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

# Simple linear regression of a fatigue metric against time (internal helper).
# Coefficients, R^2 and the two-sided slope p-value are computed in closed form
# (no summary.lm, so an essentially perfect fit does not raise a warning).
.fatigue_regress <- function(value, time_sec) {
  keep <- is.finite(value) & is.finite(time_sec)
  value <- value[keep]; time_sec <- time_sec[keep]
  n <- length(value)
  if (n < 3L || stats::sd(time_sec) == 0) {
    return(list(slope = NA_real_, intercept = NA_real_,
                r2 = NA_real_, p = NA_real_))
  }
  xbar <- mean(time_sec); ybar <- mean(value)
  Sxx <- sum((time_sec - xbar)^2)
  Sxy <- sum((time_sec - xbar) * (value - ybar))
  Syy <- sum((value - ybar)^2)
  slope <- Sxy / Sxx
  intercept <- ybar - slope * xbar
  sse <- sum((value - (intercept + slope * time_sec))^2)
  r2 <- if (Syy > 0) 1 - sse / Syy else 1
  df <- n - 2L
  p <- if (sse <= .Machine$double.eps * Syy) {
    0                       # essentially perfect fit -> slope highly significant
  } else {
    se_slope <- sqrt((sse / df) / Sxx)
    2 * stats::pt(-abs(slope / se_slope), df = df)
  }
  list(slope = slope, intercept = intercept, r2 = r2, p = p)
}

#' EMG Fatigue Regression Slope
#'
#' Fits a linear regression of a spectral fatigue metric (median or mean
#' frequency) against time across sliding windows and reports the fatigue slope
#' per channel. A negative slope indicates myoelectric fatigue (spectral
#' compression toward lower frequencies).
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param feature Metric to regress: "mdf" (median frequency) or "mnf" (mean
#'   frequency).
#' @param normalize If TRUE, also report the slope as a percentage of the
#'   initial (intercept) value per minute.
#' @param window_sec Analysis window in seconds (default: 1.0).
#' @param overlap Overlap fraction between windows (default: 0.5).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per channel, containing columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{slope_hz_per_min}{Regression slope in Hz per minute (negative under
#'       fatigue).}
#'     \item{norm_slope_pct_per_min}{Slope as a percentage of the initial value
#'       per minute (\code{NA} if \code{normalize = FALSE}).}
#'     \item{intercept_hz}{Fitted value at time zero (initial frequency), Hz.}
#'     \item{r_squared}{Coefficient of determination of the fit.}
#'     \item{p_value}{Two-sided p-value for the slope.}
#'   }
#' @seealso [emgFatigue()] for the underlying windowed frequency estimates,
#'   [emgDimitrovIndex()] for the spectral-moment fatigue index,
#'   [emgFatigueIndex()] for the initial/final ratio metric
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
#' @examples
#' # decreasing-frequency (fatiguing) signal
#' sr <- 1000; n <- 8000; t <- seq_len(n) / sr
#' f <- 90 - 3 * t                      # 90 Hz falling to ~66 Hz
#' sig <- sin(2 * pi * cumsum(f) / sr) + rnorm(n, sd = 0.1)
#' pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
#'                        samplingRate = sr)
#' emgFatigueSlope(pe, feature = "mdf")
emgFatigueSlope <- function(x, feature = c("mdf", "mnf"), normalize = TRUE,
                            window_sec = 1.0, overlap = 0.5,
                            assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  feature <- match.arg(feature)
  col <- if (feature == "mdf") "median_freq" else "mean_freq"

  fat <- emgFatigue(x, window_sec = window_sec, overlap = overlap,
                    assay_name = assay_name)
  chans <- sort(unique(fat$channel))

  rows <- lapply(chans, function(ch) {
    d <- fat[fat$channel == ch, , drop = FALSE]
    r <- .fatigue_regress(d[[col]], d$time_sec)
    slope_min <- r$slope * 60          # Hz per second -> Hz per minute
    norm_slope <- if (normalize && is.finite(r$intercept) &&
                      r$intercept != 0) {
      100 * slope_min / r$intercept
    } else {
      NA_real_
    }
    data.frame(channel = ch, slope_hz_per_min = slope_min,
               norm_slope_pct_per_min = norm_slope,
               intercept_hz = r$intercept, r_squared = r$r2,
               p_value = r$p, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' EMG Dimitrov Spectral Fatigue Index (FInsm5)
#'
#' Computes spectral moments \eqn{M(k) = \int f^k\,PSD(f)\,df} for
#' \eqn{k = -1, 0, \dots, 5} over sliding windows and the Dimitrov fatigue index
#' \eqn{FInsm5 = M(-1)/M(5)} per window. FInsm5 rises steeply with fatigue
#' because spectral compression toward low frequencies simultaneously increases
#' the low-frequency moment \eqn{M(-1)} and decreases the high-frequency moment
#' \eqn{M(5)}. The per-channel FInsm5-vs-time regression slope is attached as the
#' \code{"finsm5_slope"} attribute.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param window_sec Analysis window in seconds (default: 1.0).
#' @param overlap Overlap fraction between windows (default: 0.5).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per channel per window, containing columns
#'   \code{channel}, \code{window}, \code{time_sec}, the moments
#'   \code{m_minus1}, \code{m0}, \code{m1}, \code{m2}, \code{m3}, \code{m4},
#'   \code{m5}, and \code{finsm5}. The per-channel FInsm5 regression slope
#'   (index per minute) is available via \code{attr(result, "finsm5_slope")}.
#' @seealso [emgFatigueSlope()] for MDF/MNF regression slopes,
#'   [emgSpectralMoments()] for the M0/M1/M2 moments,
#'   [emgFatigue()] for median/mean frequency tracking
#' @references Dimitrov, G.V. et al. (2006). "Muscle fatigue during dynamic
#'   contractions assessed by new spectral indices." Medicine & Science in
#'   Sports & Exercise, 38(11), 1971-1979. doi:10.1249/01.mss.0000233794.31659.6d
#' @export
#' @examples
#' sr <- 1000; n <- 6000
#' sig <- sin(2 * pi * 80 * seq_len(n) / sr) + rnorm(n, sd = 0.1)
#' pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
#'                        samplingRate = sr)
#' di <- emgDimitrovIndex(pe, window_sec = 1)
#' head(di)
#' attr(di, "finsm5_slope")
emgDimitrovIndex <- function(x, window_sec = 1.0, overlap = 0.5,
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
      psd_half[1] <- 0  # remove DC
      df <- if (length(freqs) > 1) freqs[2] - freqs[1] else 1

      # positive-frequency mask so the k = -1 moment never divides by f = 0
      pos <- freqs > 0
      fp <- freqs[pos]; pp <- psd_half[pos]
      mom <- function(k) sum(fp^k * pp) * df
      m_neg1 <- mom(-1)
      m0 <- mom(0); m1 <- mom(1); m2 <- mom(2)
      m3 <- mom(3); m4 <- mom(4); m5 <- mom(5)
      finsm5 <- if (m5 > 0) m_neg1 / m5 else NA_real_

      results[[length(results) + 1]] <- data.frame(
        channel = ch, window = win_idx, time_sec = (start - 1) / sr,
        m_minus1 = m_neg1, m0 = m0, m1 = m1, m2 = m2,
        m3 = m3, m4 = m4, m5 = m5, finsm5 = finsm5,
        stringsAsFactors = FALSE)

      start <- start + step_samples
    }
  }

  out <- do.call(rbind, results)

  # per-channel FInsm5-vs-time slope (index per minute)
  chans <- sort(unique(out$channel))
  slope_rows <- lapply(chans, function(ch) {
    d <- out[out$channel == ch, , drop = FALSE]
    r <- .fatigue_regress(d$finsm5, d$time_sec)
    data.frame(channel = ch, slope_per_min = r$slope * 60,
               r_squared = r$r2, p_value = r$p, stringsAsFactors = FALSE)
  })
  attr(out, "finsm5_slope") <- do.call(rbind, slope_rows)
  out
}
