# ---- internal spectral / signal helpers ------------------------------------

# One-sided power spectral density (matches the package's FFT convention).
.emg_psd <- function(sig, sr) {
  n <- length(sig)
  psd <- (Mod(stats::fft(sig - mean(sig)))^2) / n
  freqs <- seq(0, sr / 2, length.out = n %/% 2 + 1)
  list(freqs = freqs, psd = psd[seq_len(length(freqs))])
}

# Power in a frequency band [flo, fhi] from a one-sided PSD.
.emg_band_power <- function(ps, flo, fhi) {
  sel <- ps$freqs >= flo & ps$freqs <= fhi
  if (!any(sel)) return(0)
  sum(ps$psd[sel])
}

# Symmetric (real-output) FFT band-pass / high-pass via a folded-frequency mask.
.emg_fft_mask <- function(sig, sr, mask_fun) {
  n <- length(sig)
  ft <- stats::fft(sig)
  freqs <- (0:(n - 1)) * sr / n
  f <- pmin(freqs, sr - freqs)                 # fold to [0, sr/2]
  Re(stats::fft(ft * mask_fun(f), inverse = TRUE)) / n
}

.emg_fft_highpass <- function(sig, sr, cutoff, trans = 5) {
  lo <- cutoff - trans
  .emg_fft_mask(sig, sr, function(f)
    ifelse(f >= cutoff, 1,
           ifelse(f <= lo, 0, 0.5 * (1 - cos(pi * (f - lo) / (cutoff - lo))))))
}

.emg_fft_bandpass <- function(sig, sr, flo, fhi) {
  .emg_fft_mask(sig, sr, function(f) as.numeric(f >= flo & f <= fhi))
}

.emg_moving_rms <- function(sig, win) {
  cs <- c(0, cumsum(sig^2))
  n <- length(sig)
  hi <- pmin(seq_len(n) + win - 1L, n)
  sqrt((cs[hi + 1L] - cs[seq_len(n)]) / (hi - seq_len(n) + 1L))
}

# ECG-contamination score: prominence of the spectral peak of the QRS-band
# envelope at heart-rate frequencies (0.6-3.5 Hz, ~36-210 bpm) relative to the
# surrounding envelope spectrum, squashed to [0, 1). Periodic QRS produces a
# sharp, prominent peak; sustained voluntary contractions and broadband noise do
# not, so this separates ECG contamination from ordinary EMG activity.
.emg_ecg_score <- function(sig, sr) {
  env <- abs(.emg_fft_bandpass(sig, sr, 8, 40))
  env <- env - mean(env)
  ps <- .emg_psd(env, sr)
  hr <- ps$freqs >= 0.6 & ps$freqs <= 3.5
  nb <- ps$freqs > 0.1 & ps$freqs <= 8
  if (!any(hr) || !any(nb)) return(0)
  base <- stats::median(ps$psd[nb])
  if (!is.finite(base) || base <= 0) return(0)
  prom <- max(ps$psd[hr]) / base
  prom / (prom + 30)                    # ~0.5 at prom = 30
}

# Simple R-peak detection on the QRS-band energy of a (possibly EMG) signal.
.emg_detect_rpeaks <- function(sig, sr) {
  energy <- .emg_fft_bandpass(sig, sr, 8, 40)^2
  win <- max(1L, as.integer(round(0.02 * sr)))
  sm <- .emg_moving_rms(sqrt(energy), win)^2
  thr <- mean(sm) + 2 * stats::sd(sm)
  min_rr <- max(1L, as.integer(round(sr * 60 / 220)))   # refractory: 220 bpm
  n <- length(sm)
  peaks <- integer(0)
  i <- 1L
  while (i <= n) {
    if (sm[i] > thr) {
      j <- i
      while (j <= n && sm[j] > thr) j <- j + 1L
      seg <- i:(j - 1L)
      peaks <- c(peaks, seg[which.max(sm[seg])])
      i <- j + min_rr
    } else {
      i <- i + 1L
    }
  }
  peaks
}

#' EMG Signal Quality Check
#'
#' Computes per-channel surface-EMG quality metrics: signal-to-noise ratio,
#' resting baseline noise, power-line interference ratio at 50 and 60 Hz,
#' clipping/saturation percentage, and an ECG-contamination flag based on the
#' periodicity of the QRS-band envelope. When the PhysioECG package is available
#' the R-peak count it detects is reported as an additional diagnostic.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param ecg_score_threshold QRS-band periodicity score above which a channel
#'   is flagged as ECG-contaminated (default: 0.5).
#' @param rms_window_ms Window (ms) for the moving-RMS used to estimate signal
#'   and rest levels (default: 100).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per channel:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{snr_db}{Signal-to-noise ratio in dB (active vs rest RMS).}
#'     \item{baseline_noise}{Resting baseline RMS (noise floor).}
#'     \item{powerline_50, powerline_60}{Fraction of total power within +/-1 Hz
#'       of 50 and 60 Hz.}
#'     \item{powerline_ratio}{The larger of \code{powerline_50}/\code{powerline_60}.}
#'     \item{clipping_pct}{Percentage of samples at the signal extremes
#'       (ADC saturation).}
#'     \item{ecg_score}{QRS-band periodicity score from 0 (none) to 1 (strong).}
#'     \item{ecg_contamination}{TRUE if \code{ecg_score} exceeds the threshold.}
#'   }
#' @seealso [emgRemoveECG()] to remove ECG contamination,
#'   [emgQCgate()] for pass/fail gating, [emgEnvelope()]
#' @references McManus, L., De Vito, G. & Lowery, M.M. (2020). "Analysis and
#'   biophysics of surface EMG for physiotherapists and kinesiologists."
#'   Frontiers in Neurology, 11, 576729. doi:10.3389/fneur.2020.576729
#' @export
#' @examples
#' set.seed(1)
#' sig <- rnorm(4000, sd = 0.02)
#' sig[1500:2500] <- rnorm(1001, sd = 0.4)          # a contraction
#' pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
#'                        samplingRate = 1000)
#' emgQualityCheck(pe)
emgQualityCheck <- function(x, ecg_score_threshold = 0.5,
                            rms_window_ms = 100, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)
  win <- max(1L, as.integer(round(rms_window_ms / 1000 * sr)))

  have_ecg <- requireNamespace("PhysioECG", quietly = TRUE)

  rows <- lapply(seq_len(n_channels), function(ch) {
    sig <- data[, ch]
    mrms <- .emg_moving_rms(sig, win)
    noise <- stats::quantile(mrms, 0.1, names = FALSE, na.rm = TRUE)
    signal_lvl <- stats::quantile(mrms, 0.9, names = FALSE, na.rm = TRUE)
    snr_db <- if (is.finite(noise) && noise > 0) {
      20 * log10(signal_lvl / noise)
    } else {
      Inf
    }

    ps <- .emg_psd(sig, sr)
    total <- sum(ps$psd)
    pl50 <- if (total > 0) .emg_band_power(ps, 49, 51) / total else 0
    pl60 <- if (total > 0) .emg_band_power(ps, 59, 61) / total else 0

    rail <- max(abs(sig), na.rm = TRUE)
    clip_pct <- if (rail > 0) {
      100 * mean(abs(sig) >= rail * (1 - 1e-6), na.rm = TRUE)
    } else {
      0
    }

    score <- .emg_ecg_score(sig, sr)

    data.frame(
      channel = ch, snr_db = snr_db, baseline_noise = noise,
      powerline_50 = pl50, powerline_60 = pl60,
      powerline_ratio = max(pl50, pl60), clipping_pct = clip_pct,
      ecg_score = score, ecg_contamination = score > ecg_score_threshold,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  # optional PhysioECG diagnostic: number of R-peaks it detects per channel
  if (have_ecg) {
    out$ecg_rpeaks <- vapply(seq_len(n_channels), function(ch) {
      pk <- tryCatch(PhysioECG::ecgDetectRpeaks(
        PhysioCore::PhysioExperiment(assays = list(raw = data[, ch, drop = FALSE]),
                                     samplingRate = sr)),
        error = function(e) NULL)
      if (is.null(pk)) NA_integer_ else
        tryCatch(length(pk$peaks[[1]]), error = function(e) NA_integer_)
    }, integer(1))
  }
  out
}

#' Remove ECG Contamination from EMG
#'
#' Removes cardiac (ECG/QRS) contamination from surface EMG using a high-pass
#' filter, QRS template subtraction, or gating (blanking) of the QRS complexes.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param method Removal method: "highpass" (default; FFT high-pass at
#'   \code{hp_cutoff}), "template" (average-QRS template subtraction), or
#'   "gating" (linear-interpolation blanking of QRS windows).
#' @param ecg_channel Optional channel index of a dedicated ECG reference used
#'   for R-peak detection in the "template"/"gating" methods; if NULL, R-peaks
#'   are detected from each EMG channel itself.
#' @param hp_cutoff High-pass cutoff in Hz for the "highpass" method
#'   (default: 30).
#' @param qrs_ms Half-width in ms of the QRS window for "template"/"gating"
#'   (default: 60).
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: "ecg_removed").
#' @return A PhysioExperiment object with an added assay \code{output_assay}
#'   containing the cleaned signal.
#' @seealso [emgQualityCheck()] to detect contamination, [emgQCgate()]
#' @references Drake, J.D.M. & Callaghan, J.P. (2006). "Elimination of
#'   electrocardiogram contamination from electromyogram signals: An evaluation
#'   of currently used removal techniques." Journal of Electromyography and
#'   Kinesiology, 16(2), 175-187. doi:10.1016/j.jelekin.2005.07.003
#' @references Willigenburg, N.W., Daffertshofer, A., Kingma, I. & van Dieen,
#'   J.H. (2012). "Removing ECG contamination from EMG recordings." Journal of
#'   Electromyography and Kinesiology, 22(3), 485-493.
#'   doi:10.1016/j.jelekin.2012.01.001
#' @export
#' @examples
#' set.seed(1)
#' sr <- 1000; n <- 5000
#' emg <- rnorm(n, sd = 0.1)
#' qrs <- numeric(n)
#' for (p in seq(500, n - 200, by = 1000)) qrs[p + -15:15] <- 3 * exp(-(-15:15)^2 / 20)
#' pe <- PhysioExperiment(assays = list(raw = matrix(emg + qrs, ncol = 1)),
#'                        samplingRate = sr)
#' clean <- emgRemoveECG(pe, method = "highpass")
emgRemoveECG <- function(x, method = c("highpass", "template", "gating"),
                         ecg_channel = NULL, hp_cutoff = 30, qrs_ms = 60,
                         assay_name = NULL, output_assay = "ecg_removed") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)
  half <- max(1L, as.integer(round(qrs_ms / 1000 * sr)))

  cleaned <- data
  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]

    if (method == "highpass") {
      cleaned[, ch] <- .emg_fft_highpass(sig, sr, hp_cutoff)
      next
    }

    det <- if (!is.null(ecg_channel)) data[, ecg_channel] else sig
    peaks <- .emg_detect_rpeaks(det, sr)
    peaks <- peaks[peaks - half >= 1L & peaks + half <= n_time]
    if (length(peaks) < 3L) {           # not enough beats: leave as-is
      cleaned[, ch] <- sig
      next
    }

    if (method == "template") {
      windows <- vapply(peaks, function(p) sig[(p - half):(p + half)],
                        numeric(2L * half + 1L))
      template <- rowMeans(windows)
      out <- sig
      for (p in peaks) out[(p - half):(p + half)] <-
        sig[(p - half):(p + half)] - template
      cleaned[, ch] <- out

    } else {                            # gating: linear-interpolate the window
      out <- sig
      for (p in peaks) {
        idx <- (p - half):(p + half)
        a <- sig[max(1L, p - half - 1L)]
        b <- sig[min(n_time, p + half + 1L)]
        out[idx] <- seq(a, b, length.out = length(idx))
      }
      cleaned[, ch] <- out
    }
  }

  dimnames(cleaned) <- dimnames(data)
  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- cleaned
  SummarizedExperiment::assays(x) <- assays
  x
}

#' EMG Quality-Control Gate
#'
#' Applies pass/fail thresholds to the [emgQualityCheck()] metrics and reports,
#' per channel, whether the signal is acceptable and, if not, why.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param min_snr Minimum acceptable SNR in dB (default: 6).
#' @param max_powerline Maximum acceptable power-line ratio (default: 0.10).
#' @param max_clipping Maximum acceptable clipping percentage (default: 1).
#' @param reject_ecg If TRUE (default), ECG-contaminated channels fail the gate.
#' @param assay_name Input assay name (default: first assay).
#' @param ... Passed to [emgQualityCheck()].
#' @return A list with \code{pass} (TRUE if every channel passes) and
#'   \code{channels}, a data.frame with columns \code{channel}, \code{pass} and
#'   \code{reasons} (a comma-separated string of failed criteria, or "").
#' @seealso [emgQualityCheck()], [emgRemoveECG()]
#' @references McManus, L., De Vito, G. & Lowery, M.M. (2020). "Analysis and
#'   biophysics of surface EMG for physiotherapists and kinesiologists."
#'   Frontiers in Neurology, 11, 576729. doi:10.3389/fneur.2020.576729
#' @export
#' @examples
#' set.seed(1)
#' sig <- rnorm(4000, sd = 0.02); sig[1500:2500] <- rnorm(1001, sd = 0.4)
#' pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
#'                        samplingRate = 1000)
#' emgQCgate(pe)$pass
emgQCgate <- function(x, min_snr = 6, max_powerline = 0.10, max_clipping = 1,
                      reject_ecg = TRUE, assay_name = NULL, ...) {
  qc <- emgQualityCheck(x, assay_name = assay_name, ...)

  rows <- lapply(seq_len(nrow(qc)), function(i) {
    reasons <- character(0)
    if (!is.na(qc$snr_db[i]) && qc$snr_db[i] < min_snr) {
      reasons <- c(reasons, sprintf("low SNR (%.1f < %.1f dB)",
                                    qc$snr_db[i], min_snr))
    }
    if (qc$powerline_ratio[i] > max_powerline) {
      reasons <- c(reasons, sprintf("power-line (%.3f > %.3f)",
                                    qc$powerline_ratio[i], max_powerline))
    }
    if (qc$clipping_pct[i] > max_clipping) {
      reasons <- c(reasons, sprintf("clipping (%.2f%% > %.2f%%)",
                                    qc$clipping_pct[i], max_clipping))
    }
    if (reject_ecg && isTRUE(qc$ecg_contamination[i])) {
      reasons <- c(reasons, "ECG contamination")
    }
    data.frame(channel = qc$channel[i], pass = length(reasons) == 0,
               reasons = paste(reasons, collapse = ", "),
               stringsAsFactors = FALSE)
  })
  channels <- do.call(rbind, rows)
  list(pass = all(channels$pass), channels = channels)
}
