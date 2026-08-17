#' EMG Envelope Extraction
#'
#' Extracts the amplitude envelope from EMG signals using various methods.
#' The RMS method computes root mean square over a sliding window. The Hilbert
#' method uses the analytic signal via the Hilbert transform. The lowpass method
#' rectifies the signal and applies a moving-average lowpass filter.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param method Envelope method: "rms" (root mean square), "hilbert"
#'   (Hilbert transform), or "lowpass" (rectification + lowpass filter).
#' @param window_ms Window size in milliseconds for RMS method (default: 50).
#' @param cutoff Cutoff frequency in Hz for lowpass method (default: 6).
#' @param assay_name Input assay name (default: first assay).
#' @param output_assay Output assay name (default: "envelope").
#' @return A PhysioExperiment object with an additional assay named
#'   \code{output_assay} containing the amplitude envelope. The envelope
#'   matrix has the same dimensions as the input (time x channels) with
#'   non-negative values representing instantaneous signal amplitude.
#' @seealso [emgAmplitudeNormalize()] for normalizing envelope values,
#'   [emgOnsetDetect()] for onset detection from envelope data,
#'   [emgFatigue()] for fatigue analysis using spectral features
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
emgEnvelope <- function(x, method = c("rms", "hilbert", "lowpass"),
                        window_ms = 50, cutoff = 6,
                        assay_name = NULL, output_assay = "envelope") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  envelope <- matrix(NA_real_, nrow = n_time, ncol = n_channels)

  if (method == "rms") {
    window_samples <- max(1L, as.integer(round(window_ms / 1000 * sr)))
    half_win <- window_samples %/% 2

    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      sig_sq <- sig^2
      cs <- c(0, cumsum(sig_sq))
      for (i in seq_len(n_time)) {
        lo <- max(1L, i - half_win)
        hi <- min(n_time, i + half_win)
        # O(1) windowed sum from the prefix-sum `cs` (cs[k+1] = sum(sig_sq[1:k]))
        envelope[i, ch] <- sqrt((cs[hi + 1] - cs[lo]) / (hi - lo + 1))
      }
    }

  } else if (method == "hilbert") {
    for (ch in seq_len(n_channels)) {
      sig <- data[, ch]
      n <- length(sig)
      ft <- fft(sig)
      h <- rep(0, n)
      if (n > 0) {
        h[1] <- 1
        if (n %% 2 == 0) {
          h[n / 2 + 1] <- 1
          h[2:(n / 2)] <- 2
        } else {
          h[2:((n + 1) / 2)] <- 2
        }
      }
      analytic <- fft(ft * h, inverse = TRUE) / n
      envelope[, ch] <- Mod(analytic)
    }

  } else if (method == "lowpass") {
    rectified <- abs(data)
    nyquist <- sr / 2
    if (cutoff >= nyquist) {
      warning("Cutoff frequency >= Nyquist, returning rectified signal", call. = FALSE)
      envelope <- rectified
    } else {
      window_samples <- max(1L, as.integer(round(sr / cutoff / 2)))
      for (ch in seq_len(n_channels)) {
        sig <- rectified[, ch]
        kernel <- rep(1 / window_samples, window_samples)
        filtered <- stats::filter(sig, kernel, sides = 2)
        envelope[, ch] <- as.numeric(filtered)
        na_idx <- which(is.na(envelope[, ch]))
        if (length(na_idx) > 0) {
          envelope[na_idx, ch] <- rectified[na_idx, ch]
        }
      }
    }
  }

  dimnames(envelope) <- dimnames(data)
  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- envelope
  SummarizedExperiment::assays(x) <- assays

  .recordProv(x, input_assay = assay_name, output_assay = output_assay,
              .package = "PhysioEMG")
}

#' Normalize EMG Amplitude
#'
#' Normalizes EMG amplitude data by a per-channel reference value.
#' \code{"mvc"} divides by the per-channel maximum of a maximum-voluntary-
#' contraction trial (percentage-of-MVC). \code{"peak"} divides by the
#' within-trial peak so each channel ranges from 0 to 1. \code{"rvc"} divides by
#' the per-channel mean amplitude of a sub-maximal reference voluntary
#' contraction trial (percentage-of-RVC). \code{"dynamic_peak"} divides by a
#' centered moving maximum, tracking a time-varying peak within the trial.
#'
#' Input is assumed to be amplitude data (a rectified signal or an
#' [emgEnvelope()] output). MVC and RVC normalization are scale-invariant to a
#' gain applied to both the signal and its reference trial.
#'
#' @param x A PhysioExperiment object (amplitude data).
#' @param method Normalization method: "mvc" (maximum voluntary contraction),
#'   "peak" (within-trial peak), "rvc" (sub-maximal reference voluntary
#'   contraction), or "dynamic_peak" (centered moving maximum).
#' @param mvc_data A PhysioExperiment containing MVC trial data (required for
#'   "mvc"). Must have the same number of channels as \code{x}.
#' @param rvc_data A PhysioExperiment containing the sub-maximal reference-task
#'   data (required for "rvc"). Must have the same number of channels as
#'   \code{x}.
#' @param rvc_window Optional window. For "rvc", a length-2 numeric
#'   \code{c(start, end)} in seconds selecting the portion of \code{rvc_data}
#'   over which the mean reference amplitude is computed (default: the whole
#'   reference trial). For "dynamic_peak", the moving-maximum window length in
#'   seconds (default: 0.5).
#' @param assay_name Assay to normalize (default: first assay).
#' @param output_assay Output assay name (default: "normalized").
#' @return A PhysioExperiment object with an additional assay named
#'   \code{output_assay} containing normalized amplitude values. For "peak" and
#'   "dynamic_peak", values are scaled to a peak. For "mvc"/"rvc", values are
#'   proportions of the reference (multiply by 100 for percentage-of-MVC /
#'   percentage-of-RVC).
#' @seealso [emgEnvelope()] for computing amplitude envelopes prior to
#'   normalization, [emgAmplitudeFeatures()] for windowed amplitude features,
#'   [emgFatigue()] for fatigue analysis,
#'   [emgOnsetDetect()] for muscle activation onset detection
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @references Yang, J.F. & Winter, D.A. (1984). "Electromyographic amplitude
#'   normalization methods: improving their sensitivity as diagnostic tools in
#'   gait analysis." Archives of Physical Medicine and Rehabilitation, 65(9),
#'   517-521.
#' @references Burden, A. & Bartlett, R. (1999). "Normalisation of EMG amplitude:
#'   an evaluation and comparison of old and new methods." Medical Engineering &
#'   Physics, 21(4), 247-257. doi:10.1016/S1350-4533(99)00054-5
#' @export
#' @examples
#' set.seed(1)
#' amp <- matrix(abs(rnorm(1000 * 2, sd = 0.3)), nrow = 1000, ncol = 2)
#' pe <- PhysioExperiment(
#'   assays = list(raw = amp),
#'   colData = S4Vectors::DataFrame(label = c("EMG1", "EMG2"),
#'                                  type = c("EMG", "EMG")),
#'   samplingRate = 1000)
#' # percentage-of-RVC using a sub-maximal reference (here, the trial itself)
#' pe_rvc <- emgAmplitudeNormalize(pe, method = "rvc", rvc_data = pe)
#' colMeans(SummarizedExperiment::assay(pe_rvc, "normalized"))  # ~ 1
emgAmplitudeNormalize <- function(x, method = c("mvc", "peak", "rvc",
                                                "dynamic_peak"),
                                   mvc_data = NULL, rvc_data = NULL,
                                   rvc_window = NULL, assay_name = NULL,
                                   output_assay = "normalized") {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_channels <- ncol(data)

  if (method == "mvc") {
    if (is.null(mvc_data)) {
      stop("mvc_data is required for MVC normalization", call. = FALSE)
    }
    stopifnot(inherits(mvc_data, "PhysioExperiment"))
    mvc_assay <- SummarizedExperiment::assay(mvc_data, assay_name)

    if (ncol(mvc_assay) != n_channels) {
      stop("mvc_data must have the same number of channels", call. = FALSE)
    }

    mvc_max <- apply(mvc_assay, 2, max, na.rm = TRUE)
    normalized <- sweep(data, 2, mvc_max, "/")

  } else if (method == "peak") {
    peak_vals <- apply(data, 2, max, na.rm = TRUE)
    normalized <- sweep(data, 2, peak_vals, "/")

  } else if (method == "rvc") {
    if (is.null(rvc_data)) {
      stop("rvc_data is required for RVC normalization", call. = FALSE)
    }
    stopifnot(inherits(rvc_data, "PhysioExperiment"))
    rvc_assay <- SummarizedExperiment::assay(rvc_data, assay_name)

    if (ncol(rvc_assay) != n_channels) {
      stop("rvc_data must have the same number of channels", call. = FALSE)
    }

    # Optionally restrict the reference to a steady sub-maximal window.
    if (!is.null(rvc_window)) {
      stopifnot(is.numeric(rvc_window), length(rvc_window) == 2L,
                rvc_window[1] < rvc_window[2])
      sr_ref <- samplingRate(rvc_data)
      lo <- max(1L, as.integer(floor(rvc_window[1] * sr_ref)) + 1L)
      hi <- min(nrow(rvc_assay), as.integer(ceiling(rvc_window[2] * sr_ref)))
      rvc_assay <- rvc_assay[lo:hi, , drop = FALSE]
    }

    # Mean amplitude of the reference task per channel (Yang & Winter 1984).
    rvc_ref <- colMeans(abs(rvc_assay), na.rm = TRUE)
    normalized <- sweep(data, 2, rvc_ref, "/")

  } else if (method == "dynamic_peak") {
    win_sec <- if (is.null(rvc_window)) 0.5 else rvc_window[1]
    stopifnot(is.numeric(win_sec), length(win_sec) >= 1L, win_sec[1] > 0)
    sr <- samplingRate(x)
    win_samp <- max(1L, as.integer(round(win_sec[1] * sr)))
    normalized <- data
    for (ch in seq_len(n_channels)) {
      mmax <- .rolling_max_centered(abs(data[, ch]), win_samp)
      mmax[mmax == 0] <- NA_real_
      normalized[, ch] <- data[, ch] / mmax
    }
  }

  dimnames(normalized) <- dimnames(data)
  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- normalized
  SummarizedExperiment::assays(x) <- assays

  x
}

# Centered moving maximum over a window of `w` samples. Internal helper for
# emgAmplitudeNormalize(method = "dynamic_peak").
.rolling_max_centered <- function(v, w) {
  n <- length(v)
  if (n == 0L) return(v)
  half <- w %/% 2L
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- if (i - half < 1L) 1L else i - half
    hi <- if (i + half > n) n else i + half
    out[i] <- max(v[lo:hi])
  }
  out
}
