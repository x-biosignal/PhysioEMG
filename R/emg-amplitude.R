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
        envelope[i, ch] <- sqrt(sum(sig_sq[lo:hi]) / (hi - lo + 1))
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

  x
}

#' Normalize EMG Amplitude
#'
#' Normalizes EMG amplitude data using MVC or peak normalization. MVC
#' normalization divides each channel by its maximum voluntary contraction
#' value, yielding percentage-of-MVC units. Peak normalization divides by the
#' within-trial peak so each channel ranges from 0 to 1.
#'
#' @param x A PhysioExperiment object.
#' @param method Normalization method: "mvc" (maximum voluntary contraction)
#'   or "peak" (peak of the trial).
#' @param mvc_data A PhysioExperiment containing MVC trial data (required for
#'   "mvc" method). Must have the same number of channels.
#' @param assay_name Assay to normalize (default: first assay).
#' @param output_assay Output assay name (default: "normalized").
#' @return A PhysioExperiment object with an additional assay named
#'   \code{output_assay} containing normalized amplitude values. For "peak"
#'   normalization, values range from 0 to 1. For "mvc" normalization, values
#'   represent proportion of maximum voluntary contraction.
#' @seealso [emgEnvelope()] for computing amplitude envelopes prior to
#'   normalization, [emgFatigue()] for fatigue analysis,
#'   [emgOnsetDetect()] for muscle activation onset detection
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
emgAmplitudeNormalize <- function(x, method = c("mvc", "peak"),
                                   mvc_data = NULL, assay_name = NULL,
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
  }

  dimnames(normalized) <- dimnames(data)
  assays <- SummarizedExperiment::assays(x)
  assays[[output_assay]] <- normalized
  SummarizedExperiment::assays(x) <- assays

  x
}
