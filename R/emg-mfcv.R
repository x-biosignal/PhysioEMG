# Parabolic sub-sample refinement of a discrete cross-correlation peak.
.parabolic_peak <- function(y_left, y_center, y_right) {
  denom <- (y_left - 2 * y_center + y_right)
  if (denom == 0) return(0)
  0.5 * (y_left - y_right) / denom
}

# Simple 1-D phase unwrap (add multiples of 2*pi to keep steps < pi).
.unwrap_phase <- function(p) {
  if (length(p) < 2L) return(p)
  d <- diff(p)
  d <- d - 2 * pi * round(d / (2 * pi))
  c(p[1], p[1] + cumsum(d))
}

# Delay (in samples) of `b` behind `a` by cross-correlation. Positive = b lags a.
.mfcv_delay_xcorr <- function(a, b, max_lag) {
  a <- a - mean(a); b <- b - mean(b)
  cc <- stats::ccf(a, b, lag.max = max_lag, plot = FALSE,
                   type = "correlation")
  acf <- as.numeric(cc$acf)
  lag <- as.numeric(cc$lag)
  i <- which.max(acf)
  peak_lag <- lag[i]           # ccf peaks at k = -delay(b behind a)
  if (i > 1L && i < length(acf)) {
    peak_lag <- peak_lag + .parabolic_peak(acf[i - 1L], acf[i], acf[i + 1L])
  }
  list(delay_samples = -peak_lag, peak_corr = acf[i])
}

# Delay (in samples) of `b` behind `a` from the cross-spectrum phase slope.
.mfcv_delay_phase <- function(a, b, sr, fmin, fmax) {
  a <- a - mean(a); b <- b - mean(b)
  n <- length(a)
  cross <- fft(a) * Conj(fft(b))          # cross-spectrum a * conj(b)
  half <- seq_len(n %/% 2)
  f <- (half - 1) * sr / n
  cs <- cross[half]
  band <- f >= fmin & f <= fmax & f > 0
  if (sum(band) < 3L) return(list(delay_samples = NA_real_, peak_corr = NA_real_))
  ph <- .unwrap_phase(Arg(cs[band]))
  w <- Mod(cs[band])                       # weight by cross-spectral magnitude
  fit <- stats::lm(ph ~ f[band], weights = w)
  slope <- unname(stats::coef(fit)[2])     # dphase/df ; phase = -2*pi*f*tau
  tau <- -slope / (2 * pi)                 # seconds; b behind a if tau > 0
  list(delay_samples = tau * sr,
       peak_corr = suppressWarnings(summary(fit)$r.squared))
}

#' Muscle-Fiber Conduction Velocity (MFCV)
#'
#' Estimates muscle-fiber conduction velocity from the propagation delay of the
#' EMG signal between pairs of electrodes aligned along the fiber direction.
#' The delay between the two channels of each pair is estimated either by
#' cross-correlation (with parabolic sub-sample interpolation) or from the slope
#' of the cross-spectrum phase, and the velocity is the inter-electrode distance
#' divided by the delay.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param electrode_pairs Channel pairs aligned along the muscle fibres, given
#'   as a list of length-2 integer vectors \code{list(c(prox, dist), ...)} or as
#'   a two-column matrix (one pair per row; column 1 proximal, column 2 distal).
#' @param ied_mm Inter-electrode distance in millimetres. A single value applied
#'   to all pairs, or one value per pair.
#' @param method Delay-estimation method: "xcorr" (cross-correlation, default)
#'   or "phase" (cross-spectrum phase slope).
#' @param band Frequency band \code{c(fmin, fmax)} in Hz used by the "phase"
#'   method (default \code{c(20, 250)}).
#' @param max_lag_ms Maximum delay searched by the "xcorr" method, in
#'   milliseconds (default: 25).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per electrode pair, containing columns:
#'   \describe{
#'     \item{pair}{Integer pair index.}
#'     \item{ch1, ch2}{Proximal and distal channel indices.}
#'     \item{delay_ms}{Estimated propagation delay in milliseconds.}
#'     \item{velocity_m_s}{Conduction velocity in metres per second
#'       (\code{ied_mm} / delay).}
#'     \item{quality}{Peak cross-correlation ("xcorr") or phase-fit R^2
#'       ("phase").}
#'   }
#' @seealso [emgFatigueSlope()] and [emgDimitrovIndex()] for spectral fatigue
#'   metrics, [emgFatigue()] for median/mean frequency tracking
#' @references Farina, D. & Merletti, R. (2000). "Comparison of algorithms for
#'   estimation of EMG variables during voluntary isometric contractions."
#'   Journal of Electromyography and Kinesiology, 10(5), 337-349.
#'   doi:10.1016/S1050-6411(00)00025-0
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
#' @examples
#' # two electrodes 10 mm apart; distal is the proximal signal delayed 3 samples
#' sr <- 2000
#' set.seed(1)
#' prox <- as.numeric(stats::filter(rnorm(4000), rep(1, 5), sides = 2))
#' prox[is.na(prox)] <- 0
#' dist <- c(rep(0, 3), prox[seq_len(length(prox) - 3)])
#' pe <- PhysioExperiment(assays = list(raw = cbind(prox, dist)),
#'                        samplingRate = sr)
#' emgMFCV(pe, electrode_pairs = list(c(1, 2)), ied_mm = 10)
emgMFCV <- function(x, electrode_pairs, ied_mm, method = c("xcorr", "phase"),
                    band = c(20, 250), max_lag_ms = 25, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  # normalise electrode_pairs to a list of length-2 integer vectors
  if (is.matrix(electrode_pairs)) {
    if (ncol(electrode_pairs) != 2L) {
      stop("electrode_pairs matrix must have 2 columns", call. = FALSE)
    }
    electrode_pairs <- lapply(seq_len(nrow(electrode_pairs)),
                              function(i) as.integer(electrode_pairs[i, ]))
  }
  stopifnot(is.list(electrode_pairs), length(electrode_pairs) >= 1L)
  if (!all(vapply(electrode_pairs, length, integer(1)) == 2L)) {
    stop("each electrode pair must have exactly 2 channels", call. = FALSE)
  }

  n_pairs <- length(electrode_pairs)
  if (length(ied_mm) == 1L) ied_mm <- rep(ied_mm, n_pairs)
  stopifnot(length(ied_mm) == n_pairs, all(ied_mm > 0))

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_channels <- ncol(data)
  max_lag <- max(1L, as.integer(round(max_lag_ms / 1000 * sr)))

  rows <- lapply(seq_len(n_pairs), function(k) {
    pr <- electrode_pairs[[k]]
    if (any(pr < 1L) || any(pr > n_channels)) {
      stop("electrode_pairs reference channels outside 1..", n_channels,
           call. = FALSE)
    }
    a <- data[, pr[1]]; b <- data[, pr[2]]
    est <- if (method == "xcorr") {
      .mfcv_delay_xcorr(a, b, max_lag)
    } else {
      .mfcv_delay_phase(a, b, sr, band[1], band[2])
    }
    delay_s <- abs(est$delay_samples) / sr
    velocity <- if (is.finite(delay_s) && delay_s > 0) {
      (ied_mm[k] / 1000) / delay_s
    } else {
      NA_real_
    }
    data.frame(pair = k, ch1 = pr[1], ch2 = pr[2],
               delay_ms = est$delay_samples / sr * 1000,
               velocity_m_s = velocity, quality = est$peak_corr,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
