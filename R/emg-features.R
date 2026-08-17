#' EMG Amplitude Features (ARV, MAV, iEMG, RMS)
#'
#' Computes time-domain amplitude features of an EMG signal over a sliding
#' window. The signal is divided into overlapping windows and, for each window
#' and channel, the requested features are computed:
#' \describe{
#'   \item{\code{arv}}{Average rectified value, \eqn{\frac{1}{N}\sum |x_i|}.}
#'   \item{\code{mav}}{Mean absolute value, \eqn{\frac{1}{N}\sum |x_i|}
#'     (identical to ARV for an unweighted window).}
#'   \item{\code{iemg}}{Integrated EMG, \eqn{\sum |x_i|\,\Delta t} with
#'     \eqn{\Delta t = 1 / f_s}; the time integral of the rectified signal,
#'     in signal-units times seconds.}
#'   \item{\code{rms}}{Root mean square, \eqn{\sqrt{\frac{1}{N}\sum x_i^2}}.}
#' }
#' The window geometry (window and step sizes in samples) matches
#' [emgFatigue()] and [emgSpectralMoments()] for identical \code{window_sec}
#' and \code{overlap}, so amplitude and spectral features align window-for-window.
#'
#' @param x A PhysioExperiment object with EMG data.
#' @param features Character vector selecting which features to return; any of
#'   \code{"arv"}, \code{"mav"}, \code{"iemg"}, \code{"rms"} (default: all four).
#' @param window_sec Analysis window in seconds (default: 1.0).
#' @param overlap Overlap fraction between consecutive windows, in
#'   \eqn{[0, 1)} (default: 0.5).
#' @param assay_name Input assay name (default: first assay).
#' @return A data.frame with one row per analysis window per channel, with
#'   columns:
#'   \describe{
#'     \item{channel}{Integer channel index.}
#'     \item{window}{Integer window index within the channel (1-based).}
#'     \item{time_sec}{Window start time in seconds.}
#'     \item{arv, mav, iemg, rms}{The requested amplitude features (only the
#'       columns named in \code{features} are present).}
#'   }
#'   If no full window fits in the signal, a 0-row data.frame with these columns
#'   is returned.
#' @seealso [emgEnvelope()] for a per-sample amplitude envelope,
#'   [emgFatigue()] for windowed spectral fatigue features,
#'   [emgSpectralMoments()] for spectral moments on the same window grid
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @references Hermens, H.J. et al. (2000). "Development of recommendations for
#'   SEMG sensors and sensor placement procedures (SENIAM)." Journal of
#'   Electromyography and Kinesiology, 10(5), 361-374.
#'   doi:10.1016/S1050-6411(00)00027-4
#' @export
#' @examples
#' # 2 s of two-channel EMG at 1000 Hz
#' set.seed(1)
#' m <- matrix(rnorm(2000 * 2, sd = 0.3), nrow = 2000, ncol = 2)
#' pe <- PhysioExperiment(
#'   assays = list(raw = m),
#'   colData = S4Vectors::DataFrame(label = c("EMG1", "EMG2"),
#'                                  type = c("EMG", "EMG")),
#'   samplingRate = 1000)
#' feats <- emgAmplitudeFeatures(pe, window_sec = 0.5, overlap = 0.5)
#' head(feats)
emgAmplitudeFeatures <- function(x, features = c("arv", "mav", "iemg", "rms"),
                                 window_sec = 1.0, overlap = 0.5,
                                 assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  features <- match.arg(features, c("arv", "mav", "iemg", "rms"),
                        several.ok = TRUE)
  # canonical order, independent of how the caller ordered `features`
  features <- intersect(c("arv", "mav", "iemg", "rms"), features)
  stopifnot(is.numeric(window_sec), length(window_sec) == 1L, window_sec > 0)
  stopifnot(is.numeric(overlap), length(overlap) == 1L,
            overlap >= 0, overlap < 1)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)
  dt <- 1 / sr

  # Window geometry identical to emgFatigue()/emgSpectralMoments().
  win_samples <- as.integer(round(window_sec * sr))
  step_samples <- as.integer(round(win_samples * (1 - overlap)))

  empty <- function() {
    out <- data.frame(channel = integer(0), window = integer(0),
                      time_sec = numeric(0), stringsAsFactors = FALSE)
    for (f in features) out[[f]] <- numeric(0)
    out
  }
  if (win_samples < 1L || win_samples > n_time) return(empty())

  results <- list()
  for (ch in seq_len(n_channels)) {
    sig <- data[, ch]
    win_idx <- 0L
    start <- 1L
    while (start + win_samples - 1 <= n_time) {
      win_idx <- win_idx + 1L
      end <- start + win_samples - 1
      segment <- sig[start:end]
      abs_seg <- abs(segment)

      row <- data.frame(channel = ch, window = win_idx,
                        time_sec = (start - 1) / sr,
                        stringsAsFactors = FALSE)
      if ("arv" %in% features)  row$arv  <- mean(abs_seg, na.rm = TRUE)
      if ("mav" %in% features)  row$mav  <- mean(abs_seg, na.rm = TRUE)
      if ("iemg" %in% features) row$iemg <- sum(abs_seg, na.rm = TRUE) * dt
      if ("rms" %in% features)  row$rms  <- sqrt(mean(segment^2, na.rm = TRUE))

      results[[length(results) + 1]] <- row
      start <- start + step_samples
    }
  }

  if (!length(results)) return(empty())
  do.call(rbind, results)
}
