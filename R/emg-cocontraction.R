# Resolve an agonist/antagonist argument (channel index or label) to an index.
.resolve_cc_channel <- function(x, ch, n_channels) {
  if (is.numeric(ch)) {
    idx <- as.integer(ch[1])
  } else {
    labels <- tryCatch(as.character(SummarizedExperiment::colData(x)$label),
                       error = function(e) NULL)
    idx <- match(as.character(ch[1]), labels)
    if (is.na(idx)) stop("channel '", ch[1], "' not found", call. = FALSE)
  }
  if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > n_channels) {
    stop("channel out of range: ", ch[1], call. = FALSE)
  }
  idx
}

# Falconer & Winter (1985) percent co-contraction: 2 * common area / total area.
.cci_falconer_winter <- function(a, b) {
  tot <- sum(a + b, na.rm = TRUE)
  if (tot <= 0) return(0)
  2 * sum(pmin(a, b), na.rm = TRUE) / tot * 100
}

# Rudolph et al. (2000): mean of (lower/higher) * (lower + higher).
.cci_rudolph <- function(a, b) {
  lower <- pmin(a, b); higher <- pmax(a, b)
  val <- ifelse(higher > 0, (lower / higher) * (lower + higher), 0)
  mean(val, na.rm = TRUE)
}

# Frost et al. (1997) variant: time-average of the instantaneous co-activation
# ratio 2 * min(a, b) / (a + b), as a percentage.
.cci_frost <- function(a, b) {
  s <- a + b
  val <- ifelse(s > 0, 2 * pmin(a, b) / s, 0)
  mean(val, na.rm = TRUE) * 100
}

#' EMG Co-Contraction Index
#'
#' Quantifies simultaneous activation of an agonist-antagonist muscle pair from
#' their amplitude envelopes, using one of three established indices. Input
#' channels are treated as (non-negative) amplitude envelopes; with
#' \code{normalize = TRUE} each is first scaled to its own peak so the two
#' muscles are amplitude-comparable.
#'
#' \describe{
#'   \item{\code{"falconer_winter"}}{Falconer & Winter (1985) percent
#'     co-contraction: \eqn{2\sum\min(a,b) / \sum(a+b) \times 100}. 100 when the
#'     two envelopes are identical, 0 when they never overlap.}
#'   \item{\code{"rudolph"}}{Rudolph et al. (2000) co-contraction index: the mean
#'     over time of \eqn{(\ell/h)(\ell+h)} with \eqn{\ell=\min(a,b)},
#'     \eqn{h=\max(a,b)}. Symmetric in the two muscles and non-negative.}
#'   \item{\code{"frost"}}{Frost et al. (1997) variant: the time-average of the
#'     instantaneous co-activation ratio \eqn{2\min(a,b)/(a+b)}, as a percentage.}
#' }
#'
#' @param x A PhysioExperiment object of EMG amplitude envelopes.
#' @param agonist,antagonist Channel index (integer) or channel label
#'   (character) of the two muscles. The indices are interchangeable for all
#'   three methods (the indices are symmetric in the pair).
#' @param method Co-contraction index: "falconer_winter" (default), "rudolph",
#'   or "frost".
#' @param normalize If TRUE (default), peak-normalize each channel to \[0, 1\]
#'   before computing the index.
#' @param window_sec Optional non-overlapping window length in seconds for a
#'   time-resolved index. If NULL (default) only the whole-signal summary is
#'   returned.
#' @param assay_name Input assay name (default: first assay).
#' @return A list with:
#'   \describe{
#'     \item{method}{The method used.}
#'     \item{summary}{The whole-signal co-contraction index (a scalar).}
#'     \item{timeseries}{If \code{window_sec} is given, a data.frame with columns
#'       \code{window}, \code{time_sec} and \code{cci} (one row per window);
#'       otherwise \code{NULL}.}
#'   }
#' @seealso [emgEnvelope()] for computing amplitude envelopes,
#'   [emgAmplitudeNormalize()] for reference normalization,
#'   [emgOnsetDetect()] for activation timing
#' @references Falconer, K. & Winter, D.A. (1985). "Quantitative assessment of
#'   co-contraction at the ankle joint in walking." Electromyography and Clinical
#'   Neurophysiology, 25(2-3), 135-149.
#' @references Rudolph, K.S., Axe, M.J. & Snyder-Mackler, L. (2000). "Dynamic
#'   stability after ACL injury: who can hop?" Knee Surgery, Sports Traumatology,
#'   Arthroscopy, 8(5), 262-269. doi:10.1007/s001670000130
#' @references Frost, G., Dowling, J., Dyson, K. & Bar-Or, O. (1997).
#'   "Cocontraction in three age groups of children during treadmill
#'   locomotion." Journal of Electromyography and Kinesiology, 7(3), 179-186.
#'   doi:10.1016/S1050-6411(97)84626-3
#' @export
#' @examples
#' set.seed(1)
#' t <- seq_len(1000)
#' a <- exp(-((t - 400) / 120)^2)            # agonist envelope
#' b <- exp(-((t - 600) / 120)^2)            # antagonist envelope
#' pe <- PhysioExperiment(
#'   assays = list(env = cbind(a, b)),
#'   colData = S4Vectors::DataFrame(label = c("TA", "GAS"),
#'                                  type = c("EMG", "EMG")),
#'   samplingRate = 1000)
#' emgCoContraction(pe, agonist = "TA", antagonist = "GAS")$summary
emgCoContraction <- function(x, agonist, antagonist,
                             method = c("falconer_winter", "rudolph", "frost"),
                             normalize = TRUE, window_sec = NULL,
                             assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  sr <- samplingRate(x)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  ai <- .resolve_cc_channel(x, agonist, n_channels)
  bi <- .resolve_cc_channel(x, antagonist, n_channels)

  a <- abs(data[, ai])                        # envelopes are non-negative
  b <- abs(data[, bi])
  if (normalize) {
    ma <- max(a, na.rm = TRUE); mb <- max(b, na.rm = TRUE)
    if (is.finite(ma) && ma > 0) a <- a / ma
    if (is.finite(mb) && mb > 0) b <- b / mb
  }

  cci_fun <- switch(method,
    falconer_winter = .cci_falconer_winter,
    rudolph = .cci_rudolph,
    frost = .cci_frost)

  summary_val <- cci_fun(a, b)

  timeseries <- NULL
  if (!is.null(window_sec)) {
    stopifnot(is.numeric(window_sec), length(window_sec) == 1L, window_sec > 0)
    win <- max(1L, as.integer(round(window_sec * sr)))
    starts <- seq(1L, n_time - win + 1L, by = win)
    if (length(starts) > 0) {
      rows <- lapply(seq_along(starts), function(k) {
        idx <- starts[k]:(starts[k] + win - 1L)
        data.frame(window = k, time_sec = (starts[k] - 1) / sr,
                   cci = cci_fun(a[idx], b[idx]), stringsAsFactors = FALSE)
      })
      timeseries <- do.call(rbind, rows)
    }
  }

  list(method = method, summary = summary_val, timeseries = timeseries)
}
