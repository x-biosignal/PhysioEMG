# EMG assessment of an ADL task: how the muscles are used to perform it.
#
# A clinical scale records THAT a reach-to-drink succeeded; EMG records HOW --
# the agonist drive (as %MVC), and crucially the antagonist co-contraction that
# marks guarded, inefficient or spastic movement (high after stroke, in
# dystonia, in the learning phase of a task). This summarises the muscle
# activation of one ADL task from an EMG envelope: the agonist/antagonist
# co-contraction index (reusing emgCoContraction), each muscle's peak and mean
# amplitude (as %MVC when a maximal-effort reference is supplied), and how much
# of the task each muscle is active for. It attaches the ICF code the task
# realises so the result joins the cross-modal ICF construct.

.emg_resolve_ch <- function(x, ch) {
  if (is.numeric(ch)) return(as.integer(ch))
  labs <- as.character(SummarizedExperiment::colData(x)$label)
  i <- match(as.character(ch), labs)
  if (is.na(i)) stop("channel '", ch, "' not found in colData$label.",
                     call. = FALSE)
  i
}

.emg_adl_icf <- function(task) {
  c(reaching = "d445", drinking = "d560", feeding = "d550",
    dressing = "d540", grooming = "d520")[[task]]
}

#' EMG muscle-activation summary of an ADL task
#'
#' Summarises how an agonist-antagonist muscle pair drives one ADL task from an
#' EMG envelope: the co-contraction index (via [emgCoContraction()]), each
#' muscle's peak and mean activation (as percent of a maximal reference when
#' `mvc_data` is given, else raw envelope units), and the fraction of the task
#' each muscle is active. Higher antagonist co-contraction means less efficient,
#' more guarded movement.
#'
#' @param x A `PhysioExperiment` of EMG envelope data (rectified / smoothed).
#' @param agonist,antagonist Channel index or `colData$label` of the pair.
#' @param mvc_data Optional `PhysioExperiment` of maximal-effort (MVC) data with
#'   the same channels; enables percent-MVC amplitudes.
#' @param task ADL task realised: `"reaching"` (d445), `"drinking"` (d560),
#'   `"feeding"` (d550), `"dressing"` (d540) or `"grooming"` (d520).
#' @param cc_method Co-contraction index method (see [emgCoContraction()]).
#' @param active_frac Activation threshold as a fraction of a muscle's own peak
#'   (default 0.2) for the active-time summary.
#' @param assay_name EMG assay to use (default: the first assay).
#' @return an `emg_adl_activation` list: `task`, `icf_code`,
#'   `cocontraction_index`, `unit` (`"%MVC"` or `"envelope"`), per-muscle
#'   `agonist_peak`/`agonist_mean`/`antagonist_peak`/`antagonist_mean` and
#'   `agonist_active_frac`/`antagonist_active_frac`.
#' @seealso [emgCoContraction()], [emgAmplitudeNormalize()]
#' @export
#' @examples
#' set.seed(1)
#' n <- 1000; t <- seq_len(n)
#' ag <- exp(-((t - 500) / 150)^2) + rnorm(n, 0, 0.02)     # agonist burst
#' an <- 0.2 * exp(-((t - 500) / 150)^2) + rnorm(n, 0, 0.02) # low antagonist
#' pe <- PhysioCore::PhysioExperiment(
#'   assays = S4Vectors::SimpleList(envelope = cbind(BIC = ag, TRI = an)),
#'   colData = S4Vectors::DataFrame(label = c("BIC", "TRI"), type = "EMG"),
#'   samplingRate = 1000)
#' emgADLActivation(pe, "BIC", "TRI", task = "drinking")$cocontraction_index
emgADLActivation <- function(x, agonist, antagonist, mvc_data = NULL,
                             task = c("reaching", "drinking", "feeding",
                                      "dressing", "grooming"),
                             cc_method = c("falconer_winter", "rudolph", "frost"),
                             active_frac = 0.2, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  task <- match.arg(task); cc_method <- match.arg(cc_method)
  if (is.null(assay_name)) assay_name <- SummarizedExperiment::assayNames(x)[1]

  cc <- emgCoContraction(x, agonist = agonist, antagonist = antagonist,
                         method = cc_method, assay_name = assay_name)

  data <- as.matrix(SummarizedExperiment::assay(x, assay_name))
  ai <- .emg_resolve_ch(x, agonist); bi <- .emg_resolve_ch(x, antagonist)
  ag <- data[, ai]; an <- data[, bi]; unit <- "envelope"
  if (!is.null(mvc_data)) {
    mvc <- as.matrix(SummarizedExperiment::assay(mvc_data, assay_name))
    ag <- ag / max(mvc[, ai], na.rm = TRUE) * 100
    an <- an / max(mvc[, bi], na.rm = TRUE) * 100
    unit <- "%MVC"
  }
  afrac <- function(v) { m <- max(v, na.rm = TRUE)
    if (!is.finite(m) || m <= 0) NA_real_ else mean(v > active_frac * m, na.rm = TRUE) }

  out <- list(
    task = task, icf_code = .emg_adl_icf(task),
    cocontraction_index = as.numeric(cc$summary), unit = unit,
    agonist_peak = max(ag, na.rm = TRUE), agonist_mean = mean(ag, na.rm = TRUE),
    antagonist_peak = max(an, na.rm = TRUE),
    antagonist_mean = mean(an, na.rm = TRUE),
    agonist_active_frac = afrac(ag), antagonist_active_frac = afrac(an))
  class(out) <- "emg_adl_activation"
  out
}

#' @export
print.emg_adl_activation <- function(x, ...) {
  cat(sprintf("<emg_adl_activation> %s (ICF %s)\n", x$task, x$icf_code))
  cat(sprintf("  co-contraction index %.3f\n", x$cocontraction_index))
  cat(sprintf("  agonist  peak %.1f / mean %.1f %s (active %.0f%%)\n",
              x$agonist_peak, x$agonist_mean, x$unit, 100 * x$agonist_active_frac))
  cat(sprintf("  antagon. peak %.1f / mean %.1f %s (active %.0f%%)\n",
              x$antagonist_peak, x$antagonist_mean, x$unit,
              100 * x$antagonist_active_frac))
  invisible(x)
}
