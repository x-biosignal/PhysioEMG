# Bridge to the PhysioHDEMG package: surface the HD-sEMG motor-unit
# decomposition entry point in the PhysioEMG namespace without hard-depending
# on PhysioHDEMG (it is only needed when the function is actually called).

#' Decompose HD-sEMG into motor units (PhysioHDEMG bridge)
#'
#' Convenience bridge that forwards to
#' \code{PhysioHDEMG::hdEMGDecompose()} for high-density surface EMG
#' motor-unit decomposition by convolutive blind source separation. The
#' PhysioHDEMG package must be installed. See
#' \code{\link[PhysioHDEMG:hdEMGDecompose]{PhysioHDEMG::hdEMGDecompose()}} for
#' the full argument list and return value.
#'
#' @param x HD-sEMG data: an `n_time x n_channels` matrix, a `PhysioExperiment`,
#'   or an `hdemg_sim` object.
#' @param ... Further arguments passed to `PhysioHDEMG::hdEMGDecompose()`.
#' @return An `hdemg_decomposition` object (see PhysioHDEMG).
#' @seealso [muscleSynergy()] for whole-muscle synergy decomposition.
#' @export
#' @examples
#' \dontrun{
#' sim <- PhysioHDEMG::make_hdemg_sim(n_units = 2, duration_sec = 3)
#' dec <- hdEMGDecompose(sim, n_units = 4)
#' }
hdEMGDecompose <- function(x, ...) {
  if (!requireNamespace("PhysioHDEMG", quietly = TRUE)) {
    stop("hdEMGDecompose() requires the PhysioHDEMG package. Install it with ",
         "install.packages('PhysioHDEMG', repos = 'https://x-biosignal.r-universe.dev').",
         call. = FALSE)
  }
  PhysioHDEMG::hdEMGDecompose(x, ...)
}
