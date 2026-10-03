#' @keywords internal
#' @importFrom stats cor embed fft prcomp rnorm sd
#' @importFrom PhysioExperiment PhysioExperiment
#' @importFrom PhysioExperiment defaultAssay samplingRate
#' @examples
#' # Simulate EMG, extract an RMS envelope, and decompose into synergies
#' pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
#' pe <- emgEnvelope(pe, method = "rms")
#' muscleSynergy(pe, n_synergies = 2, method = "nmf", seed = 1,
#'               assay_name = "envelope")$vaf
"_PACKAGE"
