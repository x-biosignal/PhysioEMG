library(testthat)
library(PhysioEMG)

# Force a headless render so device/build errors surface without a graphics
# device (as in CI).
built <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_silent(invisible(ggplot2::ggplot_build(p)))
  invisible(p)
}

test_that("plotEmgEnvelope returns a renderable ggplot", {
  pe <- make_emg(n_time = 2000, n_channels = 1)
  built(plotEmgEnvelope(pe))
  built(plotEmgEnvelope(pe, qc = FALSE))
})

test_that("plotEmgOnset returns a renderable ggplot", {
  pe <- make_emg_contraction()
  built(plotEmgOnset(pe))
  built(plotEmgOnset(pe, method = "bonato"))
})

test_that("plotEmgFatigue returns a renderable ggplot for mdf and mnf", {
  pe <- make_emg_fatigue()
  built(plotEmgFatigue(pe, feature = "mdf"))
  built(plotEmgFatigue(pe, feature = "mnf"))
})

test_that("plotSynergy renders for muscleSynergy and muscleSynergyOrder results", {
  pe <- make_emg_synergy(n_time = 400, n_channels = 6, n_synergies = 3)
  built(plotSynergy(muscleSynergy(pe, n_synergies = 3, method = "nmf", seed = 1)))
  ord <- muscleSynergyOrder(pe, max_synergies = 5, seed = 1)
  built(plotSynergy(ord))
})

test_that("plotSynergy VAF panel matches the muscleSynergyOrder VAF values", {
  pe <- make_emg_synergy(n_time = 400, n_channels = 6, n_synergies = 3)
  ord <- muscleSynergyOrder(pe, max_synergies = 5, seed = 1)
  p <- plotSynergy(ord)
  b <- ggplot2::ggplot_build(p)

  # locate the VAF panel: the layer whose y-values equal the VAF curve
  plotted <- NULL
  for (layer in b$data) {
    if (all(c("x", "y") %in% names(layer)) &&
        nrow(layer) == nrow(ord$vaf_curve) &&
        isTRUE(all.equal(sort(layer$y), sort(ord$vaf_curve$vaf),
                         tolerance = 1e-8))) {
      plotted <- layer
      break
    }
  }
  expect_false(is.null(plotted))
  expect_equal(sort(plotted$y), sort(ord$vaf_curve$vaf), tolerance = 1e-8)
})

test_that("plotCoactivation returns a renderable ggplot", {
  t <- seq_len(1000)
  a <- exp(-((t - 400) / 120)^2)
  b <- exp(-((t - 600) / 120)^2)
  pe <- PhysioExperiment(
    assays = list(env = cbind(a, b)),
    colData = S4Vectors::DataFrame(label = c("TA", "GAS")),
    samplingRate = 1000)
  built(plotCoactivation(pe, "TA", "GAS"))
  built(plotCoactivation(pe, 1, 2, method = "rudolph"))
})
