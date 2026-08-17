library(testthat)
library(PhysioEMG)

emg_pe <- function(nch = 2, n = 200, sr = 1000) {
  PhysioExperiment(
    assays = list(raw = matrix(rnorm(n * nch), ncol = nch)),
    colData = S4Vectors::DataFrame(label = paste0("ch", seq_len(nch)),
                                   type = rep("EMG", nch)),
    samplingRate = sr)
}

test_that("seniamMuscles returns the bundled catalog", {
  cat <- seniamMuscles()
  expect_s3_class(cat, "data.frame")
  expect_true(all(c("muscle", "location", "orientation", "ied_mm") %in%
                    names(cat)))
  expect_gt(nrow(cat), 20)
  expect_true(all(cat$ied_mm > 0))
  expect_true("Biceps brachii" %in% cat$muscle)
})

test_that("setEMGElectrode round-trips through getEMGElectrode and keeps validity", {
  pe <- emg_pe(2)
  pe <- setEMGElectrode(pe, muscle = c("Biceps brachii", "biceps_femoris"),
                        side = "right")

  # PhysioExperiment stays valid after writing colData
  expect_true(methods::validObject(pe))
  expect_s4_class(pe, "PhysioExperiment")

  ge <- getEMGElectrode(pe)
  expect_s3_class(ge, "data.frame")
  expect_true(all(c("channel", "muscle", "placement", "orientation",
                    "ied_mm", "reference", "side") %in% names(ge)))
  # muscle names canonicalized to the catalog spelling
  expect_equal(ge$muscle, c("Biceps brachii", "Biceps femoris"))
  # recommended IED, orientation and placement back-filled from the catalog
  expect_equal(ge$ied_mm, c(20, 20))
  expect_true(all(!is.na(ge$orientation)))
  expect_true(all(!is.na(ge$placement)))
  expect_equal(ge$side, c("right", "right"))
})

test_that("explicit ied_mm and placement override the catalog defaults", {
  pe <- setEMGElectrode(emg_pe(1), muscle = "Tibialis anterior",
                        ied_mm = 10, placement = "custom site")
  ge <- getEMGElectrode(pe)
  expect_equal(ge$ied_mm[1], 10)
  expect_equal(ge$placement[1], "custom site")
})

test_that("an unknown muscle warns with near matches", {
  pe <- emg_pe(1)
  expect_warning(setEMGElectrode(pe, muscle = "bicep"), "SENIAM catalog")
  expect_warning(setEMGElectrode(pe, muscle = "bicep"), "did you mean")

  # the provided name is still stored despite the warning
  ge <- suppressWarnings(getEMGElectrode(setEMGElectrode(pe, muscle = "bicep")))
  expect_equal(ge$muscle[1], "bicep")
})

test_that("channel targeting annotates a subset only", {
  pe <- emg_pe(3)
  pe <- setEMGElectrode(pe, muscle = "Vastus lateralis", channel = 2)
  ge <- getEMGElectrode(pe)
  expect_equal(ge$muscle[2], "Vastus lateralis")
  expect_true(is.na(ge$muscle[1]) && is.na(ge$muscle[3]))

  # labels also work as channel selectors
  pe <- setEMGElectrode(pe, muscle = "Rectus femoris", channel = "ch1")
  expect_equal(getEMGElectrode(pe)$muscle[1], "Rectus femoris")
})

test_that("the SENIAM IED is consumable by emgMFCV()", {
  sr <- 2000
  set.seed(1)
  prox <- as.numeric(stats::filter(rnorm(4000), rep(1 / 5, 5), sides = 2))
  prox[is.na(prox)] <- 0
  d <- 5
  dist <- c(rep(0, d), prox[seq_len(4000 - d)])
  pe <- PhysioExperiment(
    assays = list(raw = cbind(prox, dist)),
    colData = S4Vectors::DataFrame(label = c("prox", "dist")),
    samplingRate = sr)
  pe <- setEMGElectrode(pe, muscle = "Biceps brachii")   # 20 mm IED

  ied <- getEMGElectrode(pe)$ied_mm[1]
  expect_equal(ied, 20)
  res <- emgMFCV(pe, electrode_pairs = list(c(1, 2)), ied_mm = ied)
  expect_equal(res$velocity_m_s[1], (ied / 1000) / (d / sr), tolerance = 0.05)
})
