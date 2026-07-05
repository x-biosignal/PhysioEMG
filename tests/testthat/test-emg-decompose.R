library(testthat)
library(PhysioEMG)

test_that("muscleSynergy with NMF decomposes correctly", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)

  result <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

  expect_type(result, "list")
  expect_true(all(c("W", "H", "vaf", "method") %in% names(result)))
  expect_equal(nrow(result$W), 3)
  expect_equal(ncol(result$W), 8)
  expect_equal(nrow(result$H), 1000)
  expect_equal(ncol(result$H), 3)
  expect_gt(result$vaf, 0.5)
})

test_that("muscleSynergy with PCA works", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)

  result <- muscleSynergy(pe, n_synergies = 3, method = "pca")

  expect_type(result, "list")
  expect_equal(nrow(result$W), 3)
  expect_equal(ncol(result$W), 8)
})

test_that("synergyReconstruct reconstructs data", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)
  result <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

  recon <- synergyReconstruct(result, n_synergies = 2)

  expect_type(recon, "list")
  expect_true("reconstructed" %in% names(recon))
  expect_equal(dim(recon$reconstructed), c(1000, 8))
  expect_lt(recon$vaf, result$vaf)
})

test_that("synergyCompare computes similarity", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)
  r1 <- muscleSynergy(pe, n_synergies = 3, method = "nmf")
  set.seed(123)
  r2 <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

  sim <- synergyCompare(r1, r2)

  expect_s3_class(sim, "data.frame")
  expect_true("correlation" %in% names(sim))
  expect_true(all(abs(sim$correlation) > 0.3))
})
