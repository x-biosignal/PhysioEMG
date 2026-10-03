library(testthat)
library(PhysioEMG)

test_that("coordinationNetwork dispatches to each estimator and matches it", {
  set.seed(1)
  pe <- make_emg(n_time = 2048, n_channels = 5, sr = 1000)

  cn <- coordinationNetwork(pe, method = "coherence", freq_band = c(20, 150), nperseg = 256L)
  direct <- emgCoherenceNetwork(pe, freq_band = c(20, 150), nperseg = 256L)
  expect_equal(dim(cn$network), c(5, 5))
  expect_equal(cn$network, direct$network)             # same result as the direct call

  expect_equal(dim(coordinationNetwork(pe, "partial_coherence", freq_band = c(20, 150))$network), c(5, 5))
  expect_equal(dim(coordinationNetwork(pe, "wpli", freq_band = c(20, 150), debiased = TRUE)$network), c(5, 5))
  expect_equal(dim(coordinationNetwork(pe, "directed_gc", max_lag = 5L)$network), c(5, 5))
})

test_that("coordinationNetwork rejects an unknown method", {
  set.seed(2)
  pe <- make_emg(n_time = 1024, n_channels = 3, sr = 1000)
  expect_error(coordinationNetwork(pe, method = "not_a_method"))
})

test_that("splitHalfStability returns per-method edge + node reproducibility", {
  set.seed(3)
  pe <- make_emg(n_time = 4096, n_channels = 5, sr = 1000)

  sh <- splitHalfStability(pe, methods = c("coherence", "wpli"),
                           freq_band = c(20, 150), nperseg = 256L)
  expect_s3_class(sh, "split_half_stability")
  expect_equal(nrow(sh), 2L)
  expect_true(all(c("method", "split_half_edge_correlation",
                    "split_half_node_strength_correlation",
                    "first_half_mean_edge", "second_half_mean_edge") %in% names(sh)))
  expect_setequal(sh$method, c("coherence", "wpli"))
  expect_true(all(sh$split_half_edge_correlation >= -1 & sh$split_half_edge_correlation <= 1))
  expect_true(all(sh$split_half_node_strength_correlation >= -1 &
                  sh$split_half_node_strength_correlation <= 1))
})

test_that("splitHalfStability errors on a too-short recording", {
  set.seed(4)
  pe <- make_emg(n_time = 2, n_channels = 3, sr = 1000)
  expect_error(splitHalfStability(pe), "too short")
})
