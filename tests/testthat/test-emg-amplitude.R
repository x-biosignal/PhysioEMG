library(testthat)
library(PhysioEMG)

test_that("emgEnvelope with RMS method returns correct dimensions", {
  pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)

  result <- emgEnvelope(pe, method = "rms", window_ms = 50)

  expect_s4_class(result, "PhysioExperiment")
  expect_true("envelope" %in% SummarizedExperiment::assayNames(result))
  env <- SummarizedExperiment::assay(result, "envelope")
  expect_equal(nrow(env), 2000)
  expect_equal(ncol(env), 4)
  # Envelope should be non-negative
  expect_true(all(env >= 0, na.rm = TRUE))
})

test_that("emgEnvelope with Hilbert method works", {
  pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)

  result <- emgEnvelope(pe, method = "hilbert")

  expect_s4_class(result, "PhysioExperiment")
  env <- SummarizedExperiment::assay(result, "envelope")
  expect_true(all(env >= 0, na.rm = TRUE))
})

test_that("emgEnvelope with lowpass method works", {
  pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)

  result <- emgEnvelope(pe, method = "lowpass", cutoff = 6)

  expect_s4_class(result, "PhysioExperiment")
  expect_true("envelope" %in% SummarizedExperiment::assayNames(result))
})

test_that("emgEnvelope detects higher amplitude in burst region", {
  pe <- make_emg(n_time = 2000, n_channels = 1, sr = 1000)

  result <- emgEnvelope(pe, method = "rms", window_ms = 100)
  env <- SummarizedExperiment::assay(result, "envelope")

  # Burst is at 30-70%, rest is 0-30% and 70-100%
  burst_mean <- mean(env[600:1400, 1], na.rm = TRUE)
  rest_mean <- mean(env[c(1:500, 1500:2000), 1], na.rm = TRUE)
  expect_gt(burst_mean, rest_mean * 2)
})

test_that("emgAmplitudeNormalize with MVC works", {
  pe <- make_emg(n_time = 2000, n_channels = 2, sr = 1000)
  mvc_pe <- make_emg(n_time = 1000, n_channels = 2, sr = 1000)

  # First get envelopes
  pe <- emgEnvelope(pe, method = "rms", window_ms = 50)
  mvc_pe <- emgEnvelope(mvc_pe, method = "rms", window_ms = 50)

  result <- emgAmplitudeNormalize(pe, method = "mvc", mvc_data = mvc_pe,
                                   assay_name = "envelope")

  expect_s4_class(result, "PhysioExperiment")
  expect_true("normalized" %in% SummarizedExperiment::assayNames(result))
})

test_that("emgAmplitudeNormalize with peak works", {
  pe <- make_emg(n_time = 2000, n_channels = 2, sr = 1000)
  pe <- emgEnvelope(pe, method = "rms", window_ms = 50)

  result <- emgAmplitudeNormalize(pe, method = "peak", assay_name = "envelope")

  norm_data <- SummarizedExperiment::assay(result, "normalized")
  # Peak-normalized: max should be 1.0 per channel
  for (ch in seq_len(ncol(norm_data))) {
    expect_equal(max(norm_data[, ch], na.rm = TRUE), 1.0, tolerance = 1e-10)
  }
})
