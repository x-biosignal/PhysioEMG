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

# Build a PhysioExperiment from a raw amplitude matrix (local helper).
amp_pe <- function(m, sr = 1000) {
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  PhysioExperiment(
    assays = list(raw = m),
    colData = S4Vectors::DataFrame(
      label = paste0("EMG", seq_len(ncol(m))),
      type = rep("EMG", ncol(m))),
    samplingRate = sr)
}

test_that("RVC of the reference task by itself yields channel means ~1.0", {
  set.seed(11)
  amp <- matrix(abs(rnorm(3000 * 3, sd = 0.4)) + 0.05, nrow = 3000, ncol = 3)
  pe <- amp_pe(amp)

  res <- emgAmplitudeNormalize(pe, method = "rvc", rvc_data = pe)
  nd <- SummarizedExperiment::assay(res, "normalized")
  expect_equal(unname(colMeans(nd)), rep(1, 3), tolerance = 1e-6)
})

test_that("RVC errors when rvc_data channel count mismatches (mirrors MVC)", {
  pe <- amp_pe(matrix(abs(rnorm(2000 * 2)), 2000, 2))
  ref3 <- amp_pe(matrix(abs(rnorm(2000 * 3)), 2000, 3))

  expect_error(
    emgAmplitudeNormalize(pe, method = "rvc", rvc_data = ref3),
    "same number of channels")
  expect_error(
    emgAmplitudeNormalize(pe, method = "rvc"),
    "rvc_data is required")
})

test_that("%RVC is scale-invariant to a global gain on signal and reference", {
  set.seed(12)
  amp <- matrix(abs(rnorm(2500 * 2, sd = 0.3)) + 0.02, nrow = 2500, ncol = 2)
  ref <- matrix(abs(rnorm(1800 * 2, sd = 0.5)) + 0.02, nrow = 1800, ncol = 2)
  g <- 7.3

  n1 <- SummarizedExperiment::assay(
    emgAmplitudeNormalize(amp_pe(amp), method = "rvc",
                          rvc_data = amp_pe(ref)), "normalized")
  n2 <- SummarizedExperiment::assay(
    emgAmplitudeNormalize(amp_pe(g * amp), method = "rvc",
                          rvc_data = amp_pe(g * ref)), "normalized")
  expect_equal(n1, n2, tolerance = 1e-10)
})

test_that("RVC rvc_window restricts the reference to a sub-window", {
  # reference amplitude steps up in the second half; a window over the first
  # half yields a smaller reference mean (hence larger %RVC) than the whole trial
  ref <- matrix(c(rep(1, 1000), rep(3, 1000)), ncol = 1)
  pe <- amp_pe(matrix(rep(2, 500), ncol = 1))

  full <- SummarizedExperiment::assay(
    emgAmplitudeNormalize(pe, method = "rvc", rvc_data = amp_pe(ref)),
    "normalized")[1, 1]                      # 2 / mean(c(1,3)) = 1
  firsthalf <- SummarizedExperiment::assay(
    emgAmplitudeNormalize(pe, method = "rvc", rvc_data = amp_pe(ref),
                          rvc_window = c(0, 0.999)), "normalized")[1, 1]  # 2/1 = 2
  expect_equal(full, 1, tolerance = 1e-8)
  expect_equal(firsthalf, 2, tolerance = 1e-3)
})

test_that("dynamic_peak normalizes by a centered moving maximum", {
  set.seed(13)
  amp <- matrix(abs(rnorm(2000 * 2, sd = 0.3)) + 0.01, nrow = 2000, ncol = 2)
  pe <- amp_pe(amp)

  res <- emgAmplitudeNormalize(pe, method = "dynamic_peak", rvc_window = 0.1)
  nd <- SummarizedExperiment::assay(res, "normalized")
  expect_equal(dim(nd), dim(amp))
  # x >= 0 divided by the moving max of |x| lies in [0, 1]
  expect_true(all(nd >= -1e-9 & nd <= 1 + 1e-9))
  # the global-max sample normalizes to 1
  expect_equal(max(nd[, 1]), 1, tolerance = 1e-9)
})
