library(testthat)
library(PhysioEMG)

# Build a PhysioExperiment from a raw matrix (helper-local to this file).
pe_from_matrix <- function(m, sr = 1000) {
  if (is.null(dim(m))) m <- matrix(m, ncol = 1)
  PhysioExperiment(
    assays = list(raw = m),
    colData = S4Vectors::DataFrame(
      label = paste0("EMG", seq_len(ncol(m))),
      type = rep("EMG", ncol(m))),
    samplingRate = sr)
}

test_that("emgAmplitudeFeatures returns per-window per-channel rows + columns", {
  pe <- make_emg(n_time = 2000, n_channels = 3, sr = 1000)
  feats <- emgAmplitudeFeatures(pe, window_sec = 0.5, overlap = 0.5)

  expect_s3_class(feats, "data.frame")
  expect_true(all(c("channel", "window", "time_sec",
                    "arv", "mav", "iemg", "rms") %in% names(feats)))
  expect_setequal(unique(feats$channel), 1:3)
  # window start times increase by step within a channel
  ch1 <- feats[feats$channel == 1, ]
  expect_true(all(diff(ch1$time_sec) > 0))
})

test_that("iEMG of a constant |x|=1 signal over T seconds equals T", {
  sr <- 1000
  T_sec <- 2
  n <- T_sec * sr                       # 2000 samples -> exactly 2 s
  pe <- pe_from_matrix(rep(1, n), sr = sr)
  # single window spanning the whole signal (win_samples == n_time)
  feats <- emgAmplitudeFeatures(pe, features = "iemg",
                                window_sec = T_sec, overlap = 0)
  expect_equal(nrow(feats), 1L)
  expect_equal(feats$iemg[1], T_sec, tolerance = 1e-6)

  # also holds for a negative constant (rectified): |x| = 1
  pe2 <- pe_from_matrix(rep(-1, n), sr = sr)
  f2 <- emgAmplitudeFeatures(pe2, features = "iemg",
                             window_sec = T_sec, overlap = 0)
  expect_equal(f2$iemg[1], T_sec, tolerance = 1e-6)
})

test_that("ARV == MAV, and RMS >= ARV for a non-constant signal", {
  pe <- make_emg(n_time = 3000, n_channels = 4, sr = 1000)
  feats <- emgAmplitudeFeatures(pe, window_sec = 0.25, overlap = 0.5)

  # ARV and MAV are the same unweighted mean(|x|)
  expect_equal(feats$arv, feats$mav)

  # Cauchy-Schwarz: sqrt(mean(x^2)) >= mean(|x|), strict for non-constant x
  expect_true(all(feats$rms >= feats$arv - 1e-12))
  expect_true(any(feats$rms > feats$arv))  # burst windows are non-constant
})

test_that("windowing matches emgFatigue window/step count on identical args", {
  pe <- make_emg_fatigue(n_time = 8000, sr = 1000)
  args <- list(window_sec = 1.0, overlap = 0.5)

  feats <- emgAmplitudeFeatures(pe, window_sec = args$window_sec,
                                overlap = args$overlap)
  fat <- emgFatigue(pe, window_sec = args$window_sec, overlap = args$overlap)

  # same number of windows overall and per channel
  expect_equal(nrow(feats), nrow(fat))
  expect_equal(table(feats$channel), table(fat$channel))
  # identical window indices and start times
  expect_equal(feats$window, fat$window)
  expect_equal(feats$time_sec, fat$time_sec)

  # a second, different set of args stays consistent
  feats2 <- emgAmplitudeFeatures(pe, window_sec = 0.5, overlap = 0.25)
  fat2 <- emgFatigue(pe, window_sec = 0.5, overlap = 0.25)
  expect_equal(nrow(feats2), nrow(fat2))
  expect_equal(feats2$time_sec, fat2$time_sec)
})

test_that("feature selection returns only requested columns in canonical order", {
  pe <- make_emg(n_time = 1500, n_channels = 2, sr = 1000)

  only_rms <- emgAmplitudeFeatures(pe, features = "rms", window_sec = 0.5)
  expect_equal(names(only_rms), c("channel", "window", "time_sec", "rms"))

  # order of `features` argument does not change column order
  two <- emgAmplitudeFeatures(pe, features = c("rms", "arv"), window_sec = 0.5)
  expect_equal(names(two), c("channel", "window", "time_sec", "arv", "rms"))
})

test_that("input validation and empty-window handling", {
  pe <- make_emg(n_time = 500, n_channels = 2, sr = 1000)

  expect_error(emgAmplitudeFeatures(pe, features = "nope"))
  expect_error(emgAmplitudeFeatures(pe, window_sec = 0))
  expect_error(emgAmplitudeFeatures(pe, overlap = 1))
  expect_error(emgAmplitudeFeatures(pe, overlap = -0.1))

  # window longer than the signal -> 0-row frame with the right columns
  empty <- emgAmplitudeFeatures(pe, window_sec = 10, overlap = 0.5)
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("channel", "window", "time_sec",
                    "arv", "mav", "iemg", "rms") %in% names(empty)))
})
