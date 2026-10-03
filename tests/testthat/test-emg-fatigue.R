library(testthat)
library(PhysioEMG)

test_that("emgFatigue returns frequency metrics over time", {
  pe <- make_emg_fatigue(n_time = 10000, sr = 1000)

  result <- emgFatigue(pe, window_sec = 1.0)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("channel", "window", "time_sec", "median_freq", "mean_freq") %in% names(result)))
  expect_true(nrow(result) > 0)
  # Median frequency should decrease over time (fatigue)
  ch1 <- result[result$channel == 1, ]
  first_half <- mean(ch1$median_freq[1:(nrow(ch1) %/% 2)], na.rm = TRUE)
  second_half <- mean(ch1$median_freq[(nrow(ch1) %/% 2 + 1):nrow(ch1)], na.rm = TRUE)
  expect_gt(first_half, second_half)
})

test_that("emgFatigueIndex returns summary metric", {
  pe <- make_emg_fatigue(n_time = 10000, sr = 1000)

  idx <- emgFatigueIndex(pe)

  expect_s3_class(idx, "data.frame")
  expect_true(all(c("channel", "fatigue_index", "initial_mdf", "final_mdf") %in% names(idx)))
  # Fatigue index should be < 1 (frequency decreased)
  expect_lt(idx$fatigue_index[1], 1)
})

test_that("emgSpectralMoments returns moment values", {
  pe <- make_emg_fatigue(n_time = 10000, sr = 1000)

  moments <- emgSpectralMoments(pe, window_sec = 1.0)

  expect_s3_class(moments, "data.frame")
  expect_true(all(c("channel", "window", "m0", "m1", "m2") %in% names(moments)))
  expect_true(all(moments$m0 >= 0))
})

# ---- emgFatigueSlope --------------------------------------------------------

test_that("emgFatigueSlope MDF slope is negative and significant under fatigue", {
  pe <- make_emg_fatigue(n_time = 10000, sr = 1000)
  res <- emgFatigueSlope(pe, feature = "mdf")

  expect_s3_class(res, "data.frame")
  expect_true(all(c("channel", "slope_hz_per_min", "norm_slope_pct_per_min",
                    "intercept_hz", "r_squared", "p_value") %in% names(res)))
  expect_lt(res$slope_hz_per_min[1], 0)     # frequency falls -> fatigue
  expect_lt(res$p_value[1], 0.05)
  expect_lt(res$norm_slope_pct_per_min[1], 0)
})

test_that("emgFatigueSlope recovers an imposed linear MDF ramp within 1%", {
  sr <- 2000
  dur <- 16
  n <- sr * dur
  t <- seq_len(n) / sr
  rate <- -6                                 # Hz per second
  f0 <- 190
  f_inst <- f0 + rate * t                    # 190 -> 94 Hz
  set.seed(42)
  sig <- sin(2 * pi * cumsum(f_inst) / sr) + rnorm(n, sd = 0.02)  # chirp + noise
  pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                         samplingRate = sr)

  res <- emgFatigueSlope(pe, feature = "mdf", window_sec = 1.0, overlap = 0.5)
  imposed_hz_per_min <- rate * 60            # -360
  expect_equal(res$slope_hz_per_min[1], imposed_hz_per_min,
               tolerance = 0.01)             # within 1%
})

# ---- emgDimitrovIndex -------------------------------------------------------

test_that("emgDimitrovIndex returns per-window moments and a slope attribute", {
  pe <- make_emg_fatigue(n_time = 8000, sr = 1000)
  di <- emgDimitrovIndex(pe, window_sec = 1.0)

  expect_s3_class(di, "data.frame")
  expect_true(all(c("channel", "window", "time_sec", "m_minus1",
                    "m0", "m5", "finsm5") %in% names(di)))
  sl <- attr(di, "finsm5_slope")
  expect_s3_class(sl, "data.frame")
  expect_true(all(c("channel", "slope_per_min") %in% names(sl)))
  # fatiguing signal -> FInsm5 rises over time
  expect_gt(sl$slope_per_min[1], 0)
})

test_that("FInsm5 increases monotonically as the spectrum compresses to low freq", {
  sr <- 1000
  n <- 2000                                  # one 2 s window
  finsm5 <- vapply(c(120, 100, 80, 60, 40), function(f) {
    sig <- sin(2 * pi * f * seq_len(n) / sr)
    pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                           samplingRate = sr)
    emgDimitrovIndex(pe, window_sec = 2.0)$finsm5[1]
  }, numeric(1))
  # decreasing dominant frequency -> strictly increasing FInsm5
  expect_true(all(diff(finsm5) > 0))
})
