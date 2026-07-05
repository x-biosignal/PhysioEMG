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
