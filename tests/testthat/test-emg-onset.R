library(testthat)
library(PhysioEMG)

test_that("emgOnsetDetect with hodges_bui finds onset in burst", {
  pe <- make_emg(n_time = 2000, n_channels = 1, sr = 1000)

  result <- emgOnsetDetect(pe, method = "hodges_bui", threshold_sd = 3)

  expect_type(result, "list")
  expect_true("onsets" %in% names(result))
  expect_true("offsets" %in% names(result))
  expect_s3_class(result$onsets, "data.frame")
  expect_true(all(c("channel", "sample", "time_sec") %in% names(result$onsets)))
  expect_true(nrow(result$onsets) > 0)
  expect_true(all(result$onsets$sample >= 1 & result$onsets$sample <= 2000))
})

test_that("emgOnsetDetect with teager_kaiser works", {
  pe <- make_emg(n_time = 2000, n_channels = 1, sr = 1000)

  result <- emgOnsetDetect(pe, method = "teager_kaiser")

  expect_type(result, "list")
  expect_true(nrow(result$onsets) > 0)
})

test_that("emgOnsetDetect handles multi-channel data", {
  pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)

  result <- emgOnsetDetect(pe, method = "hodges_bui")

  expect_true(length(unique(result$onsets$channel)) > 0)
})

test_that("emgOnsetDetect detects offset after burst", {
  pe <- make_emg(n_time = 2000, n_channels = 1, sr = 1000)

  result <- emgOnsetDetect(pe, method = "hodges_bui")

  expect_true(nrow(result$offsets) > 0)
  expect_true(all(result$offsets$sample > result$onsets$sample))
})
