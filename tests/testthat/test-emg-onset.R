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

# ---- Bonato double-threshold ------------------------------------------------

test_that("Bonato onset is within +/-30 ms of a known contraction onset", {
  set.seed(101)
  pe <- make_emg_contraction(n_time = 5000, sr = 1000,
                             contraction_start = 0.3)
  gt_sample <- as.integer(5000 * 0.3)          # 1500

  res <- emgOnsetDetect(pe, method = "bonato")
  expect_true(nrow(res$onsets) >= 1)
  expect_true(all(c("channel", "sample", "time_sec", "stat") %in%
                    names(res$onsets)))
  nearest <- min(abs(res$onsets$sample - gt_sample))
  expect_lt(nearest, 30)                        # 30 samples = 30 ms at 1000 Hz
  # offsets follow onsets and carry the detector statistic
  expect_true(all(res$offsets$sample > res$onsets$sample))
  expect_true(all(res$onsets$stat > 0))
})

# ---- AGLR ------------------------------------------------------------------

test_that("AGLR and Bonato agree within 50 ms on the same clean contraction", {
  set.seed(102)
  pe <- make_emg_contraction(n_time = 5000, sr = 1000)
  gt <- as.integer(5000 * 0.3)

  b <- emgOnsetDetect(pe, method = "bonato")$onsets$sample
  a <- emgOnsetDetect(pe, method = "aglr")$onsets$sample
  expect_true(length(a) >= 1 && length(b) >= 1)
  bs <- b[which.min(abs(b - gt))]
  as_ <- a[which.min(abs(a - gt))]
  expect_lt(abs(bs - as_), 50)                  # within 50 ms
  expect_lt(abs(as_ - gt), 30)                  # AGLR also close to truth
})

# ---- rectify-and-smooth conditioning (v0.2.1) -------------------------------

# Band-limited multi-sine "EMG": rich in zero-crossings, so the rectified
# trace dips below any suprabaseline threshold every half-cycle. Without
# smoothing, no contiguous run can reach min_duration_ms and the threshold
# methods detect nothing even at high SNR.
.bandlimited_burst_pe <- function(sr = 2000, dur_sec = 4, onset_sec = 1.5,
                                  burst_sec = 1, snr_lin = 10) {
  n <- as.integer(dur_sec * sr)
  t <- seq_len(n) / sr
  multisine <- function() {
    freqs <- runif(60, 20, 450)
    ph <- runif(60, 0, 2 * pi)
    rowSums(sapply(seq_along(freqs), function(k) sin(2 * pi * freqs[k] * t + ph[k])))
  }
  base <- multisine(); burst <- multisine()
  env <- numeric(n)
  on_i <- as.integer(onset_sec * sr)
  env[on_i:min(n, on_i + as.integer(burst_sec * sr))] <- 1
  x <- base / sd(base) + (burst / sd(burst)) * env * snr_lin
  PhysioExperiment(assays = list(raw = matrix(x, ncol = 1)),
                   samplingRate = sr)
}

test_that("threshold methods detect a band-limited burst at SNR 20 dB", {
  set.seed(104)
  pe <- .bandlimited_burst_pe(snr_lin = 10)     # 20 dB

  for (mth in c("hodges_bui", "teager_kaiser")) {
    res <- emgOnsetDetect(pe, method = mth)
    expect_gt(nrow(res$onsets), 0)
    nearest <- min(abs(res$onsets$time_sec - 1.5))
    expect_lt(nearest, 0.025)                   # within 25 ms of true onset
  }
})

test_that("smoothing tolerates NA samples in the signal", {
  set.seed(106)
  pe <- .bandlimited_burst_pe(snr_lin = 10)
  a <- SummarizedExperiment::assay(pe, "raw")
  a[c(600, 5000), 1] <- NA          # one NA in baseline, one inside the burst
  pe_na <- PhysioExperiment(assays = list(raw = a), samplingRate = 2000)

  for (mth in c("hodges_bui", "teager_kaiser")) {
    res <- emgOnsetDetect(pe_na, method = mth)
    expect_gt(nrow(res$onsets), 0)
    expect_lt(min(abs(res$onsets$time_sec - 1.5)), 0.025)
  }
})

test_that("smooth_ms narrower than one sample warns and falls back", {
  set.seed(107)
  x <- matrix(rnorm(400), ncol = 1)
  pe <- PhysioExperiment(assays = list(raw = x), samplingRate = 20)

  expect_warning(emgOnsetDetect(pe, method = "hodges_bui", smooth_ms = 25),
                 "smoothing is effectively disabled")
})

test_that("smooth_ms = 0 restores the legacy unsmoothed behavior", {
  set.seed(105)
  pe <- .bandlimited_burst_pe(snr_lin = 10)

  for (mth in c("hodges_bui", "teager_kaiser")) {
    res <- emgOnsetDetect(pe, method = mth, smooth_ms = 0)
    expect_equal(nrow(res$onsets), 0)           # legacy: burst goes undetected
  }
})

# ---- false-positive rate on baseline-only noise -----------------------------

test_that("false-positive rate is < 1 per 10 s on baseline-only noise", {
  set.seed(103)
  # 60 s of pure baseline noise; default settings
  noise <- matrix(rnorm(60000, sd = 0.02), ncol = 1)
  pe <- PhysioExperiment(assays = list(raw = noise), samplingRate = 1000)

  for (mth in c("bonato", "aglr")) {
    res <- emgOnsetDetect(pe, method = mth)
    # fewer than 6 detections in 60 s  ==  rate < 1 per 10 s
    expect_lt(nrow(res$onsets), 6)
  }
})
