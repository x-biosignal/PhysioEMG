library(testthat)
library(PhysioEMG)

# Build a two-channel PhysioExperiment where the distal channel is the proximal
# signal delayed by `d` samples (a synthetic propagating action potential).
delayed_pair_pe <- function(d, n = 4000, sr = 2000, seed = 1) {
  set.seed(seed)
  prox <- as.numeric(stats::filter(rnorm(n), rep(1 / 5, 5), sides = 2))
  prox[is.na(prox)] <- 0
  dist <- c(rep(0, d), prox[seq_len(n - d)])
  PhysioExperiment(
    assays = list(raw = cbind(prox, dist)),
    colData = S4Vectors::DataFrame(label = c("prox", "dist"),
                                   type = c("EMG", "EMG")),
    samplingRate = sr)
}

test_that("MFCV recovers v = IED / delay for a synthetic delayed pair (xcorr)", {
  sr <- 2000; d <- 5; ied <- 10           # 10 mm, 5-sample delay
  pe <- delayed_pair_pe(d = d, sr = sr)

  res <- emgMFCV(pe, electrode_pairs = list(c(1, 2)), ied_mm = ied)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("pair", "ch1", "ch2", "delay_ms", "velocity_m_s",
                    "quality") %in% names(res)))

  delay_s <- d / sr
  v_true <- (ied / 1000) / delay_s        # 0.01 / 0.0025 = 4 m/s
  expect_equal(res$velocity_m_s[1], v_true, tolerance = 0.05)  # within 5%
  expect_equal(abs(res$delay_ms[1]), delay_s * 1000, tolerance = 0.05)
})

test_that("MFCV xcorr recovers a non-integer delay via sub-sample interpolation", {
  # larger delay, different IED; xcorr peak + parabolic interpolation
  sr <- 2000; d <- 8; ied <- 20
  pe <- delayed_pair_pe(d = d, n = 6000, sr = sr, seed = 3)
  res <- emgMFCV(pe, electrode_pairs = list(c(1, 2)), ied_mm = ied)
  v_true <- (ied / 1000) / (d / sr)       # 0.02 / 0.004 = 5 m/s
  expect_equal(res$velocity_m_s[1], v_true, tolerance = 0.05)
})

test_that("MFCV phase method also recovers the velocity within tolerance", {
  sr <- 2000; d <- 5; ied <- 10
  pe <- delayed_pair_pe(d = d, sr = sr, seed = 5)
  res <- emgMFCV(pe, electrode_pairs = list(c(1, 2)), ied_mm = ied,
                 method = "phase", band = c(20, 400))
  v_true <- (ied / 1000) / (d / sr)
  expect_equal(res$velocity_m_s[1], v_true, tolerance = 0.1)
})

test_that("MFCV accepts a matrix of pairs and per-pair IED, and validates input", {
  sr <- 2000
  set.seed(2)
  prox <- as.numeric(stats::filter(rnorm(4000), rep(1 / 5, 5), sides = 2))
  prox[is.na(prox)] <- 0
  d1 <- 4; d2 <- 6
  m <- cbind(prox,
             c(rep(0, d1), prox[seq_len(4000 - d1)]),
             c(rep(0, d2), prox[seq_len(4000 - d2)]))
  pe <- PhysioExperiment(assays = list(raw = m), samplingRate = sr)

  pairs <- rbind(c(1, 2), c(1, 3))
  res <- emgMFCV(pe, electrode_pairs = pairs, ied_mm = c(10, 15))
  expect_equal(nrow(res), 2L)
  expect_equal((10 / 1000) / (d1 / sr), res$velocity_m_s[1], tolerance = 0.05)
  expect_equal((15 / 1000) / (d2 / sr), res$velocity_m_s[2], tolerance = 0.05)

  expect_error(emgMFCV(pe, electrode_pairs = list(c(1, 5)), ied_mm = 10),
               "outside")
  expect_error(emgMFCV(pe, electrode_pairs = matrix(1:3, ncol = 3),
                       ied_mm = 10), "2 columns")
})
