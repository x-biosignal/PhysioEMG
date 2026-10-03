library(testthat)
library(PhysioEMG)

sr <- 1000
n <- 8000

# One-sided band power (local, matching the package's PSD convention).
band_power <- function(sig, sr, flo, fhi) {
  nn <- length(sig)
  psd <- (Mod(fft(sig - mean(sig)))^2) / nn
  freqs <- seq(0, sr / 2, length.out = nn %/% 2 + 1)
  psd <- psd[seq_len(length(freqs))]
  sum(psd[freqs >= flo & freqs <= fhi])
}

# Regular, realistic QRS train (wide enough that its energy sits below ~40 Hz).
qrs_train <- function(n, sr, bpm = 60, amp = 1.5) {
  q <- numeric(n)
  per <- as.integer(sr * 60 / bpm)
  tw <- -40:40
  for (p in seq(per, n - 60, by = per)) q[p + tw] <- amp * exp(-tw^2 / (2 * 10^2))
  q
}

pe1 <- function(sig) {
  PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)), samplingRate = sr)
}

test_that("emgQualityCheck flags injected ECG and reports the expected metrics", {
  set.seed(1)
  emg <- rnorm(n, sd = 0.15)
  contam <- emg + qrs_train(n, sr, bpm = 60)
  qc <- emgQualityCheck(pe1(contam))

  expect_s3_class(qc, "data.frame")
  expect_true(all(c("channel", "snr_db", "baseline_noise", "powerline_50",
                    "powerline_60", "powerline_ratio", "clipping_pct",
                    "ecg_score", "ecg_contamination") %in% names(qc)))
  expect_true(qc$ecg_contamination[1])
  expect_gt(qc$ecg_score[1], 0.5)
})

test_that("clean EMG is not flagged and passes the QC gate", {
  set.seed(2)
  clean <- rnorm(n, sd = 0.02)
  clean[3000:5000] <- rnorm(2001, sd = 0.5)          # a real contraction
  qc <- emgQualityCheck(pe1(clean))
  gate <- emgQCgate(pe1(clean))

  expect_false(qc$ecg_contamination[1])
  expect_lt(qc$ecg_score[1], 0.5)
  expect_true(gate$pass)
  expect_equal(gate$channels$reasons[1], "")
})

test_that("emgRemoveECG (highpass) cuts QRS-band energy while keeping EMG band", {
  set.seed(3)
  emg <- rnorm(n, sd = 0.2)
  contam <- emg + qrs_train(n, sr, bpm = 75)
  pe <- pe1(contam)

  cleaned <- emgRemoveECG(pe, method = "highpass")
  expect_s4_class(cleaned, "PhysioExperiment")
  expect_true("ecg_removed" %in% SummarizedExperiment::assayNames(cleaned))
  cl <- SummarizedExperiment::assay(cleaned, "ecg_removed")[, 1]

  qrs_reduction <- 1 - band_power(cl, sr, 5, 30) / band_power(contam, sr, 5, 30)
  emg_retention <- band_power(cl, sr, 20, 450) / band_power(contam, sr, 20, 450)
  expect_gt(qrs_reduction, 0.70)                     # >70% QRS-band removed
  expect_gt(emg_retention, 0.80)                     # >80% EMG-band retained

  # contamination is no longer flagged after removal
  cleaned_pe <- pe1(cl)
  expect_false(emgQualityCheck(cleaned_pe)$ecg_contamination[1])
})

test_that("template and gating removal run and reduce QRS-band energy", {
  set.seed(4)
  contam <- rnorm(n, sd = 0.2) + qrs_train(n, sr, bpm = 60)
  pe <- pe1(contam)
  for (m in c("template", "gating")) {
    cl <- SummarizedExperiment::assay(emgRemoveECG(pe, method = m),
                                      "ecg_removed")[, 1]
    qrs_reduction <- 1 - band_power(cl, sr, 5, 30) / band_power(contam, sr, 5, 30)
    emg_retention <- band_power(cl, sr, 20, 450) / band_power(contam, sr, 20, 450)
    expect_gt(qrs_reduction, 0.50)
    expect_gt(emg_retention, 0.80)                   # in-band EMG preserved
  }
})

test_that("power-line metric peaks at the injected line frequency", {
  for (f in c(50, 60)) {
    set.seed(f)
    tone <- rnorm(n, sd = 0.2) + 0.3 * sin(2 * pi * f * seq_len(n) / sr)
    qc <- emgQualityCheck(pe1(tone))
    if (f == 50) {
      expect_gt(qc$powerline_50[1], qc$powerline_60[1])
    } else {
      expect_gt(qc$powerline_60[1], qc$powerline_50[1])
    }
    expect_gt(qc$powerline_ratio[1], 0.1)
  }
})

test_that("emgQCgate reports failure reasons and handles multi-channel", {
  set.seed(6)
  contam <- rnorm(n, sd = 0.15) + qrs_train(n, sr, bpm = 60)   # ECG-contaminated
  clean <- rnorm(n, sd = 0.02); clean[3000:5000] <- rnorm(2001, sd = 0.5)
  pe <- PhysioExperiment(
    assays = list(raw = cbind(contam, clean)),
    colData = S4Vectors::DataFrame(label = c("bad", "good"),
                                   type = c("EMG", "EMG")),
    samplingRate = sr)

  gate <- emgQCgate(pe)
  expect_false(gate$pass)                            # one channel fails
  expect_equal(nrow(gate$channels), 2L)
  expect_true(grepl("ECG", gate$channels$reasons[1]))
  expect_true(gate$channels$pass[2])
})
