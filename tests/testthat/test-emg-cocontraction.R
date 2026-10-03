library(testthat)
library(PhysioEMG)

# Build a two-channel envelope PhysioExperiment from two vectors.
env_pe <- function(a, b, sr = 1000) {
  PhysioExperiment(
    assays = list(env = cbind(a, b)),
    colData = S4Vectors::DataFrame(label = c("AG", "ANT"),
                                   type = c("EMG", "EMG")),
    samplingRate = sr)
}

test_that("identical envelopes give Falconer-Winter CCI = 100%", {
  t <- seq_len(1000)
  a <- exp(-((t - 500) / 100)^2) + 0.05
  pe <- env_pe(a, a)

  res <- emgCoContraction(pe, agonist = 1, antagonist = 2,
                          method = "falconer_winter")
  expect_type(res, "list")
  expect_true(all(c("method", "summary", "timeseries") %in% names(res)))
  expect_equal(res$summary, 100, tolerance = 1e-8)
})

test_that("non-overlapping bursts give Falconer-Winter CCI ~ 0", {
  n <- 2000
  a <- numeric(n); b <- numeric(n)
  a[200:600] <- 1                              # agonist burst
  b[1200:1600] <- 1                            # antagonist burst (disjoint)
  pe <- env_pe(a, b)

  res <- emgCoContraction(pe, 1, 2, method = "falconer_winter")
  expect_lt(res$summary, 1e-6)                 # ~ 0

  # Frost (temporal overlap ratio) is also ~ 0 for disjoint bursts
  resf <- emgCoContraction(pe, 1, 2, method = "frost")
  expect_lt(resf$summary, 1e-6)
})

test_that("Rudolph index is symmetric under agonist/antagonist swap", {
  set.seed(3)
  a <- abs(sin(seq(0, 6, length.out = 1500))) * 5 + runif(1500, 0, 0.5)
  b <- abs(cos(seq(0, 6, length.out = 1500))) * 3 + runif(1500, 0, 0.5)
  pe <- env_pe(a, b)

  r_ab <- emgCoContraction(pe, 1, 2, method = "rudolph",
                           normalize = FALSE)$summary
  r_ba <- emgCoContraction(pe, 2, 1, method = "rudolph",
                           normalize = FALSE)$summary
  expect_equal(r_ab, r_ba, tolerance = 1e-10)
})

test_that("Rudolph index is bounded in [0, max(env)^2]", {
  set.seed(4)
  a <- abs(sin(seq(0, 8, length.out = 2000))) * 8 + runif(2000, 0, 1)
  b <- abs(cos(seq(0, 8, length.out = 2000))) * 6 + runif(2000, 0, 1)
  pe <- env_pe(a, b)
  max_env <- max(c(abs(a), abs(b)))

  res <- emgCoContraction(pe, 1, 2, method = "rudolph", normalize = FALSE,
                          window_sec = 0.2)
  expect_gte(res$summary, 0)
  expect_lte(res$summary, max_env^2)
  # every windowed value obeys the same bound
  expect_true(all(res$timeseries$cci >= 0))
  expect_true(all(res$timeseries$cci <= max_env^2))
})

test_that("time-resolved output and channel labels work", {
  t <- seq_len(3000)
  a <- exp(-((t - 1500) / 400)^2)
  b <- exp(-((t - 1500) / 400)^2) * 0.7
  pe <- env_pe(a, b)

  res <- emgCoContraction(pe, agonist = "AG", antagonist = "ANT",
                          method = "falconer_winter", window_sec = 0.5)
  expect_s3_class(res$timeseries, "data.frame")
  expect_true(all(c("window", "time_sec", "cci") %in% names(res$timeseries)))
  expect_equal(nrow(res$timeseries), 6L)       # 3 s / 0.5 s
  expect_true(all(res$timeseries$cci >= 0 & res$timeseries$cci <= 100))

  expect_error(emgCoContraction(pe, "NOPE", "ANT"), "not found")
  expect_error(emgCoContraction(pe, 1, 5), "out of range")
})
