library(testthat)
library(PhysioEMG)

test_that("muscleSynergy with NMF decomposes correctly", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)

  result <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

  expect_type(result, "list")
  expect_true(all(c("W", "H", "vaf", "method") %in% names(result)))
  expect_equal(nrow(result$W), 3)
  expect_equal(ncol(result$W), 8)
  expect_equal(nrow(result$H), 1000)
  expect_equal(ncol(result$H), 3)
  expect_gt(result$vaf, 0.5)
})

test_that("muscleSynergy with PCA works", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)

  result <- muscleSynergy(pe, n_synergies = 3, method = "pca")

  expect_type(result, "list")
  expect_equal(nrow(result$W), 3)
  expect_equal(ncol(result$W), 8)
})

test_that("synergyReconstruct reconstructs data", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)
  result <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

  recon <- synergyReconstruct(result, n_synergies = 2)

  expect_type(recon, "list")
  expect_true("reconstructed" %in% names(recon))
  expect_equal(dim(recon$reconstructed), c(1000, 8))
  expect_lt(recon$vaf, result$vaf)
})

test_that("synergyCompare computes similarity", {
  set.seed(42)
  pe <- make_emg_synergy(n_time = 1000, n_channels = 8, n_synergies = 3)
  r1 <- muscleSynergy(pe, n_synergies = 3, method = "nmf")
  set.seed(123)
  r2 <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

  sim <- synergyCompare(r1, r2)

  expect_s3_class(sim, "data.frame")
  expect_true("correlation" %in% names(sim))
  expect_true(all(abs(sim$correlation) > 0.3))
})

# Synthetic EMG built from K well-separated synergies (each drives a disjoint
# muscle group), with equal per-synergy energy so VAF(N) ~ N/K and there is a
# sharp VAF jump at N = K.
clean_synergy_pe <- function(K, n_time = 600, mper = 2, noise = 0.01, seed = 1) {
  set.seed(seed)
  M <- K * mper
  H <- matrix(0, n_time, K)
  W <- matrix(0, K, M)
  for (k in seq_len(K)) {
    c0 <- n_time * k / (K + 1)
    h <- exp(-((seq_len(n_time) - c0)^2) / (2 * (n_time / 12)^2))
    H[, k] <- h / sqrt(sum(h^2))                       # unit energy
    w <- runif(mper, 0.6, 1.4)
    W[k, ((k - 1) * mper + 1):(k * mper)] <- w / sqrt(sum(w^2))
  }
  data <- H %*% W + matrix(abs(rnorm(n_time * M, sd = noise)), n_time, M)
  PhysioExperiment(assays = list(env = data), samplingRate = 1000)
}

# ---- multi-restart muscleSynergy -------------------------------------------

test_that("multi-restart NMF VAF >= single-restart and is reproducible", {
  pe <- make_emg_synergy(n_time = 800, n_channels = 8, n_synergies = 3)

  single <- muscleSynergy(pe, n_synergies = 3, method = "nmf",
                          n_restarts = 1, seed = 42)
  multi <- muscleSynergy(pe, n_synergies = 3, method = "nmf",
                         n_restarts = 12, seed = 42)

  expect_gte(multi$vaf, single$vaf - 1e-12)     # never worse
  expect_equal(length(multi$all_vaf), 12L)
  expect_true("convergence" %in% names(multi))
  expect_true(all(diff(cummax(multi$all_vaf)) >= 0))

  # best-of-restarts reproducible under a fixed seed
  multi2 <- muscleSynergy(pe, n_synergies = 3, method = "nmf",
                          n_restarts = 12, seed = 42)
  expect_equal(multi$vaf, multi2$vaf)
  expect_equal(multi$W, multi2$W)
})

# ---- muscleSynergyOrder -----------------------------------------------------

test_that("muscleSynergyOrder selects N == K on K-synergy data", {
  for (K in 2:4) {
    pe <- clean_synergy_pe(K, seed = K)
    ord <- muscleSynergyOrder(pe, max_synergies = 6, n_restarts = 8, seed = 1)

    expect_true(all(c("vaf_curve", "selected", "selection_rule",
                      "per_channel_vaf") %in% names(ord)))
    expect_equal(ord$selected, K)
    expect_equal(ord$selection_rule, "vaf_threshold")
    # global VAF at the selected order clears the 0.90 threshold
    expect_gte(ord$vaf_curve$vaf[ord$selected], 0.90)
    # a genuine jump: order K-1 is below threshold
    expect_lt(ord$vaf_curve$vaf[K - 1], 0.90)
  }
})

test_that("muscleSynergyOrder is reproducible and reports per-channel VAF", {
  pe <- clean_synergy_pe(3, seed = 5)
  o1 <- muscleSynergyOrder(pe, max_synergies = 5, n_restarts = 6, seed = 9)
  o2 <- muscleSynergyOrder(pe, max_synergies = 5, n_restarts = 6, seed = 9)

  expect_equal(o1$selected, o2$selected)
  expect_equal(o1$vaf_curve$vaf, o2$vaf_curve$vaf)
  expect_length(o1$per_channel_vaf, 6L)         # 3 synergies x 2 muscles
  expect_true(o1$per_synergy_vaf_met)
})
