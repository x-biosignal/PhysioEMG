make_synergy_tensor <- function(M = 6, Tt = 50, K = 8, R = 2, noise = 0,
                                 seed = 1) {
  set.seed(seed)
  A0 <- matrix(runif(M * R), M, R)
  B0 <- matrix(runif(Tt * R), Tt, R)
  C0 <- matrix(runif(K * R), K, R)
  X <- array(0, c(M, Tt, K))
  for (r in seq_len(R)) for (k in seq_len(K))
    X[, , k] <- X[, , k] + (A0[, r] %o% B0[, r]) * C0[k, r]
  if (noise > 0) X <- pmax(X + array(rnorm(length(X), 0, noise), dim(X)), 0)
  list(X = X, A0 = A0, B0 = B0, C0 = C0)
}

# greedy best-match correlation between two factor matrices (permutation-free)
match_corr <- function(A, A0) {
  R <- ncol(A0); used <- logical(R); cors <- numeric(R)
  for (r in seq_len(R)) {
    cc <- sapply(seq_len(R), function(j)
      if (used[j]) -Inf else abs(stats::cor(A[, r], A0[, j])))
    j <- which.max(cc); used[j] <- TRUE; cors[r] <- cc[j]
  }
  cors
}

test_that("non-negative PARAFAC recovers planted synergies", {
  d <- make_synergy_tensor(R = 2)
  fit <- muscleSynergyTensor(d$X, n_synergies = 2, restarts = 5, seed = 1)

  expect_s3_class(fit, "synergy_tensor")
  expect_equal(dim(fit$muscle_weights), c(6, 2))
  expect_equal(dim(fit$temporal), c(50, 2))
  expect_equal(dim(fit$trial_loadings), c(8, 2))
  expect_true(all(fit$muscle_weights >= 0))
  expect_gt(fit$vaf, 0.95)
  # recovered muscle weights match the planted ones (up to permutation)
  expect_gt(min(match_corr(fit$muscle_weights, d$A0)), 0.9)
})

test_that("VAF increases with rank and order selection works", {
  d <- make_synergy_tensor(R = 3, noise = 0.01)
  sel <- muscleSynergyTensorOrder(d$X, max_synergies = 4, vaf_threshold = 0.90,
                                  restarts = 3, seed = 1)
  expect_true(all(diff(sel$vaf) >= -1e-6))       # non-decreasing VAF
  expect_lte(sel$n_synergies, 4)
  expect_gte(sel$vaf[sel$n_synergies], 0.90)
  expect_s3_class(sel$fit, "synergy_tensor")
})

test_that("input validation", {
  expect_error(muscleSynergyTensor(matrix(1, 3, 3), 2), "3-way array")
  X <- make_synergy_tensor()$X; X[1] <- -1
  expect_error(muscleSynergyTensor(X, 2), "non-negative")
})
