# Space-by-time (sNM3F) and non-negative Tucker synergies, verified on synthetic
# data built from known factors.

test_that("space-by-time recovers a low-rank space-by-time structure", {
  set.seed(1)
  M <- 8; Tn <- 50; N <- 12
  Ws <- matrix(runif(M * 2), M, 2); Wt <- matrix(runif(Tn * 2), Tn, 2)
  trials <- lapply(seq_len(N), function(n) Ws %*% matrix(runif(4), 2, 2) %*% t(Wt))
  fit <- spaceByTimeSynergy(trials, n_spatial = 2, n_temporal = 2,
                            n_restart = 3, seed = 1)
  expect_s3_class(fit, "spacetime_synergy")
  expect_gt(fit$vaf, 0.98)                            # exact low-rank -> near-perfect
  expect_equal(dim(fit$spatial), c(M, 2))
  expect_equal(dim(fit$temporal), c(Tn, 2))
  expect_length(fit$activation, N)
  expect_true(all(fit$spatial >= 0) && all(fit$temporal >= 0))
})

test_that("space-by-time VAF increases with more synergies", {
  set.seed(2)
  M <- 8; Tn <- 40; N <- 10
  Ws <- matrix(runif(M * 3), M, 3); Wt <- matrix(runif(Tn * 3), Tn, 3)
  trials <- lapply(seq_len(N), function(n) Ws %*% matrix(runif(9), 3, 3) %*% t(Wt))
  v1 <- spaceByTimeSynergy(trials, 1, 1, n_restart = 2, seed = 3)$vaf
  v3 <- spaceByTimeSynergy(trials, 3, 3, n_restart = 2, seed = 3)$vaf
  expect_gt(v3, v1)
  expect_gt(v3, 0.97)
})

test_that("space-by-time rejects negative EMG", {
  trials <- lapply(1:3, function(n) matrix(rnorm(6 * 20), 6, 20))
  expect_error(spaceByTimeSynergy(trials, 2, 2), "non-negative")
})

test_that("non-negative Tucker reconstructs a known Tucker tensor", {
  set.seed(4)
  A <- matrix(runif(6 * 2), 6, 2); B <- matrix(runif(40 * 2), 40, 2)
  C <- matrix(runif(8 * 2), 8, 2); G <- array(runif(8), c(2, 2, 2))
  X <- reconstructTucker(G, A, B, C)
  expect_equal(dim(X), c(6, 40, 8))
  fit <- tuckerSynergy(X, ranks = c(2, 2, 2), n_restart = 3, seed = 1)
  expect_s3_class(fit, "tucker_synergy")
  expect_gt(fit$vaf, 0.98)                            # recovers the structure
  expect_equal(dim(fit$core), c(2, 2, 2))
  expect_true(all(fit$spatial >= 0))
})

test_that("Tucker with a non-diagonal core beats a diagonal-core (PARAFAC-like) fit", {
  set.seed(5)
  # a core with strong OFF-diagonal coupling: PARAFAC (diagonal) cannot capture it
  A <- matrix(runif(6 * 2), 6, 2); B <- matrix(runif(30 * 2), 30, 2)
  C <- matrix(runif(10 * 2), 10, 2)
  G <- array(0, c(2, 2, 2)); G[1, 2, 1] <- 1; G[2, 1, 2] <- 1   # purely off-diagonal
  X <- reconstructTucker(G, A, B, C)
  fit <- tuckerSynergy(X, ranks = c(2, 2, 2), n_restart = 4, seed = 2)
  expect_gt(fit$vaf, 0.95)                            # Tucker captures the interaction
  expect_output(print(fit), "Tucker")
})

test_that("Tucker accepts a list of muscle x time matrices", {
  set.seed(6)
  A <- matrix(runif(5 * 2), 5, 2); B <- matrix(runif(25 * 2), 25, 2)
  C <- matrix(runif(6 * 2), 6, 2); G <- array(runif(8), c(2, 2, 2))
  X <- reconstructTucker(G, A, B, C)
  trials <- lapply(seq_len(dim(X)[3]), function(k) matrix(X[, , k], 5, 25))
  fit <- tuckerSynergy(trials, ranks = c(2, 2, 2), n_restart = 2, seed = 3)
  expect_gt(fit$vaf, 0.97)
})
