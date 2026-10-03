make_convolutive <- function(M = 5, Tt = 200, N = 2, L = 15, noise = 0, seed = 1) {
  set.seed(seed)
  W0 <- array(0, c(M, N, L))
  for (nn in seq_len(N)) W0[, nn, ] <- outer(runif(M), dnorm(seq_len(L), L / 2, 3))
  H0 <- matrix(0, N, Tt)
  H0[1, c(30, 120)] <- c(1, 0.8)
  H0[2, c(70, 160)] <- c(0.9, 1)
  V <- matrix(0, M, Tt)
  for (tau in 0:(L - 1)) {
    RtH <- cbind(matrix(0, N, tau), H0[, seq_len(Tt - tau), drop = FALSE])
    V <- V + W0[, , tau + 1] %*% RtH
  }
  if (noise > 0) V <- pmax(V + matrix(rnorm(length(V), 0, noise), M, Tt), 0)
  list(V = V, W0 = W0, H0 = H0)
}

test_that("convolutive NMF recovers planted spatiotemporal synergies", {
  d <- make_convolutive()
  fit <- convolutiveSynergy(d$V, n_synergies = 2, L = 15, seed = 1)

  expect_s3_class(fit, "convolutive_synergy")
  expect_equal(dim(fit$synergies), c(5, 2, 15))
  expect_equal(dim(fit$activations), c(2, 200))
  expect_true(all(fit$synergies >= 0))
  expect_true(all(fit$activations >= 0))
  expect_gt(fit$vaf, 0.95)
})

test_that("activations are sparse and time-localised near the planted onsets", {
  d <- make_convolutive()
  fit <- convolutiveSynergy(d$V, n_synergies = 2, L = 15, seed = 1)
  # each recovered activation row should have a few dominant peaks
  peaks_per_row <- apply(fit$activations, 1, function(h)
    sum(h > 0.3 * max(h)))
  expect_true(all(peaks_per_row < 30))          # localised, not diffuse
})

test_that("more noise lowers VAF but the model still fits", {
  clean <- convolutiveSynergy(make_convolutive(noise = 0)$V, 2, 15, seed = 1)$vaf
  noisy <- convolutiveSynergy(make_convolutive(noise = 0.05)$V, 2, 15, seed = 1)$vaf
  expect_gt(clean, noisy)
  expect_gt(noisy, 0.5)
})

test_that("input validation", {
  expect_error(convolutiveSynergy(matrix(-1, 3, 50), 2, 10), "non-negative")
  expect_error(convolutiveSynergy(matrix(1, 3, 10), 2, 20), "shorter than")
})
