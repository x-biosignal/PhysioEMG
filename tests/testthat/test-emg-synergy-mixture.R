make_two_subgroups <- function(M = 6, r = 2, per = 6, T_each = 60, seed = 1) {
  set.seed(seed)
  LA <- matrix(abs(rnorm(M * r)), M, r)
  LB <- matrix(abs(rnorm(M * r)), M, r)
  gen <- function(L) t(L %*% matrix(rnorm(r * T_each), r, T_each) +
                       matrix(rnorm(M * T_each, 0, 0.3), M, T_each))
  data <- c(replicate(per, gen(LA), simplify = FALSE),
            replicate(per, gen(LB), simplify = FALSE))
  list(data = data, truth = rep(1:2, each = per), LA = LA, LB = LB)
}

# clustering accuracy up to label permutation (2 clusters)
accuracy2 <- function(pred, truth) {
  a <- mean(pred == truth); max(a, 1 - a)
}

test_that("recovers two subgroups with distinct synergy structure", {
  d <- make_two_subgroups()
  fit <- muscleSynergyMixture(d$data, n_clusters = 2, n_factors = 2, seed = 1)

  expect_s3_class(fit, "synergy_mixture")
  expect_length(fit$synergies, 2)
  expect_equal(dim(fit$synergies[[1]]), c(6, 2))
  expect_length(fit$cluster, 12)
  expect_gt(accuracy2(fit$cluster, d$truth), 0.9)     # subgroups recovered
})

test_that("BIC prefers 2 clusters over 1 for two-subgroup data", {
  d <- make_two_subgroups()
  bic1 <- muscleSynergyMixture(d$data, 1, 2, seed = 1)$bic
  bic2 <- muscleSynergyMixture(d$data, 2, 2, seed = 1)$bic
  expect_lt(bic2, bic1)
})

test_that("responsibilities are a valid partition of unity", {
  d <- make_two_subgroups()
  fit <- muscleSynergyMixture(d$data, 2, 2, seed = 1)
  expect_equal(rowSums(fit$responsibilities), rep(1, 12), tolerance = 1e-8)
  expect_true(all(fit$responsibilities >= 0))
})

test_that("MPCA method also recovers the subgroups", {
  d <- make_two_subgroups()
  fit <- muscleSynergyMixture(d$data, 2, 2, method = "mpca", seed = 1)
  expect_equal(fit$method, "mpca")
  expect_gt(accuracy2(fit$cluster, d$truth), 0.9)
  # MPCA uses isotropic noise: each cluster's uniqueness is constant
  expect_equal(length(unique(round(fit$uniqueness[[1]], 8))), 1)
})

test_that("input validation", {
  expect_error(muscleSynergyMixture(list(matrix(1, 5, 3)), 2, 1),
               "at least")
  bad <- list(matrix(1, 5, 3), matrix(1, 5, 4))
  expect_error(muscleSynergyMixture(bad, 2, 1), "same number of muscles")
})

test_that("PhysioExperiment bridge extracts assays and fits the mixture", {
  skip_if_not_installed("PhysioCore")
  d <- make_two_subgroups()
  pes <- lapply(d$data, function(m)
    PhysioCore::PhysioExperiment(assays = list(raw = m), samplingRate = 100))
  fit <- physioSynergyMixture(pes, n_clusters = 2, n_factors = 2, seed = 1)
  expect_s3_class(fit, "synergy_mixture")
  expect_length(fit$cluster, 12)
  expect_gt(accuracy2(fit$cluster, d$truth), 0.9)
  # also accepts raw matrices
  fit2 <- physioSynergyMixture(d$data, 2, 2, seed = 1)
  expect_equal(fit2$cluster, fit$cluster)
})
