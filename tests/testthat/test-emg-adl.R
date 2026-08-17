# EMG muscle-activation summary of an ADL task (emgADLActivation).

make_pair <- function(ag_scale = 1, an_scale = 0.2, seed = 1) {
  set.seed(seed)
  n <- 1000; t <- seq_len(n)
  burst <- exp(-((t - 500) / 150)^2)
  ag <- ag_scale * burst + abs(rnorm(n, 0, 0.01))
  an <- an_scale * burst + abs(rnorm(n, 0, 0.01))
  PhysioCore::PhysioExperiment(
    assays = S4Vectors::SimpleList(envelope = cbind(BIC = ag, TRI = an)),
    colData = S4Vectors::DataFrame(label = c("BIC", "TRI"), type = "EMG"),
    samplingRate = 1000)
}

test_that("emgADLActivation summarises co-contraction and amplitude", {
  pe <- make_pair(an_scale = 0.2)
  r <- emgADLActivation(pe, "BIC", "TRI", task = "drinking")
  expect_s3_class(r, "emg_adl_activation")
  expect_equal(r$icf_code, "d560")
  expect_true(is.finite(r$cocontraction_index))
  expect_gt(r$agonist_peak, r$antagonist_peak)         # agonist dominates
  expect_equal(r$unit, "envelope")
  expect_gt(r$agonist_active_frac, 0)
})

test_that("more antagonist activity raises the co-contraction index", {
  low <- emgADLActivation(make_pair(an_scale = 0.1), "BIC", "TRI")
  high <- emgADLActivation(make_pair(an_scale = 0.9), "BIC", "TRI")
  expect_gt(high$cocontraction_index, low$cocontraction_index)
})

test_that("MVC reference yields percent-MVC amplitudes", {
  pe <- make_pair(ag_scale = 1)
  mvc <- make_pair(ag_scale = 2, an_scale = 2, seed = 9)   # stronger max effort
  r <- emgADLActivation(pe, "BIC", "TRI", mvc_data = mvc, task = "feeding")
  expect_equal(r$unit, "%MVC")
  expect_lt(r$agonist_peak, 100)                       # task < maximal effort
  expect_gt(r$agonist_peak, 0)
})

test_that("channel labels and indices both resolve", {
  pe <- make_pair()
  by_label <- emgADLActivation(pe, "BIC", "TRI")
  by_index <- emgADLActivation(pe, 1, 2)
  expect_equal(by_label$agonist_peak, by_index$agonist_peak)
  expect_error(emgADLActivation(pe, "NOPE", "TRI"), "not found")
})
