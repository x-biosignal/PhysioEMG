test_that("emgEnvelope RMS (prefix-sum) matches the naive O(n*window) reference", {
  set.seed(7)
  x <- make_emg(n_time = 5000, n_channels = 1, sr = 1000)
  env <- SummarizedExperiment::assay(
    emgEnvelope(x, method = "rms", window_ms = 50), "envelope"
  )
  sig <- SummarizedExperiment::assay(x, defaultAssay(x))[, 1]
  sig_sq <- sig^2
  ws <- max(1L, as.integer(round(50 / 1000 * 1000)))  # window in samples
  hw <- ws %/% 2
  ref <- vapply(seq_along(sig), function(i) {
    lo <- max(1L, i - hw); hi <- min(length(sig), i + hw)
    sqrt(sum(sig_sq[lo:hi]) / (hi - lo + 1))
  }, numeric(1))

  expect_equal(as.numeric(env[, 1]), ref, tolerance = 1e-10)
})
