library(testthat)
library(PhysioEMG)

test_that("emgCoherenceNetwork returns symmetric bounded network", {
  set.seed(1)
  pe <- make_emg(n_time = 2048, n_channels = 4, sr = 1000)

  result <- emgCoherenceNetwork(
    pe,
    freq_band = c(20, 150),
    nperseg = 256,
    aggregate = "mean",
    threshold = 0.2
  )

  expect_type(result, "list")
  expect_true(all(c("network", "coherence", "frequencies", "channel_names") %in% names(result)))
  expect_equal(dim(result$network), c(4, 4))
  expect_equal(result$network, t(result$network), tolerance = 1e-10)
  expect_equal(unname(diag(result$network)), rep(1, 4))
  expect_true(all(result$network >= 0 & result$network <= 1, na.rm = TRUE))

  expect_equal(dim(result$coherence)[2:3], c(4, 4))
  expect_true(is.logical(result$adjacency))
  expect_equal(dim(result$adjacency), c(4, 4))
  expect_false(any(diag(result$adjacency)))
})

test_that("emgDynamicWaveletNetwork returns valid dynamic network series", {
  set.seed(2)
  pe <- make_emg(n_time = 1500, n_channels = 3, sr = 500)

  result <- emgDynamicWaveletNetwork(
    pe,
    frequencies = seq(10, 100, by = 10),
    freq_band = c(20, 80),
    window_sec = 0.4,
    step_sec = 0.2,
    n_cycles = 6,
    smoothing_cycles = 2,
    threshold = 0.3
  )

  n_windows <- length(result$window_times)
  expect_equal(dim(result$network), c(n_windows, 3, 3))
  expect_true(all(result$network >= 0 & result$network <= 1, na.rm = TRUE))
  for (w in seq_len(n_windows)) {
    expect_equal(unname(diag(result$network[w, , ])), rep(1, 3))
  }

  expect_equal(dim(result$static_summary), c(3, 3))
  expect_true(is.logical(result$adjacency))
  expect_equal(dim(result$adjacency), c(n_windows, 3, 3))
  for (w in seq_len(n_windows)) {
    expect_false(any(diag(result$adjacency[w, , ])))
  }
})

test_that("emgInterpretNetworkKG annotates top edges with metadata and KG links", {
  net <- matrix(
    c(
      1.0, 0.8, 0.3,
      0.8, 1.0, 0.6,
      0.3, 0.6, 1.0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("EMG1", "EMG2", "EMG3"), c("EMG1", "EMG2", "EMG3"))
  )

  node_metadata <- data.frame(
    channel = c("EMG1", "EMG2", "EMG3"),
    kg_node = c("M1", "M2", "M3"),
    muscle_group = c("flexor", "extensor", "flexor"),
    stringsAsFactors = FALSE
  )
  kg_edges <- data.frame(
    node_a = c("M1", "M2"),
    node_b = c("M2", "M3"),
    relation = c("co_activate", "synergist"),
    stringsAsFactors = FALSE
  )

  result <- emgInterpretNetworkKG(
    network = net,
    node_metadata = node_metadata,
    kg_edges = kg_edges,
    threshold = 0.5,
    top_n = 3
  )

  expect_type(result, "list")
  expect_true(nrow(result$edge_table) >= 2)
  expect_true(all(c("source", "target", "weight", "kg_linked") %in% names(result$edge_table)))
  expect_true(all(result$edge_table$weight >= 0.5))
  expect_true(any(result$edge_table$kg_linked))
  expect_type(result$summary, "list")
  expect_equal(result$summary$n_channels, 3)

  dyn <- array(NA_real_, dim = c(2, 3, 3))
  dyn[1, , ] <- net
  dyn[2, , ] <- net * 0.9
  diag(dyn[2, , ]) <- 1
  dyn_res <- emgInterpretNetworkKG(dyn, threshold = 0.45, top_n = 1, window = 2)

  expect_equal(dim(dyn_res$network_matrix), c(3, 3))
  expect_equal(nrow(dyn_res$edge_table), 1)
})

test_that("emgPartialCoherenceNetwork returns bounded symmetric matrix", {
  set.seed(3)
  pe <- make_emg(n_time = 2048, n_channels = 4, sr = 1000)

  result <- emgPartialCoherenceNetwork(
    pe,
    freq_band = c(20, 200),
    nperseg = 256,
    aggregate = "mean",
    threshold = 0.2
  )

  expect_type(result, "list")
  expect_true(all(c("network", "partial_coherence", "frequencies") %in% names(result)))
  expect_equal(dim(result$network), c(4, 4))
  expect_equal(result$network, t(result$network), tolerance = 1e-10)
  expect_equal(unname(diag(result$network)), rep(1, 4))
  expect_true(all(result$network >= 0 & result$network <= 1, na.rm = TRUE))
  expect_equal(dim(result$partial_coherence)[2:3], c(4, 4))
  expect_true(is.logical(result$adjacency))
})

test_that("emgWPLINetwork returns bounded symmetric matrix", {
  set.seed(4)
  pe <- make_emg(n_time = 2048, n_channels = 4, sr = 1000)

  result <- emgWPLINetwork(
    pe,
    freq_band = c(20, 200),
    nperseg = 256,
    aggregate = "mean",
    threshold = 0.1
  )

  expect_type(result, "list")
  expect_true(all(c("network", "wpli", "frequencies") %in% names(result)))
  expect_equal(dim(result$network), c(4, 4))
  expect_equal(result$network, t(result$network), tolerance = 1e-10)
  expect_equal(unname(diag(result$network)), rep(1, 4))
  expect_true(all(result$network >= 0 & result$network <= 1, na.rm = TRUE))
  expect_equal(dim(result$wpli)[2:3], c(4, 4))
})

test_that("emgDirectedGCNetwork detects directional coupling", {
  set.seed(5)
  n <- 4000
  sr <- 1000
  x <- as.numeric(arima.sim(model = list(ar = 0.7), n = n, sd = 1))
  y <- numeric(n)
  y[1] <- rnorm(1)
  for (t in 2:n) {
    y[t] <- 0.6 * y[t - 1] + 0.5 * x[t - 1] + rnorm(1, sd = 0.8)
  }
  z <- as.numeric(arima.sim(model = list(ar = 0.5), n = n, sd = 1))

  mat <- cbind(X = x, Y = y, Z = z)
  pe <- PhysioCore::PhysioExperiment(
    assays = list(raw = mat),
    colData = S4Vectors::DataFrame(
      label = c("X", "Y", "Z"),
      type = rep("EMG", 3),
      row.names = c("X", "Y", "Z")
    ),
    samplingRate = sr
  )

  result <- emgDirectedGCNetwork(
    pe,
    max_lag = 8,
    score = "f_stat",
    p_value_cutoff = 0.05
  )

  expect_type(result, "list")
  expect_equal(dim(result$network), c(3, 3))
  expect_equal(dim(result$p_values), c(3, 3))
  expect_equal(unname(diag(result$network)), rep(0, 3))
  expect_true(is.logical(result$adjacency))

  expect_gt(result$network["X", "Y"], result$network["Y", "X"])
  expect_lt(result$p_values["X", "Y"], 0.05)
})

test_that("emgCoordinationStructure recovers modular coordination pattern", {
  net <- matrix(
    c(
      0, 0.9, 0.8, 0.1, 0.05, 0.05,
      0.9, 0, 0.85, 0.1, 0.05, 0.05,
      0.8, 0.85, 0, 0.12, 0.08, 0.05,
      0.1, 0.1, 0.12, 0, 0.88, 0.82,
      0.05, 0.05, 0.08, 0.88, 0, 0.9,
      0.05, 0.05, 0.05, 0.82, 0.9, 0
    ),
    nrow = 6,
    byrow = TRUE,
    dimnames = list(paste0("M", 1:6), paste0("M", 1:6))
  )

  res <- emgCoordinationStructure(net, threshold = 0.1, max_modules = 4)
  expect_type(res, "list")
  expect_true(all(c("node_metrics", "summary", "modules", "network") %in% names(res)))
  expect_equal(nrow(res$node_metrics), 6)
  expect_true(res$summary$n_modules >= 2)
  expect_true(is.finite(res$summary$modularity))
  expect_gt(res$summary$modularity, 0)
  expect_true(all(res$node_metrics$participation >= 0))
})

test_that("emgCoordinationStructure accepts directed/asymmetric network input", {
  dnet <- matrix(
    c(
      0, 2.0, 0.2,
      0.5, 0, 1.3,
      1.1, 0.4, 0
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
  )

  res <- emgCoordinationStructure(
    list(network = dnet),
    directed = TRUE,
    symmetrize = "max",
    normalize = TRUE
  )

  expect_equal(dim(res$network), c(3, 3))
  expect_equal(res$network, t(res$network), tolerance = 1e-12)
  expect_equal(unname(diag(res$network)), rep(0, 3))
  expect_equal(length(res$modules), 3)
})
