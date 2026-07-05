#' EMG Coherence Network Analysis
#'
#' Builds a static muscle coordination network from pairwise magnitude-squared
#' coherence between EMG channels.
#'
#' @param x A PhysioExperiment object.
#' @param freq_band Optional numeric vector \code{c(low, high)} in Hz.
#'   If NULL, uses all frequencies.
#' @param channels Integer vector of channel indices to include. If NULL, uses all.
#' @param nperseg Segment length for Welch estimation (default: 256).
#' @param noverlap Overlap length (default: \code{floor(nperseg / 2)}).
#' @param assay_name Input assay name. If NULL, uses default assay.
#' @param aggregate Aggregation across frequency bins: "mean", "max", or "median".
#' @param threshold Optional threshold for binary adjacency matrix.
#' @return A list with:
#'   \describe{
#'     \item{network}{Numeric matrix (channel x channel) of coherence strength.}
#'     \item{adjacency}{Logical matrix after thresholding, or NULL.}
#'     \item{coherence}{3D array (freq x channel x channel).}
#'     \item{frequencies}{Frequency vector (Hz) corresponding to \code{coherence}.}
#'     \item{channel_names}{Channel labels used in the network.}
#'   }
#' @seealso [emgDynamicWaveletNetwork()] for time-varying networks,
#'   [emgInterpretNetworkKG()] for annotation-aware interpretation.
#' @export
emgCoherenceNetwork <- function(x, freq_band = NULL, channels = NULL,
                                nperseg = 256L, noverlap = NULL,
                                assay_name = NULL,
                                aggregate = c("mean", "max", "median"),
                                threshold = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  aggregate <- match.arg(aggregate)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required.", call. = FALSE)
  }

  data <- .emg_extract_matrix(x, channels = channels, assay_name = assay_name)
  signal_data <- data$signal_data
  ch_names <- data$channel_names

  n_time <- nrow(signal_data)
  n_channels <- ncol(signal_data)
  nperseg <- as.integer(nperseg)
  if (nperseg < 8L) {
    stop("nperseg must be >= 8.", call. = FALSE)
  }
  if (is.null(noverlap)) {
    noverlap <- as.integer(floor(nperseg / 2L))
  } else {
    noverlap <- as.integer(noverlap)
  }
  if (noverlap < 0L || noverlap >= nperseg) {
    stop("noverlap must satisfy 0 <= noverlap < nperseg.", call. = FALSE)
  }
  if (n_time < 8L) {
    stop("Signal length is too short.", call. = FALSE)
  }

  n_freqs <- floor(min(nperseg, n_time) / 2L) + 1L
  freqs_full <- seq(0, sr / 2, length.out = n_freqs)
  freq_idx <- seq_len(n_freqs)
  if (!is.null(freq_band)) {
    if (!is.numeric(freq_band) || length(freq_band) != 2L || freq_band[1] >= freq_band[2]) {
      stop("freq_band must be numeric c(low, high) with low < high.", call. = FALSE)
    }
    freq_idx <- which(freqs_full >= freq_band[1] & freqs_full <= freq_band[2])
    if (length(freq_idx) == 0L) {
      stop("No frequency bins in requested freq_band.", call. = FALSE)
    }
  }
  freqs <- freqs_full[freq_idx]

  psd_list <- vector("list", n_channels)
  for (i in seq_len(n_channels)) {
    psd_list[[i]] <- .emg_welch_psd(signal_data[, i], nperseg, noverlap, sr)
  }

  coh_array <- array(
    NA_real_,
    dim = c(length(freq_idx), n_channels, n_channels),
    dimnames = list(
      sprintf("%.4f", freqs),
      ch_names,
      ch_names
    )
  )

  for (i in seq_len(n_channels)) {
    coh_array[, i, i] <- 1
    if (i == n_channels) next
    for (j in (i + 1L):n_channels) {
      csd_ij <- .emg_welch_csd(signal_data[, i], signal_data[, j], nperseg, noverlap, sr)
      denom <- psd_list[[i]]$psd * psd_list[[j]]$psd + .Machine$double.eps
      coh_full <- (Mod(csd_ij$csd)^2) / denom
      coh <- pmin(pmax(coh_full[freq_idx], 0), 1)
      coh_array[, i, j] <- coh
      coh_array[, j, i] <- coh
    }
  }

  agg_fun <- switch(
    aggregate,
    mean = function(v) mean(v, na.rm = TRUE),
    max = function(v) max(v, na.rm = TRUE),
    median = function(v) stats::median(v, na.rm = TRUE)
  )
  network <- apply(coh_array, c(2, 3), agg_fun)
  diag(network) <- 1

  adjacency <- NULL
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
      stop("threshold must be a single numeric value.", call. = FALSE)
    }
    adjacency <- network >= threshold
    diag(adjacency) <- FALSE
  }

  list(
    network = network,
    adjacency = adjacency,
    coherence = coh_array,
    frequencies = freqs,
    channel_names = ch_names,
    method = "welch_magnitude_squared_coherence",
    aggregate = aggregate
  )
}

#' EMG Partial Coherence Network Analysis
#'
#' Builds a static network from pairwise partial coherence, estimated from
#' the inverse cross-spectral density matrix at each frequency.
#'
#' @param x A PhysioExperiment object.
#' @param freq_band Optional numeric vector \code{c(low, high)} in Hz.
#' @param channels Integer vector of channel indices to include. If NULL, uses all.
#' @param nperseg Segment length for Welch estimation (default: 256).
#' @param noverlap Overlap length (default: \code{floor(nperseg / 2)}).
#' @param assay_name Input assay name. If NULL, uses default assay.
#' @param aggregate Aggregation across frequency bins: "mean", "max", or "median".
#' @param threshold Optional threshold for binary adjacency matrix.
#' @param ridge Ridge regularization added to spectral matrix inversion.
#' @return A list with:
#'   \describe{
#'     \item{network}{Numeric matrix (channel x channel) of partial coherence.}
#'     \item{adjacency}{Logical matrix after thresholding, or NULL.}
#'     \item{partial_coherence}{3D array (freq x channel x channel).}
#'     \item{frequencies}{Frequency vector (Hz).}
#'     \item{channel_names}{Channel labels used in the network.}
#'   }
#' @seealso [emgCoherenceNetwork()], [emgWPLINetwork()]
#' @export
emgPartialCoherenceNetwork <- function(x, freq_band = NULL, channels = NULL,
                                       nperseg = 256L, noverlap = NULL,
                                       assay_name = NULL,
                                       aggregate = c("mean", "max", "median"),
                                       threshold = NULL,
                                       ridge = 1e-6) {
  stopifnot(inherits(x, "PhysioExperiment"))
  aggregate <- match.arg(aggregate)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required.", call. = FALSE)
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("ridge must be a non-negative scalar.", call. = FALSE)
  }

  data <- .emg_extract_matrix(x, channels = channels, assay_name = assay_name)
  signal_data <- data$signal_data
  ch_names <- data$channel_names

  n_time <- nrow(signal_data)
  n_channels <- ncol(signal_data)
  nperseg <- as.integer(nperseg)
  if (nperseg < 8L) {
    stop("nperseg must be >= 8.", call. = FALSE)
  }
  if (is.null(noverlap)) {
    noverlap <- as.integer(floor(nperseg / 2L))
  } else {
    noverlap <- as.integer(noverlap)
  }
  if (noverlap < 0L || noverlap >= nperseg) {
    stop("noverlap must satisfy 0 <= noverlap < nperseg.", call. = FALSE)
  }
  if (n_time < 8L) {
    stop("Signal length is too short.", call. = FALSE)
  }

  fft_seg <- .emg_welch_fft_segments(signal_data, nperseg = nperseg, noverlap = noverlap, sr = sr)
  freqs_full <- fft_seg$frequencies
  freq_idx <- seq_along(freqs_full)
  if (!is.null(freq_band)) {
    if (!is.numeric(freq_band) || length(freq_band) != 2L || freq_band[1] >= freq_band[2]) {
      stop("freq_band must be numeric c(low, high) with low < high.", call. = FALSE)
    }
    freq_idx <- which(freqs_full >= freq_band[1] & freqs_full <= freq_band[2])
    if (length(freq_idx) == 0L) {
      stop("No frequency bins in requested freq_band.", call. = FALSE)
    }
  }
  freqs <- freqs_full[freq_idx]

  pc_array <- array(
    NA_real_,
    dim = c(length(freq_idx), n_channels, n_channels),
    dimnames = list(sprintf("%.4f", freqs), ch_names, ch_names)
  )

  for (fi_local in seq_along(freq_idx)) {
    fi <- freq_idx[[fi_local]]
    fmat <- fft_seg$fft[fi, , , drop = TRUE]
    # Frequency slice spectral matrix (channels x channels)
    s_mat <- (fmat %*% Conj(t(fmat))) / ncol(fmat)
    if (ridge > 0) {
      s_mat <- s_mat + diag(ridge, n_channels)
    }
    p_mat <- .emg_safe_solve(s_mat)
    d <- Re(diag(p_mat))
    den <- outer(d, d, "*") + .Machine$double.eps
    pc <- (Mod(p_mat)^2) / den
    pc <- pmin(pmax(Re(pc), 0), 1)
    diag(pc) <- 1
    pc_array[fi_local, , ] <- pc
  }

  agg_fun <- switch(
    aggregate,
    mean = function(v) mean(v, na.rm = TRUE),
    max = function(v) max(v, na.rm = TRUE),
    median = function(v) stats::median(v, na.rm = TRUE)
  )
  network <- apply(pc_array, c(2, 3), agg_fun)
  diag(network) <- 1

  adjacency <- NULL
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
      stop("threshold must be a single numeric value.", call. = FALSE)
    }
    adjacency <- network >= threshold
    diag(adjacency) <- FALSE
  }

  list(
    network = network,
    adjacency = adjacency,
    partial_coherence = pc_array,
    frequencies = freqs,
    channel_names = ch_names,
    method = "welch_partial_coherence",
    aggregate = aggregate
  )
}

#' EMG Weighted Phase-Lag Index (wPLI) Network
#'
#' Builds a static network from pairwise weighted phase-lag index (wPLI),
#' which is less sensitive to zero-lag coupling artifacts.
#'
#' @param x A PhysioExperiment object.
#' @param freq_band Optional numeric vector \code{c(low, high)} in Hz.
#' @param channels Integer vector of channel indices to include. If NULL, uses all.
#' @param nperseg Segment length for Welch estimation (default: 256).
#' @param noverlap Overlap length (default: \code{floor(nperseg / 2)}).
#' @param assay_name Input assay name. If NULL, uses default assay.
#' @param aggregate Aggregation across frequency bins: "mean", "max", or "median".
#' @param threshold Optional threshold for binary adjacency matrix.
#' @param debiased Logical; if TRUE, uses debiased wPLI estimator.
#' @return A list with:
#'   \describe{
#'     \item{network}{Numeric matrix (channel x channel) of wPLI values.}
#'     \item{adjacency}{Logical matrix after thresholding, or NULL.}
#'     \item{wpli}{3D array (freq x channel x channel).}
#'     \item{frequencies}{Frequency vector (Hz).}
#'     \item{channel_names}{Channel labels used in the network.}
#'   }
#' @seealso [emgCoherenceNetwork()], [emgPartialCoherenceNetwork()]
#' @export
emgWPLINetwork <- function(x, freq_band = NULL, channels = NULL,
                           nperseg = 256L, noverlap = NULL,
                           assay_name = NULL,
                           aggregate = c("mean", "max", "median"),
                           threshold = NULL,
                           debiased = FALSE) {
  stopifnot(inherits(x, "PhysioExperiment"))
  aggregate <- match.arg(aggregate)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required.", call. = FALSE)
  }

  data <- .emg_extract_matrix(x, channels = channels, assay_name = assay_name)
  signal_data <- data$signal_data
  ch_names <- data$channel_names

  n_time <- nrow(signal_data)
  n_channels <- ncol(signal_data)
  nperseg <- as.integer(nperseg)
  if (nperseg < 8L) {
    stop("nperseg must be >= 8.", call. = FALSE)
  }
  if (is.null(noverlap)) {
    noverlap <- as.integer(floor(nperseg / 2L))
  } else {
    noverlap <- as.integer(noverlap)
  }
  if (noverlap < 0L || noverlap >= nperseg) {
    stop("noverlap must satisfy 0 <= noverlap < nperseg.", call. = FALSE)
  }
  if (n_time < 8L) {
    stop("Signal length is too short.", call. = FALSE)
  }

  fft_seg <- .emg_welch_fft_segments(signal_data, nperseg = nperseg, noverlap = noverlap, sr = sr)
  freqs_full <- fft_seg$frequencies
  freq_idx <- seq_along(freqs_full)
  if (!is.null(freq_band)) {
    if (!is.numeric(freq_band) || length(freq_band) != 2L || freq_band[1] >= freq_band[2]) {
      stop("freq_band must be numeric c(low, high) with low < high.", call. = FALSE)
    }
    freq_idx <- which(freqs_full >= freq_band[1] & freqs_full <= freq_band[2])
    if (length(freq_idx) == 0L) {
      stop("No frequency bins in requested freq_band.", call. = FALSE)
    }
  }
  freqs <- freqs_full[freq_idx]

  wpli_array <- array(
    NA_real_,
    dim = c(length(freq_idx), n_channels, n_channels),
    dimnames = list(sprintf("%.4f", freqs), ch_names, ch_names)
  )

  for (i in seq_len(n_channels)) {
    wpli_array[, i, i] <- 1
    if (i == n_channels) next
    for (j in (i + 1L):n_channels) {
      vals <- rep(NA_real_, length(freq_idx))
      for (fi_local in seq_along(freq_idx)) {
        fi <- freq_idx[[fi_local]]
        x_fft <- fft_seg$fft[fi, i, ]
        y_fft <- fft_seg$fft[fi, j, ]
        im_vals <- Im(x_fft * Conj(y_fft))

        if (isTRUE(debiased)) {
          num <- (sum(im_vals, na.rm = TRUE)^2) - sum(im_vals^2, na.rm = TRUE)
          den <- (sum(abs(im_vals), na.rm = TRUE)^2) - sum(im_vals^2, na.rm = TRUE)
          if (!is.finite(den) || den <= .Machine$double.eps) {
            w <- 0
          } else {
            w <- num / den
          }
        } else {
          den <- sum(abs(im_vals), na.rm = TRUE)
          if (!is.finite(den) || den <= .Machine$double.eps) {
            w <- 0
          } else {
            w <- abs(sum(im_vals, na.rm = TRUE)) / den
          }
        }
        vals[[fi_local]] <- pmin(pmax(w, 0), 1)
      }
      wpli_array[, i, j] <- vals
      wpli_array[, j, i] <- vals
    }
  }

  agg_fun <- switch(
    aggregate,
    mean = function(v) mean(v, na.rm = TRUE),
    max = function(v) max(v, na.rm = TRUE),
    median = function(v) stats::median(v, na.rm = TRUE)
  )
  network <- apply(wpli_array, c(2, 3), agg_fun)
  diag(network) <- 1

  adjacency <- NULL
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
      stop("threshold must be a single numeric value.", call. = FALSE)
    }
    adjacency <- network >= threshold
    diag(adjacency) <- FALSE
  }

  list(
    network = network,
    adjacency = adjacency,
    wpli = wpli_array,
    frequencies = freqs,
    channel_names = ch_names,
    method = if (isTRUE(debiased)) "debiased_wpli" else "wpli",
    aggregate = aggregate
  )
}

#' Directed EMG Network via Pairwise Granger Causality
#'
#' Estimates a directed muscle coordination network using pairwise Granger
#' causality in the time domain.
#'
#' @param x A PhysioExperiment object.
#' @param channels Integer vector of channel indices to include. If NULL, uses all.
#' @param assay_name Input assay name. If NULL, uses default assay.
#' @param max_lag Lag order (in samples) for autoregressive modeling.
#' @param score Directed edge metric: "f_stat" or "delta_r2".
#' @param threshold Optional threshold for adjacency based on selected score.
#' @param p_value_cutoff Optional p-value threshold for adjacency.
#' @param standardize Logical; if TRUE, z-score each channel before GC.
#' @return A list with:
#'   \describe{
#'     \item{network}{Directed numeric matrix (source x target).}
#'     \item{p_values}{Directed matrix of GC p-values.}
#'     \item{adjacency}{Logical directed matrix, or NULL.}
#'     \item{channel_names}{Channel labels used in the network.}
#'     \item{lag}{Lag order used for modeling.}
#'   }
#' @export
emgDirectedGCNetwork <- function(x, channels = NULL, assay_name = NULL,
                                 max_lag = 10L,
                                 score = c("f_stat", "delta_r2"),
                                 threshold = NULL,
                                 p_value_cutoff = NULL,
                                 standardize = TRUE) {
  stopifnot(inherits(x, "PhysioExperiment"))
  score <- match.arg(score)

  max_lag <- as.integer(max_lag)
  if (is.na(max_lag) || max_lag < 1L) {
    stop("max_lag must be >= 1.", call. = FALSE)
  }
  if (!is.null(p_value_cutoff)) {
    if (!is.numeric(p_value_cutoff) || length(p_value_cutoff) != 1L ||
      is.na(p_value_cutoff) || p_value_cutoff <= 0 || p_value_cutoff >= 1) {
      stop("p_value_cutoff must be in (0, 1).", call. = FALSE)
    }
  }

  data <- .emg_extract_matrix(x, channels = channels, assay_name = assay_name)
  signal_data <- data$signal_data
  ch_names <- data$channel_names

  n_time <- nrow(signal_data)
  n_channels <- ncol(signal_data)
  if (n_time <= (4L * max_lag + 10L)) {
    stop("Signal length is too short for selected max_lag.", call. = FALSE)
  }

  if (isTRUE(standardize)) {
    for (j in seq_len(n_channels)) {
      s <- stats::sd(signal_data[, j], na.rm = TRUE)
      if (is.finite(s) && s > .Machine$double.eps) {
        signal_data[, j] <- (signal_data[, j] - mean(signal_data[, j], na.rm = TRUE)) / s
      } else {
        signal_data[, j] <- signal_data[, j] - mean(signal_data[, j], na.rm = TRUE)
      }
    }
  }

  net <- matrix(0, nrow = n_channels, ncol = n_channels, dimnames = list(ch_names, ch_names))
  pmat <- matrix(NA_real_, nrow = n_channels, ncol = n_channels, dimnames = list(ch_names, ch_names))

  for (i in seq_len(n_channels)) {
    for (j in seq_len(n_channels)) {
      if (i == j) next
      gc <- .emg_pairwise_granger(signal_data[, i], signal_data[, j], lag_order = max_lag)
      net[i, j] <- if (score == "f_stat") gc$f_stat else gc$delta_r2
      pmat[i, j] <- gc$p_value
    }
  }

  adjacency <- NULL
  if (!is.null(threshold) || !is.null(p_value_cutoff)) {
    cond <- matrix(TRUE, nrow = n_channels, ncol = n_channels)
    if (!is.null(threshold)) {
      if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
        stop("threshold must be a single numeric value.", call. = FALSE)
      }
      cond <- cond & (net >= threshold)
    }
    if (!is.null(p_value_cutoff)) {
      cond <- cond & (pmat <= p_value_cutoff)
    }
    diag(cond) <- FALSE
    adjacency <- cond
  }

  list(
    network = net,
    p_values = pmat,
    adjacency = adjacency,
    channel_names = ch_names,
    lag = max_lag,
    score = score,
    method = "pairwise_granger_causality"
  )
}

#' EMG Coordination Structure Summary from Network Topology
#'
#' Quantifies higher-order muscle coordination structure from a weighted
#' network matrix, including module structure, efficiency, and node roles.
#'
#' @param network Numeric square matrix (channels x channels) or a list
#'   containing \code{$network}.
#' @param threshold Optional edge-weight threshold. Values below threshold are
#'   set to zero before topology estimation.
#' @param n_modules Optional number of modules. If NULL, chooses a value
#'   automatically by maximizing weighted modularity over candidates.
#' @param max_modules Maximum number of candidate modules for automatic search.
#' @param directed Logical; set TRUE if \code{network} is directed/asymmetric.
#' @param symmetrize Method to convert directed matrices to undirected form:
#'   "mean", "max", or "min".
#' @param normalize Logical; if TRUE, rescales weights to \code{[0, 1]}.
#' @return A list with:
#'   \describe{
#'     \item{network}{Processed undirected weighted network matrix.}
#'     \item{node_metrics}{Data.frame of node-level topology features.}
#'     \item{modules}{Named integer vector of module assignments.}
#'     \item{summary}{Global network topology summary.}
#'   }
#' @seealso [emgCoherenceNetwork()], [emgDynamicWaveletNetwork()]
#' @export
emgCoordinationStructure <- function(network,
                                     threshold = NULL,
                                     n_modules = NULL,
                                     max_modules = 6L,
                                     directed = FALSE,
                                     symmetrize = c("mean", "max", "min"),
                                     normalize = TRUE) {
  symmetrize <- match.arg(symmetrize)

  if (is.list(network) && "network" %in% names(network)) {
    net <- network$network
  } else {
    net <- network
  }
  if (!is.matrix(net) || !is.numeric(net) || nrow(net) != ncol(net)) {
    stop("network must be a numeric square matrix, or a list containing $network.", call. = FALSE)
  }
  if (nrow(net) < 2L) {
    stop("network must contain at least 2 channels.", call. = FALSE)
  }
  if (!all(is.finite(net) | is.na(net))) {
    stop("network contains non-finite values.", call. = FALSE)
  }

  n <- nrow(net)
  ch_names <- rownames(net)
  if (is.null(ch_names)) {
    ch_names <- paste0("EMG", seq_len(n))
  }

  net <- as.matrix(net)
  net[is.na(net)] <- 0
  if (isTRUE(directed) || !isTRUE(all.equal(net, t(net)))) {
    net <- switch(
      symmetrize,
      mean = (net + t(net)) / 2,
      max = pmax(net, t(net)),
      min = pmin(net, t(net))
    )
  }
  net <- pmax(net, 0)
  diag(net) <- 0

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) || threshold < 0) {
      stop("threshold must be a non-negative numeric scalar.", call. = FALSE)
    }
    net[net < threshold] <- 0
  }

  if (isTRUE(normalize)) {
    mx <- max(net, na.rm = TRUE)
    if (is.finite(mx) && mx > .Machine$double.eps) {
      net <- net / mx
    }
  }

  strengths <- rowSums(net)
  max_strength <- max(strengths, na.rm = TRUE)
  strength_norm <- if (is.finite(max_strength) && max_strength > 0) strengths / max_strength else rep(0, n)
  clustering <- .emg_weighted_clustering(net)

  dist_full <- .emg_shortest_path_lengths(net)
  eff_vals <- 1 / dist_full
  eff_vals[!is.finite(eff_vals)] <- NA_real_
  diag(eff_vals) <- NA_real_
  global_eff <- mean(eff_vals, na.rm = TRUE)
  if (!is.finite(global_eff)) global_eff <- 0

  local_eff <- numeric(n)
  for (i in seq_len(n)) {
    nei <- which(net[i, ] > 0)
    if (length(nei) < 2L) {
      local_eff[i] <- 0
      next
    }
    sub_net <- net[nei, nei, drop = FALSE]
    d_sub <- .emg_shortest_path_lengths(sub_net)
    e_sub <- 1 / d_sub
    e_sub[!is.finite(e_sub)] <- NA_real_
    diag(e_sub) <- NA_real_
    v <- mean(e_sub, na.rm = TRUE)
    local_eff[i] <- if (is.finite(v)) v else 0
  }

  modules <- rep(1L, n)
  names(modules) <- ch_names
  modularity <- NA_real_

  max_modules <- as.integer(max_modules)
  if (is.na(max_modules) || max_modules < 2L) {
    max_modules <- min(6L, n - 1L)
  }
  if (!is.null(n_modules)) {
    n_modules <- as.integer(n_modules)
    if (is.na(n_modules) || n_modules < 1L) {
      stop("n_modules must be >= 1 when provided.", call. = FALSE)
    }
  }

  if (sum(net) > .Machine$double.eps && n >= 3L) {
    if (is.null(n_modules)) {
      k_max <- min(max_modules, n - 1L)
      if (k_max >= 2L) {
        modules <- .emg_detect_modules(net, k_candidates = 2L:k_max)
      }
    } else if (n_modules >= 2L && n_modules <= n) {
      modules <- .emg_detect_modules(net, k_candidates = n_modules)
    }
    modularity <- .emg_weighted_modularity(net, modules)
  }

  uniq_mod <- sort(unique(as.integer(modules)))
  participation <- numeric(n)
  within_z <- numeric(n)
  for (i in seq_len(n)) {
    ki <- strengths[[i]]
    if (!is.finite(ki) || ki <= .Machine$double.eps) {
      participation[i] <- 0
      next
    }
    frac_sq <- 0
    for (m in uniq_mod) {
      idx <- which(modules == m)
      kis <- sum(net[i, idx], na.rm = TRUE)
      frac_sq <- frac_sq + (kis / ki)^2
    }
    participation[i] <- max(0, 1 - frac_sq)
  }
  for (m in uniq_mod) {
    idx <- which(modules == m)
    if (length(idx) == 0L) next
    kin <- rowSums(net[idx, idx, drop = FALSE])
    s <- stats::sd(kin, na.rm = TRUE)
    if (!is.finite(s) || s <= .Machine$double.eps) {
      within_z[idx] <- 0
    } else {
      within_z[idx] <- (kin - mean(kin, na.rm = TRUE)) / s
    }
  }

  node_metrics <- data.frame(
    channel = ch_names,
    strength = as.numeric(strengths),
    normalized_strength = as.numeric(strength_norm),
    clustering = as.numeric(clustering),
    local_efficiency = as.numeric(local_eff),
    participation = as.numeric(participation),
    module = as.integer(modules),
    within_module_z = as.numeric(within_z),
    stringsAsFactors = FALSE
  )

  summary <- list(
    n_nodes = n,
    n_edges = sum(net[upper.tri(net)] > 0, na.rm = TRUE),
    density = mean(net[upper.tri(net)] > 0, na.rm = TRUE),
    mean_edge_weight = mean(net[upper.tri(net)], na.rm = TRUE),
    global_efficiency = as.numeric(global_eff),
    mean_clustering = mean(clustering, na.rm = TRUE),
    modularity = as.numeric(modularity),
    n_modules = length(unique(modules))
  )

  list(
    network = net,
    node_metrics = node_metrics,
    modules = modules,
    summary = summary,
    method = "coordination_topology_summary"
  )
}

#' Dynamic Wavelet Coherence Network for EMG
#'
#' Builds a time-varying coordination network by computing wavelet coherence
#' between channel pairs and aggregating coherence within sliding windows.
#'
#' @param x A PhysioExperiment object.
#' @param frequencies Numeric vector of wavelet center frequencies (Hz).
#' @param freq_band Optional numeric vector \code{c(low, high)} for aggregation.
#' @param channels Integer vector of channel indices to include. If NULL, uses all.
#' @param window_sec Sliding window length in seconds.
#' @param step_sec Sliding window step in seconds.
#' @param n_cycles Number of Morlet cycles (default: 7).
#' @param smoothing_cycles Smoothing width in cycles (default: 3).
#' @param assay_name Input assay name. If NULL, uses default assay.
#' @param aggregate Aggregation across time-frequency bins: "mean", "max", or "median".
#' @param threshold Optional threshold for binary adjacency network per window.
#' @param respect_coi Logical; if TRUE, masks frequencies below COI before aggregation.
#' @return A list with:
#'   \describe{
#'     \item{network}{3D array (window x channel x channel).}
#'     \item{adjacency}{Logical 3D array thresholded from \code{network}, or NULL.}
#'     \item{window_times}{Window center times in seconds.}
#'     \item{static_summary}{Mean network across windows.}
#'     \item{frequencies}{Frequency vector used for wavelet transform.}
#'     \item{coi}{Cone-of-influence frequency at each time sample.}
#'   }
#' @seealso [emgCoherenceNetwork()] for static spectral network,
#'   [emgInterpretNetworkKG()] for annotation-aware interpretation.
#' @export
emgDynamicWaveletNetwork <- function(x,
                                     frequencies = seq(5, 120, by = 5),
                                     freq_band = NULL,
                                     channels = NULL,
                                     window_sec = 0.5,
                                     step_sec = 0.1,
                                     n_cycles = 7,
                                     smoothing_cycles = 3,
                                     assay_name = NULL,
                                     aggregate = c("mean", "max", "median"),
                                     threshold = NULL,
                                     respect_coi = TRUE) {
  stopifnot(inherits(x, "PhysioExperiment"))
  aggregate <- match.arg(aggregate)

  sr <- samplingRate(x)
  if (is.na(sr) || sr <= 0) {
    stop("Valid sampling rate is required.", call. = FALSE)
  }
  if (!is.numeric(frequencies) || length(frequencies) < 2L || any(!is.finite(frequencies)) || any(frequencies <= 0)) {
    stop("frequencies must be numeric values > 0.", call. = FALSE)
  }
  frequencies <- sort(unique(as.numeric(frequencies)))

  if (!is.numeric(window_sec) || length(window_sec) != 1L || window_sec <= 0) {
    stop("window_sec must be > 0.", call. = FALSE)
  }
  if (!is.numeric(step_sec) || length(step_sec) != 1L || step_sec <= 0) {
    stop("step_sec must be > 0.", call. = FALSE)
  }
  if (!is.numeric(n_cycles) || length(n_cycles) != 1L || n_cycles <= 0) {
    stop("n_cycles must be > 0.", call. = FALSE)
  }
  if (!is.numeric(smoothing_cycles) || length(smoothing_cycles) != 1L || smoothing_cycles <= 0) {
    stop("smoothing_cycles must be > 0.", call. = FALSE)
  }

  data <- .emg_extract_matrix(x, channels = channels, assay_name = assay_name)
  signal_data <- data$signal_data
  ch_names <- data$channel_names
  n_time <- nrow(signal_data)
  n_channels <- ncol(signal_data)

  window_samples <- as.integer(round(window_sec * sr))
  step_samples <- as.integer(round(step_sec * sr))
  if (window_samples < 2L || window_samples > n_time) {
    stop("window_sec produced an invalid window length.", call. = FALSE)
  }
  if (step_samples < 1L) {
    stop("step_sec produced an invalid step length.", call. = FALSE)
  }

  starts <- seq.int(1L, n_time - window_samples + 1L, by = step_samples)
  n_windows <- length(starts)
  window_times <- ((starts - 1L) + (window_samples / 2)) / sr

  freq_idx <- seq_along(frequencies)
  if (!is.null(freq_band)) {
    if (!is.numeric(freq_band) || length(freq_band) != 2L || freq_band[1] >= freq_band[2]) {
      stop("freq_band must be numeric c(low, high) with low < high.", call. = FALSE)
    }
    freq_idx <- which(frequencies >= freq_band[1] & frequencies <= freq_band[2])
    if (length(freq_idx) == 0L) {
      stop("No frequencies in requested freq_band.", call. = FALSE)
    }
  }

  transforms <- vector("list", n_channels)
  for (i in seq_len(n_channels)) {
    transforms[[i]] <- .emg_morlet_wavelet_transform(signal_data[, i], frequencies, n_cycles, sr)
  }

  coi <- .emg_compute_coi(n_time, sr, n_cycles)

  agg_fun <- switch(
    aggregate,
    mean = function(v) mean(v, na.rm = TRUE),
    max = function(v) max(v, na.rm = TRUE),
    median = function(v) stats::median(v, na.rm = TRUE)
  )

  network <- array(
    NA_real_,
    dim = c(n_windows, n_channels, n_channels),
    dimnames = list(
      sprintf("%.4f", window_times),
      ch_names,
      ch_names
    )
  )

  for (i in seq_len(n_channels)) {
    network[, i, i] <- 1
    if (i == n_channels) next
    for (j in (i + 1L):n_channels) {
      w_x <- transforms[[i]]
      w_y <- transforms[[j]]
      w_xy <- w_x * Conj(w_y)
      p_xx <- Mod(w_x)^2
      p_yy <- Mod(w_y)^2

      coh_tf <- matrix(NA_real_, nrow = n_time, ncol = length(frequencies))

      for (fi in seq_along(frequencies)) {
        f0 <- frequencies[fi]
        sigma_t <- smoothing_cycles / f0 * sr
        half_width <- min(ceiling(3 * sigma_t), floor(n_time / 2))
        kernel_idx <- seq(-half_width, half_width)
        kernel <- exp(-0.5 * (kernel_idx / sigma_t)^2)
        kernel <- kernel / sum(kernel)

        smooth_xy <- .emg_convolve_mirror(w_xy[, fi], kernel)
        smooth_xx <- .emg_convolve_mirror(p_xx[, fi], kernel)
        smooth_yy <- .emg_convolve_mirror(p_yy[, fi], kernel)
        denom <- smooth_xx * smooth_yy
        valid <- denom > .Machine$double.eps

        coh_i <- rep(NA_real_, n_time)
        coh_i[valid] <- (Mod(smooth_xy[valid])^2) / denom[valid]
        coh_i <- pmin(pmax(coh_i, 0), 1)

        if (isTRUE(respect_coi)) {
          coh_i[frequencies[fi] < coi] <- NA_real_
        }
        coh_tf[, fi] <- coh_i
      }

      series_ij <- rep(NA_real_, n_windows)
      for (w in seq_len(n_windows)) {
        idx_w <- starts[w]:(starts[w] + window_samples - 1L)
        vals <- coh_tf[idx_w, freq_idx, drop = FALSE]
        if (all(is.na(vals))) {
          series_ij[w] <- NA_real_
        } else {
          series_ij[w] <- agg_fun(vals)
        }
      }

      network[, i, j] <- series_ij
      network[, j, i] <- series_ij
    }
  }

  adjacency <- NULL
  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
      stop("threshold must be a single numeric value.", call. = FALSE)
    }
    adjacency <- network >= threshold
    for (w in seq_len(n_windows)) {
      diag(adjacency[w, , ]) <- FALSE
    }
  }

  static_summary <- apply(network, c(2, 3), function(v) mean(v, na.rm = TRUE))
  diag(static_summary) <- 1

  list(
    network = network,
    adjacency = adjacency,
    window_times = window_times,
    static_summary = static_summary,
    frequencies = frequencies,
    freq_indices = freq_idx,
    coi = coi,
    channel_names = ch_names,
    method = "wavelet_coherence_dynamic_network",
    aggregate = aggregate
  )
}

#' Interpret EMG Network with Knowledge-Graph Metadata
#'
#' Adds metadata/KG context to high-weight network edges. This function does
#' not require a specific backend and works with user-provided tables exported
#' from PhysioAnnotationHub/physioKG workflows.
#'
#' @param network A square matrix (channels x channels) or 3D array
#'   (window x channels x channels).
#' @param node_metadata Optional data.frame containing at least `channel`.
#'   Additional columns (e.g. `kg_node`, `muscle_name`, `muscle_group`) are
#'   carried into the edge table.
#' @param kg_edges Optional data.frame of KG links. Expected columns are
#'   `node_a`, `node_b`, and optional `relation`.
#' @param threshold Optional edge threshold. Default is 75th percentile of
#'   upper-triangle weights.
#' @param top_n Maximum number of edges returned after thresholding.
#' @param window Optional window index when `network` is a 3D array.
#'   If NULL, uses the mean across windows.
#' @return A list with:
#'   \describe{
#'     \item{edge_table}{Ranked edge table with optional metadata/KG annotations.}
#'     \item{threshold}{Applied threshold value.}
#'     \item{network_matrix}{Matrix used for interpretation.}
#'     \item{summary}{List of summary statistics.}
#'     \item{kg_relation_summary}{Relation counts from matched KG links, or NULL.}
#'   }
#' @export
emgInterpretNetworkKG <- function(network,
                                  node_metadata = NULL,
                                  kg_edges = NULL,
                                  threshold = NULL,
                                  top_n = 20L,
                                  window = NULL) {
  if (length(dim(network)) == 3L) {
    if (is.null(window)) {
      net <- apply(network, c(2, 3), function(v) mean(v, na.rm = TRUE))
    } else {
      window <- as.integer(window)
      if (!is.finite(window) || window < 1L || window > dim(network)[1]) {
        stop("window is out of range.", call. = FALSE)
      }
      net <- network[window, , , drop = TRUE]
    }
  } else if (is.matrix(network)) {
    net <- network
  } else {
    stop("network must be a matrix or 3D array.", call. = FALSE)
  }

  if (!is.numeric(net) || nrow(net) != ncol(net)) {
    stop("network matrix must be square numeric.", call. = FALSE)
  }

  ch_names <- rownames(net)
  if (is.null(ch_names)) {
    ch_names <- paste0("Ch", seq_len(nrow(net)))
    rownames(net) <- ch_names
  }
  if (is.null(colnames(net))) {
    colnames(net) <- ch_names
  }

  ut <- which(upper.tri(net), arr.ind = TRUE)
  edge_table <- data.frame(
    source = rownames(net)[ut[, 1]],
    target = colnames(net)[ut[, 2]],
    weight = net[ut],
    stringsAsFactors = FALSE
  )
  edge_table <- edge_table[order(edge_table$weight, decreasing = TRUE), , drop = FALSE]
  rownames(edge_table) <- NULL

  if (is.null(threshold)) {
    threshold <- stats::quantile(edge_table$weight, probs = 0.75, na.rm = TRUE, names = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("threshold must be a single numeric value.", call. = FALSE)
  }

  edge_table <- edge_table[edge_table$weight >= threshold, , drop = FALSE]
  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1L) {
    stop("top_n must be >= 1.", call. = FALSE)
  }
  if (nrow(edge_table) > top_n) {
    edge_table <- edge_table[seq_len(top_n), , drop = FALSE]
  }

  if (!is.null(node_metadata)) {
    if (!is.data.frame(node_metadata)) {
      stop("node_metadata must be a data.frame.", call. = FALSE)
    }
    nm <- node_metadata
    if (!("channel" %in% names(nm))) {
      if (!is.null(rownames(nm))) {
        nm$channel <- rownames(nm)
      } else {
        stop("node_metadata must contain `channel` column.", call. = FALSE)
      }
    }

    nm_source <- nm
    names(nm_source) <- paste0("source_", names(nm_source))
    nm_target <- nm
    names(nm_target) <- paste0("target_", names(nm_target))

    edge_table <- merge(edge_table, nm_source,
                        by.x = "source", by.y = "source_channel",
                        all.x = TRUE, sort = FALSE)
    edge_table <- merge(edge_table, nm_target,
                        by.x = "target", by.y = "target_channel",
                        all.x = TRUE, sort = FALSE)
    edge_table <- edge_table[order(edge_table$weight, decreasing = TRUE), , drop = FALSE]
    rownames(edge_table) <- NULL
  }

  kg_relation_summary <- NULL
  if (!is.null(kg_edges) && nrow(edge_table) > 0) {
    if (!is.data.frame(kg_edges)) {
      stop("kg_edges must be a data.frame.", call. = FALSE)
    }
    kge <- kg_edges
    if (!all(c("node_a", "node_b") %in% names(kge))) {
      if (all(c("source", "target") %in% names(kge))) {
        kge <- data.frame(
          node_a = as.character(kge$source),
          node_b = as.character(kge$target),
          relation = if ("relation" %in% names(kge)) as.character(kge$relation) else "linked",
          stringsAsFactors = FALSE
        )
      } else {
        stop("kg_edges must contain node_a/node_b (or source/target).", call. = FALSE)
      }
    }
    if (!("relation" %in% names(kge))) {
      kge$relation <- "linked"
    }

    if (all(c("source_kg_node", "target_kg_node") %in% names(edge_table))) {
      edge_table$kg_linked <- FALSE
      edge_table$kg_relations <- NA_character_
      edge_table$kg_n_links <- 0L

      for (i in seq_len(nrow(edge_table))) {
        a <- edge_table$source_kg_node[i]
        b <- edge_table$target_kg_node[i]
        if (is.na(a) || is.na(b)) next
        hit <- (kge$node_a == a & kge$node_b == b) | (kge$node_a == b & kge$node_b == a)
        if (any(hit)) {
          rels <- unique(as.character(kge$relation[hit]))
          edge_table$kg_linked[i] <- TRUE
          edge_table$kg_relations[i] <- paste(rels, collapse = ";")
          edge_table$kg_n_links[i] <- sum(hit)
        }
      }

      if (any(edge_table$kg_linked, na.rm = TRUE)) {
        rel_all <- unlist(strsplit(edge_table$kg_relations[edge_table$kg_linked], ";", fixed = TRUE))
        kg_relation_summary <- as.data.frame(table(rel_all), stringsAsFactors = FALSE)
        names(kg_relation_summary) <- c("relation", "count")
        kg_relation_summary <- kg_relation_summary[order(kg_relation_summary$count, decreasing = TRUE), , drop = FALSE]
        rownames(kg_relation_summary) <- NULL
      }
    }
  }

  summary <- list(
    n_channels = nrow(net),
    n_candidate_edges = nrow(ut),
    threshold = as.numeric(threshold),
    n_selected_edges = nrow(edge_table)
  )

  list(
    edge_table = edge_table,
    threshold = as.numeric(threshold),
    network_matrix = net,
    summary = summary,
    kg_relation_summary = kg_relation_summary
  )
}

# ---- Internal helpers --------------------------------------------------------

#' @noRd
.emg_extract_matrix <- function(x, channels = NULL, assay_name = NULL) {
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  dims <- dim(data)
  if (length(dims) == 2) {
    signal_data <- data
  } else if (length(dims) >= 3) {
    signal_data <- data[, , 1]
  } else {
    stop("Unsupported assay dimension.", call. = FALSE)
  }

  if (is.null(channels)) channels <- seq_len(ncol(signal_data))
  channels <- as.integer(channels)
  if (any(is.na(channels)) || any(channels < 1L) || any(channels > ncol(signal_data))) {
    stop("Invalid channels.", call. = FALSE)
  }
  signal_data <- signal_data[, channels, drop = FALSE]

  ch_names <- colnames(signal_data)
  if (is.null(ch_names)) {
    ch_names <- paste0("EMG", channels)
  }

  list(signal_data = signal_data, channel_names = ch_names)
}

#' @noRd
.emg_welch_psd <- function(x, nperseg, noverlap, sr) {
  n <- length(x)
  if (n < nperseg) {
    nperseg <- n
    noverlap <- 0L
  }
  step <- nperseg - noverlap
  if (step <= 0L) stop("Invalid step length in Welch PSD.", call. = FALSE)
  n_segments <- floor((n - noverlap) / step)
  n_segments <- max(1L, n_segments)

  n_freqs <- floor(nperseg / 2) + 1
  psd <- numeric(n_freqs)
  window <- 0.5 * (1 - cos(2 * pi * seq(0, nperseg - 1) / (nperseg - 1)))
  used <- 0L

  for (seg in seq_len(n_segments)) {
    start <- (seg - 1L) * step + 1L
    end <- start + nperseg - 1L
    if (end > n) break
    xx <- x[start:end] * window
    fft_x <- stats::fft(xx)
    psd <- psd + (Mod(fft_x[seq_len(n_freqs)])^2)
    used <- used + 1L
  }
  if (used == 0L) stop("Welch PSD could not use any segment.", call. = FALSE)
  psd <- psd / used
  psd <- psd / (sr * sum(window^2))

  list(psd = psd, frequencies = seq(0, sr / 2, length.out = n_freqs))
}

#' @noRd
.emg_welch_csd <- function(x, y, nperseg, noverlap, sr) {
  n <- min(length(x), length(y))
  x <- x[seq_len(n)]
  y <- y[seq_len(n)]
  if (n < nperseg) {
    nperseg <- n
    noverlap <- 0L
  }
  step <- nperseg - noverlap
  if (step <= 0L) stop("Invalid step length in Welch CSD.", call. = FALSE)
  n_segments <- floor((n - noverlap) / step)
  n_segments <- max(1L, n_segments)

  n_freqs <- floor(nperseg / 2) + 1
  csd <- complex(real = rep(0, n_freqs), imaginary = rep(0, n_freqs))
  window <- 0.5 * (1 - cos(2 * pi * seq(0, nperseg - 1) / (nperseg - 1)))
  used <- 0L

  for (seg in seq_len(n_segments)) {
    start <- (seg - 1L) * step + 1L
    end <- start + nperseg - 1L
    if (end > n) break
    xx <- x[start:end] * window
    yy <- y[start:end] * window
    fft_x <- stats::fft(xx)
    fft_y <- stats::fft(yy)
    csd <- csd + fft_x[seq_len(n_freqs)] * Conj(fft_y[seq_len(n_freqs)])
    used <- used + 1L
  }
  if (used == 0L) stop("Welch CSD could not use any segment.", call. = FALSE)
  csd <- csd / used
  csd <- csd / (sr * sum(window^2))

  list(csd = csd, frequencies = seq(0, sr / 2, length.out = n_freqs))
}

#' @noRd
.emg_morlet_wavelet_transform <- function(x, frequencies, n_cycles, sr) {
  n <- length(x)
  n_fft <- stats::nextn(n + max(n, 1024L), factors = 2L)
  fft_x <- stats::fft(c(x, rep(0, n_fft - n)))
  fft_freqs <- seq(0, sr * (1 - 1 / n_fft), length.out = n_fft)

  out <- matrix(complex(0), nrow = n, ncol = length(frequencies))

  for (fi in seq_along(frequencies)) {
    f0 <- frequencies[fi]
    sigma_f <- f0 / n_cycles
    gauss <- exp(-0.5 * ((fft_freqs - f0) / sigma_f)^2)
    gauss[fft_freqs > sr / 2] <- 0
    conv <- stats::fft(fft_x * gauss, inverse = TRUE) / n_fft
    out[, fi] <- conv[seq_len(n)]
  }
  out
}

#' @noRd
.emg_convolve_mirror <- function(x, kernel) {
  n <- length(x)
  k <- length(kernel)
  hw <- floor(k / 2L)
  if (n <= 1L || k <= 1L) return(x)

  idx <- seq_len(n)
  left <- idx[seq.int(hw, 1L, by = -1L)]
  right <- idx[seq.int(n, n - hw + 1L, by = -1L)]
  xp <- c(x[left], x, x[right])
  if (is.complex(xp)) {
    out_re <- stats::filter(Re(xp), kernel, sides = 2)
    out_im <- stats::filter(Im(xp), kernel, sides = 2)
    out <- out_re[(hw + 1L):(hw + n)] + (1i * out_im[(hw + 1L):(hw + n)])
  } else {
    out <- stats::filter(xp, kernel, sides = 2)
    out <- out[(hw + 1L):(hw + n)]
  }
  as.vector(out)
}

#' @noRd
.emg_compute_coi <- function(n, sr, n_cycles) {
  t <- seq(0, (n - 1) / sr, length.out = n)
  edge_dist <- pmin(t, t[n] - t)
  sqrt(2) * n_cycles / (2 * pi * pmax(edge_dist, 1 / sr))
}

#' @noRd
.emg_welch_fft_segments <- function(signal_data, nperseg, noverlap, sr) {
  n <- nrow(signal_data)
  n_channels <- ncol(signal_data)
  if (n < nperseg) {
    nperseg <- n
    noverlap <- 0L
  }
  step <- nperseg - noverlap
  if (step <= 0L) stop("Invalid step length for Welch FFT segments.", call. = FALSE)

  starts <- seq.int(1L, n - nperseg + 1L, by = step)
  if (length(starts) < 1L) {
    starts <- 1L
    nperseg <- n
    noverlap <- 0L
  }

  n_freqs <- floor(nperseg / 2L) + 1L
  window <- 0.5 * (1 - cos(2 * pi * seq(0, nperseg - 1L) / (nperseg - 1L)))
  fft_arr <- array(complex(real = 0, imaginary = 0),
                   dim = c(n_freqs, n_channels, length(starts)))

  used <- 0L
  for (s in starts) {
    idx <- s:(s + nperseg - 1L)
    if (max(idx) > n) next
    seg <- signal_data[idx, , drop = FALSE]
    seg <- sweep(seg, 1L, window, FUN = "*")
    fft_seg <- stats::mvfft(seg)
    used <- used + 1L
    fft_arr[, , used] <- fft_seg[seq_len(n_freqs), , drop = FALSE]
  }
  if (used == 0L) {
    stop("No valid Welch segments were found.", call. = FALSE)
  }
  if (used < dim(fft_arr)[3]) {
    fft_arr <- fft_arr[, , seq_len(used), drop = FALSE]
  }

  list(
    fft = fft_arr,
    frequencies = seq(0, sr / 2, length.out = n_freqs),
    nperseg = nperseg,
    noverlap = noverlap
  )
}

#' @noRd
.emg_safe_solve <- function(mat) {
  out <- try(solve(mat), silent = TRUE)
  if (!inherits(out, "try-error")) {
    return(out)
  }
  jitter <- 1e-8
  for (k in 1:6) {
    out <- try(solve(mat + diag(jitter, nrow(mat))), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
    jitter <- jitter * 10
    if (jitter > 1e-2) break
  }
  solve(mat + diag(1e-2, nrow(mat)))
}

#' @noRd
.emg_pairwise_granger <- function(x, y, lag_order) {
  n <- min(length(x), length(y))
  x <- as.numeric(x[seq_len(n)])
  y <- as.numeric(y[seq_len(n)])
  p <- as.integer(lag_order)
  if (n <= (2L * p + 5L)) {
    stop("Insufficient samples for GC estimation.", call. = FALSE)
  }

  y_emb <- embed(y, p + 1L)
  x_emb <- embed(x, p + 1L)
  target <- y_emb[, 1]
  y_lags <- y_emb[, 2:(p + 1L), drop = FALSE]
  x_lags <- x_emb[, 2:(p + 1L), drop = FALSE]

  xr <- cbind(1, y_lags)
  xf <- cbind(1, y_lags, x_lags)

  fit_r <- .emg_lm_fit_rss(xr, target)
  fit_f <- .emg_lm_fit_rss(xf, target)

  rss_r <- fit_r$rss
  rss_f <- fit_f$rss
  df1 <- p
  df2 <- length(target) - ncol(xf)

  if (!is.finite(rss_r) || !is.finite(rss_f) || rss_f <= .Machine$double.eps || df2 <= 0L) {
    return(list(f_stat = NA_real_, p_value = NA_real_, delta_r2 = NA_real_))
  }

  f_stat <- ((rss_r - rss_f) / df1) / (rss_f / df2)
  if (!is.finite(f_stat) || f_stat < 0) {
    f_stat <- 0
  }
  p_val <- 1 - stats::pf(f_stat, df1 = df1, df2 = df2)

  tss <- sum((target - mean(target))^2)
  if (!is.finite(tss) || tss <= .Machine$double.eps) {
    delta_r2 <- NA_real_
  } else {
    r2_r <- 1 - (rss_r / tss)
    r2_f <- 1 - (rss_f / tss)
    delta_r2 <- max(r2_f - r2_r, 0)
  }

  list(
    f_stat = as.numeric(f_stat),
    p_value = as.numeric(p_val),
    delta_r2 = as.numeric(delta_r2)
  )
}

#' @noRd
.emg_lm_fit_rss <- function(x, y) {
  fit <- stats::lm.fit(x = x, y = y)
  rss <- sum(fit$residuals^2)
  list(rss = rss)
}

#' @noRd
.emg_weighted_clustering <- function(net) {
  n <- nrow(net)
  out <- numeric(n)
  for (i in seq_len(n)) {
    nei <- which(net[i, ] > 0)
    k <- length(nei)
    if (k < 2L) {
      out[i] <- 0
      next
    }
    tri_sum <- 0
    for (a in seq_len(k - 1L)) {
      j <- nei[[a]]
      for (b in (a + 1L):k) {
        h <- nei[[b]]
        wjh <- net[j, h]
        if (wjh <= 0) next
        tri_sum <- tri_sum + (net[i, j] * net[i, h] * wjh)^(1 / 3)
      }
    }
    out[i] <- (2 * tri_sum) / (k * (k - 1))
  }
  out
}

#' @noRd
.emg_shortest_path_lengths <- function(net) {
  n <- nrow(net)
  dist <- matrix(Inf, nrow = n, ncol = n)
  diag(dist) <- 0

  idx <- which(net > 0, arr.ind = TRUE)
  if (nrow(idx) > 0) {
    dist[idx] <- 1 / pmax(net[idx], .Machine$double.eps)
  }

  for (k in seq_len(n)) {
    via_k <- outer(dist[, k], dist[k, ], "+")
    dist <- pmin(dist, via_k)
  }
  dist
}

#' @noRd
.emg_weighted_modularity <- function(net, modules) {
  m <- sum(net) / 2
  if (!is.finite(m) || m <= .Machine$double.eps) {
    return(NA_real_)
  }
  k <- rowSums(net)
  b <- net - (outer(k, k) / (2 * m))
  same <- outer(modules, modules, "==")
  as.numeric(sum(b[same]) / (2 * m))
}

#' @noRd
.emg_detect_modules <- function(net, k_candidates) {
  n <- nrow(net)
  if (n < 3L) {
    out <- rep(1L, n)
    names(out) <- rownames(net)
    return(out)
  }
  k_candidates <- sort(unique(as.integer(k_candidates)))
  k_candidates <- k_candidates[is.finite(k_candidates) & k_candidates >= 2L & k_candidates <= n]
  if (length(k_candidates) == 0L) {
    out <- rep(1L, n)
    names(out) <- rownames(net)
    return(out)
  }

  sim <- net
  mx <- max(sim, na.rm = TRUE)
  if (is.finite(mx) && mx > .Machine$double.eps) {
    sim <- sim / mx
  }
  sim <- pmin(pmax(sim, 0), 1)
  diag(sim) <- 1
  dmat <- 1 - sim
  diag(dmat) <- 0
  hc <- stats::hclust(stats::as.dist(dmat), method = "average")

  best_q <- -Inf
  best <- rep(1L, n)
  for (k in k_candidates) {
    g <- stats::cutree(hc, k = k)
    q <- .emg_weighted_modularity(net, g)
    if (is.finite(q) && q > best_q) {
      best_q <- q
      best <- as.integer(g)
    }
  }
  names(best) <- rownames(net)
  best
}
