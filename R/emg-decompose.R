#' Muscle Synergy Decomposition
#'
#' Decomposes multi-channel EMG into muscle synergies using matrix factorization.
#'
#' @param x A PhysioExperiment object with multi-channel EMG.
#' @param n_synergies Number of synergies to extract.
#' @param method Decomposition method: "nmf" (non-negative matrix factorization),
#'   "pca" (principal component analysis), or "ica" (independent component analysis).
#' @param n_restarts Number of random NMF restarts; the fit with the lowest
#'   reconstruction error (highest VAF) is returned (default: 1). Ignored for
#'   "pca".
#' @param max_iter Maximum iterations for NMF (default: 200).
#' @param tol Convergence tolerance for NMF (default: 1e-4).
#' @param seed Optional integer seed for reproducible random initialization
#'   (NMF restarts and ICA). If set, the returned best-of-restarts solution is
#'   deterministic.
#' @param assay_name Input assay name (default: first assay).
#' @return A list with:
#'   \itemize{
#'     \item \code{W}: Synergy weight matrix (n_synergies x channels)
#'     \item \code{H}: Activation pattern matrix (time x n_synergies)
#'     \item \code{vaf}: Variance accounted for (0-1)
#'     \item \code{method}: Method used
#'     \item \code{n_restarts}: Number of restarts used
#'     \item \code{all_vaf}: VAF of each NMF restart (NA for pca/ica)
#'     \item \code{convergence}: For NMF, a list with the best restart's
#'       \code{iterations} and \code{converged} flag (NULL otherwise)
#'     \item \code{original_data}: Original data matrix for reconstruction
#'   }
#' @seealso [muscleSynergyOrder()] for model-order (synergy count) selection,
#'   [synergyReconstruct()] for reconstructing data from synergies,
#'   [synergyCompare()] for comparing synergy solutions
#' @references Lee, D.D. & Seung, H.S. (1999). "Learning the parts of objects by
#'   non-negative matrix factorization." Nature, 401(6755), 788-791.
#'   doi:10.1038/44565
#' @references Tresch, M.C., Cheung, V.C.K. & d'Avella, A. (2006). "Matrix
#'   factorization algorithms for the identification of muscle synergies."
#'   Journal of Neurophysiology, 95(4), 2199-2212. doi:10.1152/jn.00222.2005
#' @export
muscleSynergy <- function(x, n_synergies, method = c("nmf", "pca", "ica"),
                           n_restarts = 1L, max_iter = 200L, tol = 1e-4,
                           seed = NULL, assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  if (n_synergies > n_channels) {
    stop("n_synergies cannot exceed number of channels", call. = FALSE)
  }

  convergence <- NULL
  all_vaf <- NA_real_

  if (method == "nmf") {
    data_nn <- pmax(data, 0)
    sst_nn <- sum((data_nn - mean(data_nn))^2)

    best <- NULL
    all_vaf <- numeric(max(1L, as.integer(n_restarts)))
    for (rs in seq_along(all_vaf)) {
      fit <- .nmf_multiplicative(data_nn, n_synergies, max_iter, tol)
      all_vaf[rs] <- if (sst_nn > 0) 1 - fit$error / sst_nn else 1
      if (is.null(best) || fit$error < best$error) best <- fit
    }
    W <- best$W
    H <- best$H
    convergence <- list(iterations = best$iter, converged = best$converged)
    original_data <- data_nn

  } else if (method == "pca") {
    pca_result <- prcomp(data, center = TRUE, scale. = FALSE, rank. = n_synergies)
    W <- t(pca_result$rotation[, seq_len(n_synergies)])
    H <- pca_result$x[, seq_len(n_synergies)]
    original_data <- data

  } else if (method == "ica") {
    centered <- scale(data, center = TRUE, scale = FALSE)
    pca_result <- prcomp(centered, rank. = n_synergies)
    whitened <- pca_result$x[, seq_len(n_synergies)]

    n <- nrow(whitened)
    p <- ncol(whitened)
    Wica <- matrix(rnorm(p * p), p, p)
    svd_w <- svd(Wica)
    Wica <- svd_w$u %*% t(svd_w$v)

    for (iter in seq_len(max_iter)) {
      W_old <- Wica
      for (i in seq_len(p)) {
        wx <- whitened %*% Wica[i, ]
        gwx <- tanh(wx)
        g_prime <- 1 - gwx^2
        Wica[i, ] <- colMeans(whitened * as.vector(gwx)) - mean(g_prime) * Wica[i, ]
      }
      svd_w <- svd(Wica)
      Wica <- svd_w$u %*% t(svd_w$v)
      if (max(abs(abs(rowSums(Wica * W_old)) - 1)) < tol) break
    }

    H <- whitened %*% t(Wica)
    W <- Wica %*% t(pca_result$rotation[, seq_len(n_synergies)])
    original_data <- data
  }

  reconstruction <- H %*% W
  ss_total <- sum((original_data - mean(original_data))^2)
  ss_resid <- sum((original_data - reconstruction)^2)
  vaf <- 1 - ss_resid / ss_total

  list(
    W = W,
    H = H,
    vaf = vaf,
    method = method,
    n_restarts = if (method == "nmf") length(all_vaf) else 1L,
    all_vaf = all_vaf,
    convergence = convergence,
    original_data = original_data
  )
}

# Single Lee-Seung multiplicative NMF run with random initialization.
.nmf_multiplicative <- function(data_nn, k, max_iter, tol) {
  n_time <- nrow(data_nn)
  n_ch <- ncol(data_nn)
  W <- matrix(abs(rnorm(k * n_ch, sd = 0.5)), nrow = k, ncol = n_ch)
  H <- matrix(abs(rnorm(n_time * k, sd = 0.5)), nrow = n_time, ncol = k)

  prev_error <- Inf
  converged <- FALSE
  last_iter <- 0L
  for (iter in seq_len(max_iter)) {
    last_iter <- iter
    num_H <- data_nn %*% t(W)
    den_H <- H %*% (W %*% t(W))
    den_H[den_H < 1e-10] <- 1e-10
    H <- H * num_H / den_H
    H[!is.finite(H)] <- 1e-10
    H <- pmax(H, 1e-10)

    num_W <- t(H) %*% data_nn
    den_W <- (t(H) %*% H) %*% W
    den_W[den_W < 1e-10] <- 1e-10
    W <- W * num_W / den_W
    W[!is.finite(W)] <- 1e-10
    W <- pmax(W, 1e-10)

    error <- sum((data_nn - H %*% W)^2)
    if (!is.finite(error)) break
    rel_change <- abs(prev_error - error) /
      (abs(prev_error) + abs(error) + 1e-10)
    if (is.finite(rel_change) && rel_change < tol) {
      converged <- TRUE
      break
    }
    prev_error <- error
  }
  list(W = W, H = H, error = sum((data_nn - H %*% W)^2),
       iter = last_iter, converged = converged)
}

# Kneedle-style knee: point of maximum perpendicular distance from the chord
# joining the first and last points of the VAF curve.
.vaf_knee <- function(ns, vafs) {
  if (length(ns) <= 2L) return(ns[which.max(vafs)])
  x1 <- ns[1]; y1 <- vafs[1]
  x2 <- ns[length(ns)]; y2 <- vafs[length(ns)]
  den <- sqrt((y2 - y1)^2 + (x2 - x1)^2)
  if (den == 0) return(ns[which.max(vafs)])
  d <- abs((y2 - y1) * ns - (x2 - x1) * vafs + x2 * y1 - y2 * x1) / den
  ns[which.max(d)]
}

# Per-channel variance accounted for by a reconstruction.
.per_channel_vaf <- function(original, recon) {
  vapply(seq_len(ncol(original)), function(j) {
    o <- original[, j]
    sst <- sum((o - mean(o))^2)
    if (sst <= 0) return(1)
    1 - sum((o - recon[, j])^2) / sst
  }, numeric(1))
}

#' Select Muscle Synergy Model Order
#'
#' Sweeps the number of synergies from 1 to \code{max_synergies}, fits a
#' multi-restart NMF at each, and selects the model order from the
#' variance-accounted-for (VAF) curve. The selected order is the smallest number
#' of synergies whose global VAF reaches \code{vaf_threshold}; if the threshold
#' is never reached, the knee (elbow) of the VAF curve is used instead.
#'
#' @param x A PhysioExperiment object with multi-channel EMG.
#' @param max_synergies Maximum number of synergies to test (capped at the
#'   number of channels).
#' @param vaf_threshold Global VAF threshold for selection (default: 0.90).
#' @param per_synergy_vaf Per-channel VAF each muscle should reach at the
#'   selected order; reported as a quality check (default: 0.75).
#' @param n_restarts Number of NMF restarts per synergy count (default: 10).
#' @param max_iter Maximum NMF iterations (default: 200).
#' @param tol NMF convergence tolerance (default: 1e-4).
#' @param seed Optional integer seed making the whole sweep reproducible.
#' @param assay_name Input assay name (default: first assay).
#' @return A list with:
#'   \describe{
#'     \item{vaf_curve}{A data.frame with columns \code{n_synergies} and
#'       \code{vaf}.}
#'     \item{selected}{The selected number of synergies.}
#'     \item{selection_rule}{"vaf_threshold" or "knee".}
#'     \item{vaf_threshold, per_synergy_vaf}{The thresholds used.}
#'     \item{per_channel_vaf}{Per-channel VAF at the selected order.}
#'     \item{per_synergy_vaf_met}{TRUE if every channel's VAF meets
#'       \code{per_synergy_vaf}.}
#'     \item{fit}{The [muscleSynergy()] result at the selected order.}
#'   }
#' @seealso [muscleSynergy()] for the underlying decomposition
#' @references Cheung, V.C.K. et al. (2005). "Central and sensory contributions
#'   to the activation and organization of muscle synergies during natural motor
#'   behaviors." Journal of Neuroscience, 25(27), 6419-6434.
#'   doi:10.1523/JNEUROSCI.4904-04.2005
#' @references Tresch, M.C., Cheung, V.C.K. & d'Avella, A. (2006). "Matrix
#'   factorization algorithms for the identification of muscle synergies."
#'   Journal of Neurophysiology, 95(4), 2199-2212. doi:10.1152/jn.00222.2005
#' @export
#' @examples
#' # three synergies driving disjoint muscle pairs
#' set.seed(1)
#' tt <- seq_len(500)
#' H <- sapply(c(150, 250, 350), function(c0) exp(-((tt - c0) / 40)^2))
#' W <- rbind(c(1, 1, 0, 0, 0, 0), c(0, 0, 1, 1, 0, 0), c(0, 0, 0, 0, 1, 1))
#' data <- H %*% W + matrix(abs(rnorm(500 * 6, sd = 0.02)), 500, 6)
#' pe <- PhysioExperiment(assays = list(env = data), samplingRate = 1000)
#' ord <- muscleSynergyOrder(pe, max_synergies = 5, seed = 1)
#' ord$selected
#' ord$vaf_curve
muscleSynergyOrder <- function(x, max_synergies, vaf_threshold = 0.90,
                               per_synergy_vaf = 0.75, n_restarts = 10L,
                               max_iter = 200L, tol = 1e-4, seed = NULL,
                               assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_channels <- ncol(data)
  max_synergies <- min(as.integer(max_synergies), n_channels)
  stopifnot(max_synergies >= 1L)
  if (!is.null(seed)) set.seed(seed)

  ns <- seq_len(max_synergies)
  vafs <- numeric(length(ns))
  fits <- vector("list", length(ns))
  for (i in ns) {
    fit <- muscleSynergy(x, n_synergies = i, method = "nmf",
                         n_restarts = n_restarts, max_iter = max_iter,
                         tol = tol, seed = NULL, assay_name = assay_name)
    vafs[i] <- fit$vaf
    fits[[i]] <- fit
  }
  vaf_curve <- data.frame(n_synergies = ns, vaf = vafs,
                          stringsAsFactors = FALSE)

  reached <- which(vafs >= vaf_threshold)
  if (length(reached) > 0) {
    selected <- reached[1]
    rule <- "vaf_threshold"
  } else {
    selected <- .vaf_knee(ns, vafs)
    rule <- "knee"
  }

  sel <- fits[[selected]]
  pcv <- .per_channel_vaf(sel$original_data, sel$H %*% sel$W)

  list(vaf_curve = vaf_curve, selected = selected, selection_rule = rule,
       vaf_threshold = vaf_threshold, per_synergy_vaf = per_synergy_vaf,
       per_channel_vaf = pcv, per_synergy_vaf_met = all(pcv >= per_synergy_vaf),
       fit = sel)
}

#' Reconstruct Data from Synergies
#'
#' Reconstructs EMG data using a subset of synergies.
#'
#' @param synergy_result Result from \code{\link{muscleSynergy}}.
#' @param n_synergies Number of synergies to use for reconstruction.
#' @return A list with:
#'   \itemize{
#'     \item \code{reconstructed}: Reconstructed data matrix
#'     \item \code{vaf}: VAF of the reconstruction
#'   }
#' @seealso [muscleSynergy()] for computing the initial decomposition,
#'   [synergyCompare()] for comparing synergy solutions
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @export
synergyReconstruct <- function(synergy_result, n_synergies) {
  stopifnot(is.list(synergy_result))
  stopifnot(all(c("W", "H", "original_data") %in% names(synergy_result)))

  W <- synergy_result$W[seq_len(n_synergies), , drop = FALSE]
  H <- synergy_result$H[, seq_len(n_synergies), drop = FALSE]
  original <- synergy_result$original_data

  reconstructed <- H %*% W

  ss_total <- sum((original - mean(original))^2)
  ss_resid <- sum((original - reconstructed)^2)
  vaf <- 1 - ss_resid / ss_total

  list(reconstructed = reconstructed, vaf = vaf)
}

#' Compare Two Synergy Results
#'
#' Computes pairwise correlation between synergy weight vectors from two
#' decompositions. Uses best-match pairing.
#'
#' @param result1 First result from \code{\link{muscleSynergy}}.
#' @param result2 Second result from \code{\link{muscleSynergy}}.
#' @return A data.frame with columns: synergy1, synergy2, correlation.
#' @seealso [muscleSynergy()] for computing synergy decompositions,
#'   [synergyReconstruct()] for reconstructing data from synergies
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
synergyCompare <- function(result1, result2) {
  W1 <- result1$W
  W2 <- result2$W
  n1 <- nrow(W1)
  n2 <- nrow(W2)

  cor_matrix <- matrix(NA_real_, nrow = n1, ncol = n2)
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      cor_matrix[i, j] <- cor(W1[i, ], W2[j, ])
    }
  }

  results <- list()
  used1 <- logical(n1)
  used2 <- logical(n2)
  n_pairs <- min(n1, n2)

  for (p in seq_len(n_pairs)) {
    best_cor <- -Inf
    best_i <- 0
    best_j <- 0
    for (i in seq_len(n1)) {
      if (used1[i]) next
      for (j in seq_len(n2)) {
        if (used2[j]) next
        if (abs(cor_matrix[i, j]) > best_cor) {
          best_cor <- abs(cor_matrix[i, j])
          best_i <- i
          best_j <- j
        }
      }
    }
    used1[best_i] <- TRUE
    used2[best_j] <- TRUE
    results[[p]] <- data.frame(
      synergy1 = best_i, synergy2 = best_j,
      correlation = cor_matrix[best_i, best_j],
      stringsAsFactors = FALSE)
  }

  do.call(rbind, results)
}
