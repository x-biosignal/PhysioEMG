#' Muscle Synergy Decomposition
#'
#' Decomposes multi-channel EMG into muscle synergies using matrix factorization.
#'
#' @param x A PhysioExperiment object with multi-channel EMG.
#' @param n_synergies Number of synergies to extract.
#' @param method Decomposition method: "nmf" (non-negative matrix factorization),
#'   "pca" (principal component analysis), or "ica" (independent component analysis).
#' @param max_iter Maximum iterations for NMF (default: 200).
#' @param tol Convergence tolerance for NMF (default: 1e-4).
#' @param assay_name Input assay name (default: first assay).
#' @return A list with:
#'   \itemize{
#'     \item \code{W}: Synergy weight matrix (n_synergies x channels)
#'     \item \code{H}: Activation pattern matrix (time x n_synergies)
#'     \item \code{vaf}: Variance accounted for (0-1)
#'     \item \code{method}: Method used
#'     \item \code{original_data}: Original data matrix for reconstruction
#'   }
#' @seealso [synergyReconstruct()] for reconstructing data from synergies,
#'   [synergyCompare()] for comparing synergy solutions,
#'   [emgEnvelope()] for amplitude envelope extraction
#' @references De Luca, C.J. (1997). "The use of surface electromyography in
#'   biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
#'   doi:10.1123/jab.13.2.135
#' @references Merletti, R. & Parker, P.A. (2004). "Electromyography:
#'   Physiology, Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
#'   doi:10.1002/0471678384
#' @export
muscleSynergy <- function(x, n_synergies, method = c("nmf", "pca", "ica"),
                           max_iter = 200L, tol = 1e-4,
                           assay_name = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  method <- match.arg(method)

  if (is.null(assay_name)) assay_name <- defaultAssay(x)
  data <- SummarizedExperiment::assay(x, assay_name)
  n_time <- nrow(data)
  n_channels <- ncol(data)

  if (n_synergies > n_channels) {
    stop("n_synergies cannot exceed number of channels", call. = FALSE)
  }

  if (method == "nmf") {
    data_nn <- pmax(data, 0)

    W <- matrix(abs(rnorm(n_synergies * n_channels, sd = 0.5)),
                nrow = n_synergies, ncol = n_channels)
    H <- matrix(abs(rnorm(n_time * n_synergies, sd = 0.5)),
                nrow = n_time, ncol = n_synergies)

    prev_error <- Inf
    for (iter in seq_len(max_iter)) {
      numerator_H <- data_nn %*% t(W)
      denominator_H <- H %*% (W %*% t(W))
      denominator_H[denominator_H < 1e-10] <- 1e-10
      H <- H * numerator_H / denominator_H
      H[!is.finite(H)] <- 1e-10
      H <- pmax(H, 1e-10)

      numerator_W <- t(H) %*% data_nn
      denominator_W <- (t(H) %*% H) %*% W
      denominator_W[denominator_W < 1e-10] <- 1e-10
      W <- W * numerator_W / denominator_W
      W[!is.finite(W)] <- 1e-10
      W <- pmax(W, 1e-10)

      reconstruction <- H %*% W
      error <- sum((data_nn - reconstruction)^2)
      if (!is.finite(error)) break
      rel_change <- abs(prev_error - error) / (abs(prev_error) + abs(error) + 1e-10)
      if (is.finite(rel_change) && rel_change < tol) break
      prev_error <- error
    }

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
    original_data = original_data
  )
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
