# Non-negative tensor factorisation of muscle synergies -------------------
# Spatial NMF (emg-decompose.R) factorises a muscle x time matrix. Many designs
# are three-way -- muscle x time x (trial | condition | subject) -- and a tensor
# decomposition keeps that structure: each synergy has a single muscle weighting
# and temporal profile, modulated per trial by a loading. This file implements
# non-negative CANDECOMP/PARAFAC (Cichocki et al. 2009) by multiplicative
# updates, the multiway generalisation of muscle-synergy NMF.

#' @keywords internal
#' @noRd
.khatri_rao <- function(A, B) {
  # column-wise Kronecker: (nrow(A)*nrow(B)) x ncol, column r = kron(A[,r], B[,r])
  R <- ncol(A)
  vapply(seq_len(R), function(r) kronecker(A[, r], B[, r]),
         numeric(nrow(A) * nrow(B)))
}

#' Muscle-synergy tensor decomposition (non-negative PARAFAC)
#'
#' Non-negative CANDECOMP/PARAFAC of a three-way EMG tensor
#' `muscle x time x trial`. Each of the `n_synergies` components is a rank-one
#' term with a muscle weight vector, a temporal activation profile, and a
#' per-trial loading — the multiway analogue of a muscle synergy that captures
#' how the same spatial/temporal modules are reused and rescaled across trials or
#' conditions.
#'
#' @param tensor A non-negative numeric array with dimensions
#'   `muscle x time x trial`.
#' @param n_synergies Number of synergies (rank) to extract.
#' @param max_iter Maximum multiplicative-update iterations (default 500).
#' @param tol Relative reconstruction-error tolerance for convergence
#'   (default 1e-6).
#' @param restarts Random restarts; the best (lowest error) is kept (default 5).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `"synergy_tensor"`: a list with `muscle_weights`
#'   (`muscle x synergy`), `temporal` (`time x synergy`), `trial_loadings`
#'   (`trial x synergy`), `vaf` (variance accounted for), `iterations`, and
#'   `n_synergies`.
#'
#' @references Cichocki A, Zdunek R, Phan AH, Amari S (2009). Nonnegative Matrix
#'   and Tensor Factorizations. Wiley.
#'
#' @examples
#' set.seed(1)
#' M <- 6; Tt <- 50; K <- 8; R <- 2
#' A0 <- matrix(runif(M * R), M, R); B0 <- matrix(runif(Tt * R), Tt, R)
#' C0 <- matrix(runif(K * R), K, R)
#' X <- array(0, c(M, Tt, K))
#' for (r in 1:R) for (k in 1:K) X[, , k] <- X[, , k] +
#'   (A0[, r] %o% B0[, r]) * C0[k, r]
#' fit <- muscleSynergyTensor(X, n_synergies = 2, restarts = 3, seed = 1)
#' fit$vaf
#' @export
muscleSynergyTensor <- function(tensor, n_synergies, max_iter = 500,
                                tol = 1e-6, restarts = 5, seed = NULL) {
  X <- tensor
  if (length(dim(X)) != 3L)
    stop("`tensor` must be a 3-way array (muscle x time x trial).", call. = FALSE)
  if (any(X < 0, na.rm = TRUE))
    stop("`tensor` must be non-negative.", call. = FALSE)
  R <- n_synergies
  d <- dim(X); I <- d[1]; J <- d[2]; K <- d[3]
  eps <- .Machine$double.eps

  # mode unfoldings (column order matches the Khatri-Rao products below)
  X1 <- matrix(X, nrow = I)                       # I x (J*K)
  X2 <- matrix(aperm(X, c(2, 1, 3)), nrow = J)    # J x (I*K)
  X3 <- matrix(aperm(X, c(3, 1, 2)), nrow = K)    # K x (I*J)
  normX <- sum(X^2)

  best <- NULL
  for (rs in seq_len(restarts)) {
    if (!is.null(seed)) set.seed(seed + rs)
    A <- matrix(stats::runif(I * R), I, R)
    B <- matrix(stats::runif(J * R), J, R)
    C <- matrix(stats::runif(K * R), K, R)
    prev_err <- Inf; iter <- 0L
    for (it in seq_len(max_iter)) {
      iter <- it
      krCB <- .khatri_rao(C, B)                    # (K*J) x R -> matches X1 cols
      A <- A * (X1 %*% krCB) / (A %*% (t(krCB) %*% krCB) + eps)
      krCA <- .khatri_rao(C, A)                    # (K*I) x R -> matches X2 cols
      B <- B * (X2 %*% krCA) / (B %*% (t(krCA) %*% krCA) + eps)
      krBA <- .khatri_rao(B, A)                    # (J*I) x R -> matches X3 cols
      C <- C * (X3 %*% krBA) / (C %*% (t(krBA) %*% krBA) + eps)

      resid <- sum((X1 - A %*% t(.khatri_rao(C, B)))^2)
      err <- resid / normX
      if (abs(prev_err - err) < tol) break
      prev_err <- err
    }
    resid <- sum((X1 - A %*% t(.khatri_rao(C, B)))^2)
    if (is.null(best) || resid < best$resid)
      best <- list(A = A, B = B, C = C, resid = resid, iter = iter)
  }

  # normalise: unit-sum muscle weights, push scale into trial loadings
  scaleA <- colSums(best$A); scaleA[scaleA == 0] <- 1
  scaleB <- sqrt(colSums(best$B^2)); scaleB[scaleB == 0] <- 1
  A <- sweep(best$A, 2, scaleA, "/")
  B <- sweep(best$B, 2, scaleB, "/")
  C <- sweep(best$C, 2, 1 / (scaleA * scaleB), "/")

  structure(
    list(muscle_weights = A, temporal = B, trial_loadings = C,
         vaf = 1 - best$resid / normX, iterations = best$iter,
         n_synergies = R),
    class = "synergy_tensor")
}

#' @export
print.synergy_tensor <- function(x, ...) {
  cat("<Muscle-synergy tensor (non-negative PARAFAC)>\n")
  cat(sprintf("  synergies: %d\n", x$n_synergies))
  cat(sprintf("  muscles x time x trials: %d x %d x %d\n",
              nrow(x$muscle_weights), nrow(x$temporal), nrow(x$trial_loadings)))
  cat(sprintf("  VAF: %.3f  (%d iterations)\n", x$vaf, x$iterations))
  invisible(x)
}

#' Select the number of tensor synergies by a VAF threshold
#'
#' Fits non-negative PARAFAC for increasing rank and returns the smallest number
#' of synergies whose reconstruction VAF reaches `vaf_threshold`.
#'
#' @param tensor A non-negative `muscle x time x trial` array.
#' @param max_synergies Largest rank to try (default 6).
#' @param vaf_threshold VAF target (default 0.90).
#' @param ... Passed to [muscleSynergyTensor()].
#' @return A list with `n_synergies` (selected), `vaf` (per rank), and `fit`
#'   (the selected `"synergy_tensor"`).
#' @export
muscleSynergyTensorOrder <- function(tensor, max_synergies = 6,
                                     vaf_threshold = 0.90, ...) {
  vafs <- numeric(max_synergies)
  fits <- vector("list", max_synergies)
  for (r in seq_len(max_synergies)) {
    fits[[r]] <- muscleSynergyTensor(tensor, n_synergies = r, ...)
    vafs[r] <- fits[[r]]$vaf
  }
  sel <- which(vafs >= vaf_threshold)[1]
  if (is.na(sel)) sel <- max_synergies
  list(n_synergies = sel, vaf = vafs, fit = fits[[sel]])
}
