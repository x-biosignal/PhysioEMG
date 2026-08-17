# Time-varying (convolutive) muscle synergies -----------------------------
# Spatial NMF (emg-decompose.R) gives time-invariant synergies; the tensor
# decomposition (emg-tensor.R) adds a trial mode. This file adds shift-invariant
# time-varying synergies (d'Avella & Bizzi 2003), extracted by convolutive NMF
# (Smaragdis 2004; O'Grady & Pearlmutter 2006): each synergy is a short
# muscle x time spatiotemporal template that can be recruited at any latency.

#' @keywords internal
#' @noRd
.shift_right <- function(A, tau) {
  if (tau == 0L) return(A)
  n <- ncol(A)
  cbind(matrix(0, nrow(A), tau), A[, seq_len(n - tau), drop = FALSE])
}

#' @keywords internal
#' @noRd
.shift_left <- function(A, tau) {
  if (tau == 0L) return(A)
  n <- ncol(A)
  cbind(A[, (tau + 1L):n, drop = FALSE], matrix(0, nrow(A), tau))
}

#' Time-varying (convolutive) muscle synergies
#'
#' Extracts shift-invariant spatiotemporal muscle synergies by convolutive NMF.
#' Each of the `n_synergies` components is a `muscle x L` template that the model
#' can place at any time with any amplitude, so a synergy captures a fixed
#' *spatiotemporal* pattern of muscle activation (d'Avella time-varying
#' synergies), not just a spatial weighting.
#'
#' @param emg A non-negative `muscle x time` matrix (rows = muscles).
#' @param n_synergies Number of synergies.
#' @param L Temporal length of each synergy (samples).
#' @param max_iter Maximum multiplicative-update iterations (default 300).
#' @param tol Relative-error convergence tolerance (default 1e-6).
#' @param restarts Random restarts; the best is kept (default 3).
#' @param seed Optional integer seed.
#'
#' @return An object of class `"convolutive_synergy"`: a list with `synergies`
#'   (a `muscle x synergy x L` array), `activations` (`synergy x time`), `vaf`,
#'   `iterations`, `n_synergies`, and `L`.
#'
#' @references
#' d'Avella A, Saltiel P, Bizzi E (2003). Combinations of muscle synergies in the
#' construction of a natural motor behavior. \emph{Nat Neurosci} 6:300-308.
#' Smaragdis P (2004). Non-negative matrix factor deconvolution. \emph{ICA}.
#'
#' @examples
#' set.seed(1)
#' M <- 5; Tt <- 200; N <- 2; L <- 15
#' W0 <- array(0, c(M, N, L))
#' for (nn in 1:N) W0[, nn, ] <- outer(runif(M), dnorm(1:L, L / 2, 3))
#' H0 <- matrix(0, N, Tt); H0[1, c(30, 120)] <- 1; H0[2, c(70, 160)] <- 1
#' V <- matrix(0, M, Tt)
#' for (tau in 0:(L - 1)) V <- V +
#'   W0[, , tau + 1] %*% cbind(matrix(0, N, tau), H0[, 1:(Tt - tau), drop = FALSE])
#' fit <- convolutiveSynergy(V, n_synergies = 2, L = 15, seed = 1)
#' fit$vaf
#' @export
convolutiveSynergy <- function(emg, n_synergies, L, max_iter = 300,
                               tol = 1e-6, restarts = 3, seed = NULL) {
  V <- as.matrix(emg)
  if (any(V < 0, na.rm = TRUE))
    stop("`emg` must be non-negative (rectified/enveloped).", call. = FALSE)
  M <- nrow(V); Tt <- ncol(V); N <- n_synergies
  if (L >= Tt) stop("`L` must be shorter than the recording.", call. = FALSE)
  eps <- .Machine$double.eps
  normV <- sum(V^2)

  recon <- function(W, H) {
    Lam <- matrix(0, M, Tt)
    for (tau in 0:(L - 1)) Lam <- Lam + W[, , tau + 1] %*% .shift_right(H, tau)
    Lam
  }

  best <- NULL
  for (rs in seq_len(restarts)) {
    if (!is.null(seed)) set.seed(seed + rs)
    W <- array(stats::runif(M * N * L), c(M, N, L))
    H <- matrix(stats::runif(N * Tt), N, Tt)
    prev <- Inf; iter <- 0L
    for (it in seq_len(max_iter)) {
      iter <- it
      Lam <- pmax(recon(W, H), eps)
      for (tau in 0:(L - 1)) {                       # update each W_tau
        RtH <- .shift_right(H, tau)
        W[, , tau + 1] <- W[, , tau + 1] *
          (V %*% t(RtH)) / (Lam %*% t(RtH) + eps)
      }
      Lam <- pmax(recon(W, H), eps)
      num <- matrix(0, N, Tt); den <- matrix(0, N, Tt)
      for (tau in 0:(L - 1)) {                       # update H
        Wt <- W[, , tau + 1]
        num <- num + .shift_left(t(Wt) %*% V, tau)
        den <- den + .shift_left(t(Wt) %*% Lam, tau)
      }
      H <- H * num / (den + eps)

      err <- sum((V - recon(W, H))^2) / normV
      if (abs(prev - err) < tol) break
      prev <- err
    }
    resid <- sum((V - recon(W, H))^2)
    if (is.null(best) || resid < best$resid)
      best <- list(W = W, H = H, resid = resid, iter = iter)
  }

  structure(
    list(synergies = best$W, activations = best$H,
         vaf = 1 - best$resid / normV, iterations = best$iter,
         n_synergies = N, L = L),
    class = "convolutive_synergy")
}

#' @export
print.convolutive_synergy <- function(x, ...) {
  cat("<Time-varying (convolutive) muscle synergies>\n")
  cat(sprintf("  synergies: %d (each %d muscles x %d samples)\n",
              x$n_synergies, dim(x$synergies)[1], x$L))
  cat(sprintf("  time:      %d samples\n", ncol(x$activations)))
  cat(sprintf("  VAF:       %.3f  (%d iterations)\n", x$vaf, x$iterations))
  invisible(x)
}
