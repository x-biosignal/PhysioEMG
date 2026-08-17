# Space-by-time and hierarchical (Tucker) muscle-synergy models.
#
# The package already has spatial NMF synergies (muscleSynergy), a non-negative
# PARAFAC tensor model (muscleSynergyTensor) and time-varying/convolutive
# synergies. Two families complete the picture:
#   * space-by-time decomposition (Delis et al. 2014, sNM3F) -- factorises each
#     trial's muscle x time EMG into SHARED spatial synergies and SHARED temporal
#     synergies combined by a trial-specific activation matrix, separating "which
#     muscles" from "when".
#   * non-negative Tucker (hierarchical) decomposition -- generalises the PARAFAC
#     tensor model with a CORE tensor that lets any spatial synergy combine with
#     any temporal and any trial mode, capturing interactions PARAFAC's diagonal
#     core cannot.
# Dependency-free base R, multiplicative non-negative updates.

# Coerce input to a list of muscle x time matrices + a [muscle,time,trial] array.
.emg_trials <- function(trials) {
  if (is.array(trials) && length(dim(trials)) == 3L) {
    d <- dim(trials)
    lst <- lapply(seq_len(d[3]), function(k) matrix(trials[, , k], d[1], d[2]))
    arr <- trials
  } else if (is.list(trials)) {
    lst <- lapply(trials, as.matrix)
    d <- c(nrow(lst[[1]]), ncol(lst[[1]]), length(lst))
    arr <- array(0, d); for (k in seq_along(lst)) arr[, , k] <- lst[[k]]
  } else stop("`trials` must be a list of muscle x time matrices or a [muscle, time, trial] array.", call. = FALSE)
  if (any(unlist(lst) < 0)) stop("EMG must be non-negative (rectified/enveloped).", call. = FALSE)
  list(list = lst, array = arr, M = d[1], T = d[2], N = d[3])
}

#' Space-by-time muscle synergy decomposition (sNM3F)
#'
#' Decomposes a set of trials' muscle x time EMG into shared SPATIAL synergies
#' (which muscles group together), shared TEMPORAL synergies (when they are
#' active) and a per-trial activation matrix linking them (Delis et al. 2014).
#'
#' @param trials A list of `muscle x time` non-negative matrices (one per trial),
#'   or a `[muscle, time, trial]` array.
#' @param n_spatial,n_temporal Numbers of spatial and temporal synergies.
#' @param max_iter,tol Multiplicative-update iterations and convergence tolerance.
#' @param n_restart Random restarts (best kept; default 5).
#' @param seed Optional RNG seed.
#' @return a `spacetime_synergy`: `spatial` (`muscle x n_spatial`), `temporal`
#'   (`time x n_temporal`), `activation` (list of `n_spatial x n_temporal` per
#'   trial), `vaf`, `iterations`.
#' @references Delis I, et al. (2014) Front Comput Neurosci 8:118.
#' @seealso [muscleSynergy()], [tuckerSynergy()]
#' @export
#' @examples
#' set.seed(1)
#' Ws <- matrix(stats::runif(6 * 2), 6, 2); Wt <- matrix(stats::runif(50 * 2), 50, 2)
#' trials <- lapply(1:10, function(n) Ws %*% (matrix(stats::runif(4), 2, 2)) %*% t(Wt))
#' fit <- spaceByTimeSynergy(trials, 2, 2, n_restart = 2)
#' fit$vaf
spaceByTimeSynergy <- function(trials, n_spatial, n_temporal, max_iter = 500L,
                               tol = 1e-6, n_restart = 5L, seed = NULL) {
  dat <- .emg_trials(trials); Es <- dat$list
  M <- dat$M; Tn <- dat$T; N <- dat$N
  Ns <- n_spatial; Nt <- n_temporal
  normX <- sum(vapply(Es, function(e) sum(e^2), numeric(1)))
  eps <- 1e-10
  if (!is.null(seed)) set.seed(seed)
  best <- NULL
  for (rs in seq_len(n_restart)) {
    Ws <- matrix(stats::runif(M * Ns), M, Ns); Wt <- matrix(stats::runif(Tn * Nt), Tn, Nt)
    As <- lapply(seq_len(N), function(n) matrix(stats::runif(Ns * Nt), Ns, Nt))
    prev <- Inf
    for (it in seq_len(max_iter)) {
      WstWs <- crossprod(Ws); WttWt <- crossprod(Wt)
      # activations
      for (n in seq_len(N)) {
        num <- crossprod(Ws, Es[[n]]) %*% Wt
        den <- WstWs %*% As[[n]] %*% WttWt
        As[[n]] <- As[[n]] * num / (den + eps)
      }
      # spatial
      numW <- matrix(0, M, Ns); denW <- matrix(0, M, Ns)
      for (n in seq_len(N)) {
        Bn <- As[[n]] %*% t(Wt)                       # Ns x T
        numW <- numW + Es[[n]] %*% t(Bn)
        denW <- denW + Ws %*% tcrossprod(Bn)
      }
      Ws <- Ws * numW / (denW + eps)
      # temporal
      numT <- matrix(0, Tn, Nt); denT <- matrix(0, Tn, Nt)
      for (n in seq_len(N)) {
        Cn <- crossprod(As[[n]], t(Ws))               # Nt x M
        numT <- numT + crossprod(Es[[n]], t(Cn))
        denT <- denT + Wt %*% tcrossprod(Cn)
      }
      Wt <- Wt * numT / (denT + eps)
      if (it %% 25L == 0L || it == max_iter) {
        resid <- sum(vapply(seq_len(N), function(n)
          sum((Es[[n]] - Ws %*% As[[n]] %*% t(Wt))^2), numeric(1)))
        if (is.finite(prev) && abs(prev - resid) < tol * (prev + eps)) break
        prev <- resid
      }
    }
    resid <- sum(vapply(seq_len(N), function(n)
      sum((Es[[n]] - Ws %*% As[[n]] %*% t(Wt))^2), numeric(1)))
    if (is.null(best) || resid < best$resid)
      best <- list(Ws = Ws, Wt = Wt, As = As, resid = resid, iter = it)
  }
  structure(list(spatial = best$Ws, temporal = best$Wt, activation = best$As,
                 vaf = 1 - best$resid / normX, iterations = best$iter,
                 n_spatial = Ns, n_temporal = Nt, n_trials = N),
            class = "spacetime_synergy")
}

#' @export
print.spacetime_synergy <- function(x, ...) {
  cat(sprintf("Space-by-time synergies -- %d spatial x %d temporal, %d trials\n",
              x$n_spatial, x$n_temporal, x$n_trials))
  cat(sprintf("  VAF = %.3f (%d iterations)\n", x$vaf, x$iterations))
  invisible(x)
}

# --- tensor helpers (self-consistent unfold/fold/ttm) -----------------------
.tns_unfold <- function(X, n) {
  d <- dim(X); perm <- c(n, setdiff(seq_along(d), n))
  matrix(aperm(X, perm), d[n])
}
.tns_fold <- function(mat, n, dims) {
  perm <- c(n, setdiff(seq_along(dims), n))
  aperm(array(mat, dims[perm]), order(perm))
}
.tns_ttm <- function(X, U, n) {                       # mode-n product X x_n U
  d <- dim(X); Y <- U %*% .tns_unfold(X, n)
  dd <- d; dd[n] <- nrow(U)
  .tns_fold(Y, n, dd)
}

#' Non-negative Tucker (hierarchical) muscle-synergy decomposition
#'
#' A hierarchical tensor synergy model: factorises the `muscle x time x trial`
#' EMG tensor into non-negative spatial, temporal and trial factor matrices plus
#' a CORE tensor that weights every combination of them. Unlike PARAFAC
#' ([muscleSynergyTensor()]), the core is not diagonal, so a spatial synergy can
#' pair with several temporal synergies -- capturing interactions PARAFAC cannot.
#'
#' @param trials A `[muscle, time, trial]` array or a list of `muscle x time`
#'   matrices.
#' @param ranks Length-3 factor ranks `c(spatial, temporal, trial)`.
#' @param max_iter,tol Iterations and convergence tolerance.
#' @param n_restart Random restarts (default 3).
#' @param seed Optional RNG seed.
#' @return a `tucker_synergy`: `spatial` (`muscle x P`), `temporal` (`time x Q`),
#'   `trial` (`trial x R`), `core` (`P x Q x R`), `vaf`, `iterations`.
#' @references Kim YD, Choi S (2007) non-negative Tucker; Tucker LR (1966).
#' @seealso [muscleSynergyTensor()], [spaceByTimeSynergy()]
#' @export
#' @examples
#' set.seed(1)
#' A <- matrix(stats::runif(6*2),6,2); B <- matrix(stats::runif(40*2),40,2); C <- matrix(stats::runif(8*2),8,2)
#' G <- array(stats::runif(2*2*2), c(2,2,2))
#' X <- reconstructTucker(G, A, B, C)
#' tuckerSynergy(X, ranks = c(2, 2, 2), n_restart = 1)$vaf
tuckerSynergy <- function(trials, ranks, max_iter = 300L, tol = 1e-6,
                          n_restart = 3L, seed = NULL) {
  dat <- .emg_trials(trials); X <- dat$array
  I <- dat$M; J <- dat$T; K <- dat$N
  P <- ranks[1]; Q <- ranks[2]; R <- ranks[3]
  normX <- sum(X^2); eps <- 1e-10
  X1 <- .tns_unfold(X, 1); X2 <- .tns_unfold(X, 2); X3 <- .tns_unfold(X, 3)
  if (!is.null(seed)) set.seed(seed)
  best <- NULL
  for (rs in seq_len(n_restart)) {
    A <- matrix(stats::runif(I * P), I, P); B <- matrix(stats::runif(J * Q), J, Q)
    C <- matrix(stats::runif(K * R), K, R); G <- array(stats::runif(P * Q * R), c(P, Q, R))
    prev <- Inf
    for (it in seq_len(max_iter)) {
      ZA <- .tns_unfold(.tns_ttm(.tns_ttm(G, B, 2), C, 3), 1)   # P x (J*K)
      A <- A * (X1 %*% t(ZA)) / (A %*% tcrossprod(ZA) + eps)
      ZB <- .tns_unfold(.tns_ttm(.tns_ttm(G, A, 1), C, 3), 2)   # Q x (I*K)
      B <- B * (X2 %*% t(ZB)) / (B %*% tcrossprod(ZB) + eps)
      ZC <- .tns_unfold(.tns_ttm(.tns_ttm(G, A, 1), B, 2), 3)   # R x (I*J)
      C <- C * (X3 %*% t(ZC)) / (C %*% tcrossprod(ZC) + eps)
      # normalise factor columns into the core (removes the scale indeterminacy
      # that otherwise lets the factors diverge to Inf/NaN over many iterations)
      dm <- function(v) { m <- matrix(0, length(v), length(v)); diag(m) <- v; m }
      nrm <- function(F) { s <- sqrt(colSums(F^2)); s[s < eps] <- 1; s }
      sa <- nrm(A); A <- sweep(A, 2L, sa, "/"); G <- .tns_ttm(G, dm(sa), 1)
      sb <- nrm(B); B <- sweep(B, 2L, sb, "/"); G <- .tns_ttm(G, dm(sb), 2)
      sc <- nrm(C); C <- sweep(C, 2L, sc, "/"); G <- .tns_ttm(G, dm(sc), 3)
      num <- .tns_ttm(.tns_ttm(.tns_ttm(X, t(A), 1), t(B), 2), t(C), 3)
      den <- .tns_ttm(.tns_ttm(.tns_ttm(G, crossprod(A), 1), crossprod(B), 2),
                      crossprod(C), 3)
      G <- G * num / (den + eps)
      if (it %% 20L == 0L || it == max_iter) {
        Xhat <- .tns_ttm(.tns_ttm(.tns_ttm(G, A, 1), B, 2), C, 3)
        resid <- sum((X - Xhat)^2)
        if (is.finite(prev) && abs(prev - resid) < tol * (prev + eps)) break
        prev <- resid
      }
    }
    Xhat <- .tns_ttm(.tns_ttm(.tns_ttm(G, A, 1), B, 2), C, 3)
    resid <- sum((X - Xhat)^2)
    if (is.null(best) || resid < best$resid)
      best <- list(A = A, B = B, C = C, G = G, resid = resid, iter = it)
  }
  structure(list(spatial = best$A, temporal = best$B, trial = best$C,
                 core = best$G, vaf = 1 - best$resid / normX, iterations = best$iter,
                 ranks = ranks), class = "tucker_synergy")
}

#' @export
print.tucker_synergy <- function(x, ...) {
  cat(sprintf("Non-negative Tucker synergies -- ranks (%d, %d, %d)\n",
              x$ranks[1], x$ranks[2], x$ranks[3]))
  cat(sprintf("  VAF = %.3f (%d iterations)\n", x$vaf, x$iterations))
  invisible(x)
}

#' Reconstruct a Tucker model (helper)
#'
#' Rebuilds the `muscle x time x trial` tensor from a Tucker core and factors;
#' exposed mainly for examples and testing.
#'
#' @param core A `P x Q x R` core array.
#' @param spatial,temporal,trial Factor matrices.
#' @return The reconstructed `[muscle, time, trial]` array.
#' @export
reconstructTucker <- function(core, spatial, temporal, trial) {
  .tns_ttm(.tns_ttm(.tns_ttm(core, spatial, 1), temporal, 2), trial, 3)
}
