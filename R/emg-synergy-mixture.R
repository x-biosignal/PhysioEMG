# Mixture-based muscle synergy analysis -----------------------------------
# Spatial NMF assumes one global synergy set for everyone. In heterogeneous
# populations, subjects may fall into subgroups with distinct synergy structure.
# muscleSynergyMixture() models each subject's muscle activation with a mixture
# of factor analysers: cluster k has its own loading matrix Lambda_k (that
# subgroup's synergies), and subjects are assigned to the cluster whose synergy
# structure fits them best -- the single-function distillation of the method in
# the synergyMixR package (Matsui, https://github.com/matsui-lab/synergyMixR).
# Algorithm: mixture of factor analysers EM (Ghahramani & Hinton 1996).

#' @keywords internal
#' @noRd
.mvn_logdens <- function(Xc, Sig) {
  # log N(x; 0, Sig) for each row of the already-centred Xc (rows = samples)
  M <- ncol(Xc)
  L <- tryCatch(chol(Sig), error = function(e) chol(Sig + diag(1e-6, M)))
  logdet <- 2 * sum(log(diag(L)))
  z <- backsolve(L, t(Xc), transpose = TRUE)     # M x T
  -0.5 * (M * log(2 * pi) + logdet + colSums(z^2))
}

#' @keywords internal
#' @noRd
.ppca <- function(S, r) {
  # closed-form probabilistic PCA (Tipping & Bishop 1999): isotropic noise
  M <- nrow(S)
  eg <- eigen(S, symmetric = TRUE)
  lam <- pmax(eg$values, 1e-8)
  sigma2 <- if (M > r) mean(lam[(r + 1):M]) else 1e-6
  U <- eg$vectors[, seq_len(r), drop = FALSE]
  W <- U %*% diag(sqrt(pmax(lam[seq_len(r)] - sigma2, 1e-8)), r, r)
  list(Lambda = W, Psi = rep(sigma2, M))
}

#' @keywords internal
#' @noRd
.fa_em <- function(S, Lambda, Psi, iters = 15L) {
  M <- nrow(S); r <- ncol(Lambda)
  for (it in seq_len(iters)) {
    C <- Lambda %*% t(Lambda) + diag(Psi, M)
    Cinv <- tryCatch(solve(C), error = function(e) solve(C + diag(1e-6, M)))
    beta <- t(Lambda) %*% Cinv                    # r x M
    Ezz <- diag(r) - beta %*% Lambda + beta %*% S %*% t(beta)
    Lambda <- (S %*% t(beta)) %*% solve(Ezz)
    Psi <- pmax(diag(S - Lambda %*% beta %*% S), 1e-6)
  }
  list(Lambda = Lambda, Psi = Psi)
}

#' Mixture-based muscle synergy analysis
#'
#' Fits a mixture of factor analysers to per-subject muscle activation data to
#' identify subgroups of subjects with distinct muscle-synergy structure, rather
#' than forcing a single global synergy model on a heterogeneous population.
#' Each cluster's factor loading matrix is that subgroup's synergy set.
#'
#' @param list_of_data A list of non-empty numeric matrices, one per subject,
#'   each `time x muscle` (rows = samples/time, columns = muscles; all subjects
#'   share the same set of muscles/columns).
#' @param n_clusters Number of subgroups (mixture components).
#' @param n_factors Number of synergies (factors) per subgroup.
#' @param method `"mfa"` (mixture of factor analysers, per-muscle diagonal noise;
#'   default) or `"mpca"` (mixture of probabilistic PCA, isotropic noise).
#' @param max_iter Maximum EM iterations (default 100).
#' @param tol Log-likelihood convergence tolerance (default 1e-4).
#' @param n_init Random/k-means restarts; the best log-likelihood is kept
#'   (default 5).
#' @param seed Optional integer seed.
#'
#' @return An object of class `"synergy_mixture"`: a list with `cluster`
#'   (subgroup label per subject), `synergies` (length-`n_clusters` list of
#'   `muscle x factor` loading matrices), `uniqueness` (per-cluster diagonal
#'   noise), `mean` (per-cluster mean), `proportions` (mixing weights),
#'   `responsibilities` (`subject x cluster`), `loglik`, and `bic`.
#'
#' @references
#' Matsui Y. Mixture-Based Muscle Synergy Analysis. SSRN preprint.
#' \url{https://ssrn.com/abstract=6798399}
#'
#' Ghahramani Z, Hinton GE (1996). The EM algorithm for mixtures of factor
#' analyzers. Technical Report CRG-TR-96-1, University of Toronto.
#'
#' Matsui Y. synergyMixR: Mixture-Based Muscle Synergy Analysis.
#' \url{https://github.com/matsui-lab/synergyMixR}
#'
#' @examples
#' set.seed(1)
#' M <- 6; r <- 2
#' LA <- matrix(abs(rnorm(M * r)), M, r)          # subgroup A synergies
#' LB <- matrix(abs(rnorm(M * r)), M, r)          # subgroup B synergies
#' gen <- function(L, n = 60) t(L %*% matrix(rnorm(r * n), r, n) +
#'                              matrix(rnorm(M * n, 0, 0.3), M, n))
#' data <- c(replicate(6, gen(LA), simplify = FALSE),
#'           replicate(6, gen(LB), simplify = FALSE))
#' fit <- muscleSynergyMixture(data, n_clusters = 2, n_factors = 2, seed = 1)
#' fit$cluster
#' @export
muscleSynergyMixture <- function(list_of_data, n_clusters, n_factors,
                                 method = c("mfa", "mpca"),
                                 max_iter = 100, tol = 1e-4, n_init = 5,
                                 seed = NULL) {
  method <- match.arg(method)
  if (!is.list(list_of_data) || length(list_of_data) < n_clusters)
    stop("`list_of_data` must be a list with at least `n_clusters` subjects.",
         call. = FALSE)
  X <- lapply(list_of_data, as.matrix)
  M <- ncol(X[[1]])
  if (any(vapply(X, ncol, integer(1)) != M))
    stop("all subjects must have the same number of muscles (columns).",
         call. = FALSE)
  N <- length(X); K <- n_clusters; r <- n_factors
  Ti <- vapply(X, nrow, integer(1))
  subj_mean <- t(vapply(X, colMeans, numeric(M)))

  init_params <- function(z) {
    mu <- vector("list", K); Lam <- vector("list", K); Psi <- vector("list", K)
    for (k in seq_len(K)) {
      idx <- which(z == k); if (!length(idx)) idx <- sample(N, 1)
      Xk <- do.call(rbind, X[idx])
      mu[[k]] <- colMeans(Xk)
      Xc <- sweep(Xk, 2, mu[[k]])
      pc <- svd(Xc / sqrt(max(nrow(Xc) - 1, 1)), nu = 0, nv = r)
      Lam[[k]] <- pc$v %*% diag(pc$d[seq_len(r)], r, r)
      Psi[[k]] <- pmax(apply(Xc, 2, stats::var) - rowSums(Lam[[k]]^2), 1e-3)
    }
    list(mu = mu, Lambda = Lam, Psi = Psi, pi = rep(1 / K, K))
  }

  best <- NULL
  for (init in seq_len(n_init)) {
    if (!is.null(seed)) set.seed(seed + init)
    z0 <- tryCatch(stats::kmeans(subj_mean, K, nstart = 3)$cluster,
                   error = function(e) sample(seq_len(K), N, replace = TRUE))
    if (init > 1) z0 <- sample(seq_len(K), N, replace = TRUE)
    p <- init_params(z0)
    prev_ll <- -Inf; R <- NULL
    for (iter in seq_len(max_iter)) {
      # E-step: subject-level log responsibilities
      logR <- matrix(0, N, K)
      for (k in seq_len(K)) {
        Sig <- p$Lambda[[k]] %*% t(p$Lambda[[k]]) + diag(p$Psi[[k]], M)
        for (i in seq_len(N))
          logR[i, k] <- log(p$pi[k] + 1e-16) +
            sum(.mvn_logdens(sweep(X[[i]], 2, p$mu[[k]]), Sig))
      }
      mx <- apply(logR, 1, max)
      W <- exp(logR - mx); R <- W / rowSums(W)
      ll <- sum(mx + log(rowSums(W)))

      # M-step
      p$pi <- pmax(colMeans(R), 1e-6); p$pi <- p$pi / sum(p$pi)
      for (k in seq_len(K)) {
        wsum <- sum(R[, k] * Ti)
        mu_k <- Reduce(`+`, lapply(seq_len(N),
                       function(i) R[i, k] * colSums(X[[i]]))) / wsum
        S <- matrix(0, M, M)
        for (i in seq_len(N)) {
          Xc <- sweep(X[[i]], 2, mu_k)
          S <- S + R[i, k] * crossprod(Xc)
        }
        S <- S / wsum
        fa <- if (method == "mpca") .ppca(S, r)
              else .fa_em(S, p$Lambda[[k]], p$Psi[[k]])
        p$mu[[k]] <- mu_k; p$Lambda[[k]] <- fa$Lambda; p$Psi[[k]] <- fa$Psi
      }
      if (abs(ll - prev_ll) < tol) break
      prev_ll <- ll
    }
    if (is.null(best) || ll > best$ll)
      best <- list(p = p, R = R, ll = ll)
  }

  noise_par <- if (method == "mpca") 1 else M
  npar <- K * (M + M * r - r * (r - 1) / 2 + noise_par) + (K - 1)
  bic <- -2 * best$ll + npar * log(N)
  structure(
    list(cluster = apply(best$R, 1, which.max),
         synergies = best$p$Lambda, uniqueness = best$p$Psi,
         mean = best$p$mu, proportions = best$p$pi,
         responsibilities = best$R, loglik = best$ll, bic = bic,
         n_clusters = K, n_factors = r, method = method),
    class = "synergy_mixture")
}

#' @export
print.synergy_mixture <- function(x, ...) {
  cat("<Mixture-based muscle synergy analysis>\n")
  cat(sprintf("  subgroups: %d (%d synergies each)\n",
              x$n_clusters, x$n_factors))
  cat(sprintf("  sizes:     %s\n",
              paste(tabulate(x$cluster, x$n_clusters), collapse = " ")))
  cat(sprintf("  method: %s   logLik: %.1f   BIC: %.1f\n",
              toupper(x$method), x$loglik, x$bic))
  invisible(x)
}

#' Mixture-based muscle synergy analysis on PhysioExperiment objects
#'
#' Thin entry point that runs [muscleSynergyMixture()] directly on ecosystem
#' data: it reads each subject's `time x muscle` activation from a
#' PhysioExperiment (or accepts matrices) and fits the subgroup mixture.
#'
#' @param pe_list A list with one element per subject, each a `PhysioExperiment`
#'   (its assay is read) or a `time x muscle` matrix.
#' @param n_clusters,n_factors Number of subgroups and synergies (factors).
#' @param method `"mfa"` or `"mpca"` (see [muscleSynergyMixture()]).
#' @param assay_name Assay to read from each PhysioExperiment (default: its
#'   default assay).
#' @param ... Further arguments passed to [muscleSynergyMixture()].
#' @return A `"synergy_mixture"` object.
#' @examples
#' set.seed(1)
#' mk <- function(L) {
#'   pe <- PhysioCore::PhysioExperiment(
#'     assays = list(raw = t(L %*% matrix(rnorm(2 * 50), 2, 50) +
#'                           matrix(rnorm(6 * 50, 0, .3), 6, 50))),
#'     samplingRate = 100)
#'   pe
#' }
#' LA <- matrix(abs(rnorm(12)), 6, 2); LB <- matrix(abs(rnorm(12)), 6, 2)
#' pes <- c(replicate(5, mk(LA), simplify = FALSE),
#'          replicate(5, mk(LB), simplify = FALSE))
#' fit <- physioSynergyMixture(pes, n_clusters = 2, n_factors = 2, seed = 1)
#' fit$cluster
#' @export
physioSynergyMixture <- function(pe_list, n_clusters, n_factors,
                                 method = c("mfa", "mpca"),
                                 assay_name = NULL, ...) {
  method <- match.arg(method)
  if (!is.list(pe_list) || !length(pe_list))
    stop("`pe_list` must be a non-empty list.", call. = FALSE)
  list_of_data <- lapply(pe_list, function(pe) {
    if (inherits(pe, "PhysioExperiment")) {
      an <- if (is.null(assay_name)) defaultAssay(pe) else assay_name
      as.matrix(SummarizedExperiment::assay(pe, an))
    } else {
      as.matrix(pe)
    }
  })
  muscleSynergyMixture(list_of_data, n_clusters, n_factors,
                       method = method, ...)
}
