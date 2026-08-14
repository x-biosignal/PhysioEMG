# Muscle-synergy tensor decomposition (non-negative PARAFAC)

Non-negative CANDECOMP/PARAFAC of a three-way EMG tensor
`muscle x time x trial`. Each of the `n_synergies` components is a
rank-one term with a muscle weight vector, a temporal activation
profile, and a per-trial loading — the multiway analogue of a muscle
synergy that captures how the same spatial/temporal modules are reused
and rescaled across trials or conditions.

## Usage

``` r
muscleSynergyTensor(
  tensor,
  n_synergies,
  max_iter = 500,
  tol = 1e-06,
  restarts = 5,
  seed = NULL
)
```

## Arguments

- tensor:

  A non-negative numeric array with dimensions `muscle x time x trial`.

- n_synergies:

  Number of synergies (rank) to extract.

- max_iter:

  Maximum multiplicative-update iterations (default 500).

- tol:

  Relative reconstruction-error tolerance for convergence (default
  1e-6).

- restarts:

  Random restarts; the best (lowest error) is kept (default 5).

- seed:

  Optional integer seed for reproducibility.

## Value

An object of class `"synergy_tensor"`: a list with `muscle_weights`
(`muscle x synergy`), `temporal` (`time x synergy`), `trial_loadings`
(`trial x synergy`), `vaf` (variance accounted for), `iterations`, and
`n_synergies`.

## References

Cichocki A, Zdunek R, Phan AH, Amari S (2009). Nonnegative Matrix and
Tensor Factorizations. Wiley.

## Examples

``` r
set.seed(1)
M <- 6; Tt <- 50; K <- 8; R <- 2
A0 <- matrix(runif(M * R), M, R); B0 <- matrix(runif(Tt * R), Tt, R)
C0 <- matrix(runif(K * R), K, R)
X <- array(0, c(M, Tt, K))
for (r in 1:R) for (k in 1:K) X[, , k] <- X[, , k] +
  (A0[, r] %o% B0[, r]) * C0[k, r]
fit <- muscleSynergyTensor(X, n_synergies = 2, restarts = 3, seed = 1)
fit$vaf
#> [1] 0.9999953
```
