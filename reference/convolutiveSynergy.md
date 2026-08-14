# Time-varying (convolutive) muscle synergies

Extracts shift-invariant spatiotemporal muscle synergies by convolutive
NMF. Each of the `n_synergies` components is a `muscle x L` template
that the model can place at any time with any amplitude, so a synergy
captures a fixed *spatiotemporal* pattern of muscle activation (d'Avella
time-varying synergies), not just a spatial weighting.

## Usage

``` r
convolutiveSynergy(
  emg,
  n_synergies,
  L,
  max_iter = 300,
  tol = 1e-06,
  restarts = 3,
  seed = NULL
)
```

## Arguments

- emg:

  A non-negative `muscle x time` matrix (rows = muscles).

- n_synergies:

  Number of synergies.

- L:

  Temporal length of each synergy (samples).

- max_iter:

  Maximum multiplicative-update iterations (default 300).

- tol:

  Relative-error convergence tolerance (default 1e-6).

- restarts:

  Random restarts; the best is kept (default 3).

- seed:

  Optional integer seed.

## Value

An object of class `"convolutive_synergy"`: a list with `synergies` (a
`muscle x synergy x L` array), `activations` (`synergy x time`), `vaf`,
`iterations`, `n_synergies`, and `L`.

## References

d'Avella A, Saltiel P, Bizzi E (2003). Combinations of muscle synergies
in the construction of a natural motor behavior. *Nat Neurosci*
6:300-308. Smaragdis P (2004). Non-negative matrix factor deconvolution.
*ICA*.

## Examples

``` r
set.seed(1)
M <- 5; Tt <- 200; N <- 2; L <- 15
W0 <- array(0, c(M, N, L))
for (nn in 1:N) W0[, nn, ] <- outer(runif(M), dnorm(1:L, L / 2, 3))
H0 <- matrix(0, N, Tt); H0[1, c(30, 120)] <- 1; H0[2, c(70, 160)] <- 1
V <- matrix(0, M, Tt)
for (tau in 0:(L - 1)) V <- V +
  W0[, , tau + 1] %*% cbind(matrix(0, N, tau), H0[, 1:(Tt - tau), drop = FALSE])
fit <- convolutiveSynergy(V, n_synergies = 2, L = 15, seed = 1)
fit$vaf
#> [1] 0.9996941
```
