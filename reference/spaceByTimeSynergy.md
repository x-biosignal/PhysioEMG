# Space-by-time muscle synergy decomposition (sNM3F)

Decomposes a set of trials' muscle x time EMG into shared SPATIAL
synergies (which muscles group together), shared TEMPORAL synergies
(when they are active) and a per-trial activation matrix linking them
(Delis et al. 2014).

## Usage

``` r
spaceByTimeSynergy(
  trials,
  n_spatial,
  n_temporal,
  max_iter = 500L,
  tol = 1e-06,
  n_restart = 5L,
  seed = NULL
)
```

## Arguments

- trials:

  A list of `muscle x time` non-negative matrices (one per trial), or a
  `[muscle, time, trial]` array.

- n_spatial, n_temporal:

  Numbers of spatial and temporal synergies.

- max_iter, tol:

  Multiplicative-update iterations and convergence tolerance.

- n_restart:

  Random restarts (best kept; default 5).

- seed:

  Optional RNG seed.

## Value

a `spacetime_synergy`: `spatial` (`muscle x n_spatial`), `temporal`
(`time x n_temporal`), `activation` (list of `n_spatial x n_temporal`
per trial), `vaf`, `iterations`.

## References

Delis I, et al. (2014) Front Comput Neurosci 8:118.

## See also

[`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md),
[`tuckerSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/tuckerSynergy.md)

## Examples

``` r
set.seed(1)
Ws <- matrix(stats::runif(6 * 2), 6, 2); Wt <- matrix(stats::runif(50 * 2), 50, 2)
trials <- lapply(1:10, function(n) Ws %*% (matrix(stats::runif(4), 2, 2)) %*% t(Wt))
fit <- spaceByTimeSynergy(trials, 2, 2, n_restart = 2)
fit$vaf
#> [1] 0.9999996
```
