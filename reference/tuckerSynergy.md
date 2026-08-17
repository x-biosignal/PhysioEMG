# Non-negative Tucker (hierarchical) muscle-synergy decomposition

A hierarchical tensor synergy model: factorises the
`muscle x time x trial` EMG tensor into non-negative spatial, temporal
and trial factor matrices plus a CORE tensor that weights every
combination of them. Unlike PARAFAC
([`muscleSynergyTensor()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyTensor.md)),
the core is not diagonal, so a spatial synergy can pair with several
temporal synergies – capturing interactions PARAFAC cannot.

## Usage

``` r
tuckerSynergy(
  trials,
  ranks,
  max_iter = 300L,
  tol = 1e-06,
  n_restart = 3L,
  seed = NULL
)
```

## Arguments

- trials:

  A `[muscle, time, trial]` array or a list of `muscle x time` matrices.

- ranks:

  Length-3 factor ranks `c(spatial, temporal, trial)`.

- max_iter, tol:

  Iterations and convergence tolerance.

- n_restart:

  Random restarts (default 3).

- seed:

  Optional RNG seed.

## Value

a `tucker_synergy`: `spatial` (`muscle x P`), `temporal` (`time x Q`),
`trial` (`trial x R`), `core` (`P x Q x R`), `vaf`, `iterations`.

## References

Kim YD, Choi S (2007) non-negative Tucker; Tucker LR (1966).

## See also

[`muscleSynergyTensor()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyTensor.md),
[`spaceByTimeSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/spaceByTimeSynergy.md)

## Examples

``` r
set.seed(1)
A <- matrix(stats::runif(6*2),6,2); B <- matrix(stats::runif(40*2),40,2); C <- matrix(stats::runif(8*2),8,2)
G <- array(stats::runif(2*2*2), c(2,2,2))
X <- reconstructTucker(G, A, B, C)
tuckerSynergy(X, ranks = c(2, 2, 2), n_restart = 1)$vaf
#> [1] 0.9994916
```
