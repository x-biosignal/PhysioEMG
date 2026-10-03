# Reconstruct a Tucker model (helper)

Rebuilds the `muscle x time x trial` tensor from a Tucker core and
factors; exposed mainly for examples and testing.

## Usage

``` r
reconstructTucker(core, spatial, temporal, trial)
```

## Arguments

- core:

  A `P x Q x R` core array.

- spatial, temporal, trial:

  Factor matrices.

## Value

The reconstructed `[muscle, time, trial]` array.

## Examples

``` r
set.seed(1)
A <- matrix(stats::runif(6 * 2), 6, 2)
B <- matrix(stats::runif(40 * 2), 40, 2)
C <- matrix(stats::runif(8 * 2), 8, 2)
G <- array(stats::runif(2 * 2 * 2), c(2, 2, 2))
X <- reconstructTucker(G, A, B, C)
dim(X)  # 6 x 40 x 8
#> [1]  6 40  8
```
