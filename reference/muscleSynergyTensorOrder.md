# Select the number of tensor synergies by a VAF threshold

Fits non-negative PARAFAC for increasing rank and returns the smallest
number of synergies whose reconstruction VAF reaches `vaf_threshold`.

## Usage

``` r
muscleSynergyTensorOrder(tensor, max_synergies = 6, vaf_threshold = 0.9, ...)
```

## Arguments

- tensor:

  A non-negative `muscle x time x trial` array.

- max_synergies:

  Largest rank to try (default 6).

- vaf_threshold:

  VAF target (default 0.90).

- ...:

  Passed to
  [`muscleSynergyTensor()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyTensor.md).

## Value

A list with `n_synergies` (selected), `vaf` (per rank), and `fit` (the
selected `"synergy_tensor"`).

## Examples

``` r
set.seed(1)
M <- 6; Tt <- 40; K <- 8; R <- 2
A0 <- matrix(runif(M * R), M, R); B0 <- matrix(runif(Tt * R), Tt, R)
C0 <- matrix(runif(K * R), K, R)
X <- array(0, c(M, Tt, K))
for (r in 1:R) for (k in 1:K) X[, , k] <- X[, , k] +
  (A0[, r] %o% B0[, r]) * C0[k, r]
ord <- muscleSynergyTensorOrder(X, max_synergies = 3, seed = 1)
ord$n_synergies
#> [1] 1
```
