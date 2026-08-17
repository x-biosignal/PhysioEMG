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
