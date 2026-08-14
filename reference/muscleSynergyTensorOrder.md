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
