# Reconstruct Data from Synergies

Reconstructs EMG data using a subset of synergies.

## Usage

``` r
synergyReconstruct(synergy_result, n_synergies)
```

## Arguments

- synergy_result:

  Result from
  [`muscleSynergy`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md).

- n_synergies:

  Number of synergies to use for reconstruction.

## Value

A list with:

- `reconstructed`: Reconstructed data matrix

- `vaf`: VAF of the reconstruction

## References

De Luca, C.J. (1997). "The use of surface electromyography in
biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
doi:10.1123/jab.13.2.135

## See also

[`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
for computing the initial decomposition,
[`synergyCompare()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyCompare.md)
for comparing synergy solutions

## Examples

``` r
pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
syn <- muscleSynergy(pe, n_synergies = 3, method = "nmf", seed = 1)
rec <- synergyReconstruct(syn, n_synergies = 2)
rec$vaf
#> [1] 0.5195073
```
