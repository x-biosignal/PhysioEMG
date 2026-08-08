# Compare Two Synergy Results

Computes pairwise correlation between synergy weight vectors from two
decompositions. Uses best-match pairing.

## Usage

``` r
synergyCompare(result1, result2)
```

## Arguments

- result1:

  First result from
  [`muscleSynergy`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md).

- result2:

  Second result from
  [`muscleSynergy`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md).

## Value

A data.frame with columns: synergy1, synergy2, correlation.

## References

De Luca, C.J. (1997). "The use of surface electromyography in
biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
doi:10.1123/jab.13.2.135

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
for computing synergy decompositions,
[`synergyReconstruct()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyReconstruct.md)
for reconstructing data from synergies
