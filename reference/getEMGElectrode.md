# Get SENIAM Electrode Metadata

Retrieves the SENIAM electrode/muscle metadata previously attached with
[`setEMGElectrode()`](https://x-biosignal.github.io/PhysioEMG/reference/setEMGElectrode.md).

## Usage

``` r
getEMGElectrode(x)
```

## Arguments

- x:

  A PhysioExperiment object.

## Value

A data.frame with one row per channel and columns `channel`, `muscle`,
`placement`, `orientation`, `ied_mm`, `reference` and `side` (`NA` where
unset).

## See also

[`setEMGElectrode()`](https://x-biosignal.github.io/PhysioEMG/reference/setEMGElectrode.md),
[`seniamMuscles()`](https://x-biosignal.github.io/PhysioEMG/reference/seniamMuscles.md)

## Examples

``` r
pe <- PhysioExperiment(assays = list(raw = matrix(rnorm(100), ncol = 1)),
                       samplingRate = 1000)
pe <- setEMGElectrode(pe, muscle = "Tibialis anterior")
getEMGElectrode(pe)
#>   channel            muscle
#> 1       1 Tibialis anterior
#>                                                                              placement
#> 1 At 1/3 on the line between the tip of the fibula and the tip of the medial malleolus
#>       orientation ied_mm reference side
#> 1 Along the tibia     20      <NA> <NA>
```
