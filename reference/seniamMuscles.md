# SENIAM Muscle Catalog

Returns the bundled catalog of SENIAM-recommended surface-EMG electrode
locations (Hermens et al. 2000): for each muscle the recommended sensor
location, the electrode orientation (fibre direction) and the
recommended inter-electrode distance (mm).

## Usage

``` r
seniamMuscles()
```

## Value

A data.frame with columns `muscle`, `location`, `orientation` and
`ied_mm`.

## References

Hermens, H.J., Freriks, B., Disselhorst-Klug, C. & Rau, G. (2000).
"Development of recommendations for SEMG sensors and sensor placement
procedures." Journal of Electromyography and Kinesiology, 10(5),
361-374. doi:10.1016/S1050-6411(00)00027-4

## See also

[`setEMGElectrode()`](https://x-biosignal.github.io/PhysioEMG/reference/setEMGElectrode.md)
to attach electrode metadata to a PhysioExperiment,
[`getEMGElectrode()`](https://x-biosignal.github.io/PhysioEMG/reference/getEMGElectrode.md)
to retrieve it

## Examples

``` r
head(seniamMuscles())
#>                    muscle
#> 1    Trapezius descendens
#> 2 Trapezius transversalis
#> 3     Trapezius ascendens
#> 4     Deltoideus anterior
#> 5       Deltoideus medius
#> 6    Deltoideus posterior
#>                                                               location
#> 1       At 1/2 on the line from the acromion to the C7 spinous process
#> 2 At 1/2 between the medial border of the scapula and the spine, at T3
#> 3               At 2/3 on the line from the spine (T8) to the acromion
#> 4                 One finger width distal and anterior to the acromion
#> 5        From the acromion to the lateral epicondyle at greatest bulge
#> 6             About two finger widths behind the angle of the acromion
#>                  orientation ied_mm
#> 1          Towards the spine     20
#> 2          Towards the spine     20
#> 3       Towards the acromion     20
#> 4          Towards the thumb     20
#> 5 Vertically over the muscle     20
#> 6  Towards the little finger     20
```
