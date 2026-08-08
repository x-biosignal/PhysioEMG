# Set SENIAM Electrode Metadata

Attaches standardized SENIAM electrode/muscle metadata to the channels
of a PhysioExperiment, written as `seniam_*` columns of `colData`.
Muscle names are validated against
[`seniamMuscles()`](https://x-biosignal.github.io/PhysioEMG/reference/seniamMuscles.md);
a recognized muscle also back-fills its recommended electrode location,
orientation and inter-electrode distance when those are not supplied.

## Usage

``` r
setEMGElectrode(
  x,
  muscle,
  placement = NULL,
  ied_mm = NULL,
  reference = NULL,
  side = NULL,
  channel = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG channels.

- muscle:

  Muscle name(s) (character). Length 1 (applied to all targeted
  channels) or one per targeted channel.

- placement:

  Optional electrode location text; back-filled from the catalog when
  `NULL`.

- ied_mm:

  Optional inter-electrode distance in mm; back-filled from the catalog
  (SENIAM default 20 mm) when `NULL`.

- reference:

  Optional reference-electrode description.

- side:

  Optional body side ("left"/"right"/"L"/"R"/...).

- channel:

  Channels to annotate: integer indices, channel labels, or `NULL`
  (default) for all channels.

## Value

The PhysioExperiment with `seniam_muscle`, `seniam_placement`,
`seniam_orientation`, `seniam_ied_mm`, `seniam_reference` and
`seniam_side` added to `colData`.

## References

Hermens, H.J. et al. (2000). "Development of recommendations for SEMG
sensors and sensor placement procedures." Journal of Electromyography
and Kinesiology, 10(5), 361-374. doi:10.1016/S1050-6411(00)00027-4

## See also

[`getEMGElectrode()`](https://x-biosignal.github.io/PhysioEMG/reference/getEMGElectrode.md),
[`seniamMuscles()`](https://x-biosignal.github.io/PhysioEMG/reference/seniamMuscles.md),
[`emgMFCV()`](https://x-biosignal.github.io/PhysioEMG/reference/emgMFCV.md)
which can use the stored inter-electrode distance

## Examples

``` r
pe <- PhysioExperiment(
  assays = list(raw = matrix(rnorm(200), ncol = 2)),
  colData = S4Vectors::DataFrame(label = c("ch1", "ch2")),
  samplingRate = 1000)
pe <- setEMGElectrode(pe, muscle = c("Biceps brachii", "Triceps brachii"),
                      side = "right")
getEMGElectrode(pe)
#>   channel          muscle
#> 1       1  Biceps brachii
#> 2       2 Triceps brachii
#>                                                                                       placement
#> 1       On the line between the medial acromion and the fossa cubit at 1/3 from the fossa cubit
#> 2 At 1/2 on the line between the posterior acromion and the olecranon at 2 finger widths medial
#>                                    orientation ied_mm reference  side
#> 1 In the direction of the line to the acromion     20      <NA> right
#> 2 In the direction of the line to the acromion     20      <NA> right
```
