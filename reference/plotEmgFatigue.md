# Plot EMG Spectral-Fatigue Regression

Scatters the windowed median or mean frequency against time with the
fitted fatigue-regression line and slope annotation.

## Usage

``` r
plotEmgFatigue(
  x,
  feature = c("mdf", "mnf"),
  window_sec = 1,
  overlap = 0.5,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- feature:

  "mdf" (median frequency) or "mnf" (mean frequency).

- window_sec, overlap:

  Windowing passed to
  [`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md).

- assay_name:

  Input assay name (default: first assay).

## Value

A ggplot object (faceted by channel).

## References

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md),
[`emgFatigueSlope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueSlope.md)

## Examples

``` r
pe <- make_emg_fatigue()
plotEmgFatigue(pe)
```
