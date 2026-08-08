# EMG Spectral Moments

Computes spectral moments (M0, M1, M2) over sliding windows. M0 is total
power, M1 is the first spectral moment (related to mean frequency), and
M2 is the second moment (related to bandwidth). These can be combined to
derive the mean frequency (M1/M0) and spectral bandwidth.

## Usage

``` r
emgSpectralMoments(x, window_sec = 1, overlap = 0.5, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- window_sec:

  Analysis window in seconds (default: 1.0).

- overlap:

  Overlap fraction (default: 0.5).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per channel per window, containing columns:

- channel:

  Integer channel index.

- window:

  Integer window number (1-indexed).

- m0:

  Zeroth spectral moment (total power).

- m1:

  First spectral moment (frequency-weighted power).

- m2:

  Second spectral moment (frequency-squared-weighted power).

## References

De Luca, C.J. (1984). "Myoelectrical manifestations of localized
muscular fatigue in humans." Critical Reviews in Biomedical Engineering,
11(4), 251-279.

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for median and mean frequency tracking,
[`emgFatigueIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueIndex.md)
for summary fatigue metric,
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for time-domain amplitude analysis
