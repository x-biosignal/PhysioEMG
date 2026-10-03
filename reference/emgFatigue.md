# EMG Fatigue Analysis

Tracks median and mean frequency over time to assess muscle fatigue.
Decreasing median frequency indicates fatigue due to reduced motor unit
conduction velocity. The signal is divided into overlapping windows and
the power spectral density is computed via FFT for each window.

## Usage

``` r
emgFatigue(x, window_sec = 1, overlap = 0.5, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- window_sec:

  Analysis window in seconds (default: 1.0).

- overlap:

  Overlap fraction between windows (default: 0.5).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per channel per window, containing columns:

- channel:

  Integer channel index.

- window:

  Integer window number (1-indexed).

- time_sec:

  Start time of the window in seconds.

- median_freq:

  Median frequency (Hz) at which 50 percent of the spectral power is
  below.

- mean_freq:

  Power-weighted mean frequency (Hz).

- rms_amplitude:

  Root mean square amplitude of the window.

## References

De Luca, C.J. (1984). "Myoelectrical manifestations of localized
muscular fatigue in humans." Critical Reviews in Biomedical Engineering,
11(4), 251-279.

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigueIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueIndex.md)
for a summary fatigue metric,
[`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
for spectral moment analysis,
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for amplitude envelope extraction

## Examples

``` r
pe <- make_emg_fatigue(n_time = 10000, sr = 1000)
fat <- emgFatigue(pe, window_sec = 1.0)
head(fat)
#>   channel window time_sec median_freq mean_freq rms_amplitude
#> 1       1      1      0.0          80  157.8732     0.2815291
#> 2       1      2      0.5          80  157.0243     0.2885262
#> 3       1      3      1.0          76  153.3099     0.2915280
#> 4       1      4      1.5          76  152.7389     0.2865469
#> 5       1      5      2.0          72  150.8788     0.2903286
#> 6       1      6      2.5          72  150.1204     0.2960479
```
