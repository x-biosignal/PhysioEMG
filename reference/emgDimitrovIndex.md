# EMG Dimitrov Spectral Fatigue Index (FInsm5)

Computes spectral moments \\M(k) = \int f^k\\PSD(f)\\df\\ for \\k = -1,
0, \dots, 5\\ over sliding windows and the Dimitrov fatigue index
\\FInsm5 = M(-1)/M(5)\\ per window. FInsm5 rises steeply with fatigue
because spectral compression toward low frequencies simultaneously
increases the low-frequency moment \\M(-1)\\ and decreases the
high-frequency moment \\M(5)\\. The per-channel FInsm5-vs-time
regression slope is attached as the `"finsm5_slope"` attribute.

## Usage

``` r
emgDimitrovIndex(x, window_sec = 1, overlap = 0.5, assay_name = NULL)
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

A data.frame with one row per channel per window, containing columns
`channel`, `window`, `time_sec`, the moments `m_minus1`, `m0`, `m1`,
`m2`, `m3`, `m4`, `m5`, and `finsm5`. The per-channel FInsm5 regression
slope (index per minute) is available via
`attr(result, "finsm5_slope")`.

## References

Dimitrov, G.V. et al. (2006). "Muscle fatigue during dynamic
contractions assessed by new spectral indices." Medicine & Science in
Sports & Exercise, 38(11), 1971-1979.
doi:10.1249/01.mss.0000233794.31659.6d

## See also

[`emgFatigueSlope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueSlope.md)
for MDF/MNF regression slopes,
[`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
for the M0/M1/M2 moments,
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for median/mean frequency tracking

## Examples

``` r
sr <- 1000; n <- 6000
sig <- sin(2 * pi * 80 * seq_len(n) / sr) + rnorm(n, sd = 0.1)
pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                       samplingRate = sr)
di <- emgDimitrovIndex(pe, window_sec = 1)
head(di)
#>   channel window time_sec m_minus1       m0       m1      m2        m3
#> 1       1      1      0.0 3.193003 254.9597 21307.52 2044424 296705025
#> 2       1      2      0.5 3.293293 261.6253 21828.79 2075671 293183280
#> 3       1      3      1.0 3.273150 260.8329 21698.37 2034646 276126944
#> 4       1      4      1.5 3.199904 255.8276 21265.91 1995946 273688041
#> 5       1      5      2.0 3.239246 259.6633 21603.78 2033593 280097577
#> 6       1      6      2.5 3.246536 259.1679 21590.98 2043891 285489784
#>            m4           m5       finsm5
#> 1 78479797921 2.958729e+13 1.079181e-13
#> 2 75143755582 2.782835e+13 1.183431e-13
#> 3 67438546902 2.424274e+13 1.350156e-13
#> 4 68213770241 2.494587e+13 1.282739e-13
#> 5 69773955425 2.536914e+13 1.276845e-13
#> 6 72159097919 2.643225e+13 1.228248e-13
attr(di, "finsm5_slope")
#>   channel slope_per_min r_squared   p_value
#> 1       1 -1.349253e-13 0.1386665 0.2594059
```
