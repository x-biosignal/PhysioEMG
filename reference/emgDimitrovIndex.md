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
#> 1       1      1      0.0 3.175804 254.6519 21325.03 2056167 299493123
#> 2       1      2      0.5 3.159098 252.3385 21070.24 2009293 283962046
#> 3       1      3      1.0 3.193799 254.9841 21316.86 2055038 301903375
#> 4       1      4      1.5 3.197576 256.2773 21477.77 2085573 311515124
#> 5       1      5      2.0 3.151104 252.5249 21052.55 1993245 277620319
#> 6       1      6      2.5 3.133642 250.1177 20883.55 1989383 281089669
#>            m4           m5       finsm5
#> 1 78405813142 2.894247e+13 1.097282e-13
#> 2 71952596330 2.614338e+13 1.208374e-13
#> 3 80725383077 3.051720e+13 1.046557e-13
#> 4 84513358334 3.210148e+13 9.960836e-14
#> 5 69773757053 2.544662e+13 1.238319e-13
#> 6 71707359747 2.639359e+13 1.187274e-13
attr(di, "finsm5_slope")
#>   channel slope_per_min r_squared    p_value
#> 1       1  1.611849e-13 0.2895953 0.08769738
```
