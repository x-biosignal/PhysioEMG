# EMG Fatigue Regression Slope

Fits a linear regression of a spectral fatigue metric (median or mean
frequency) against time across sliding windows and reports the fatigue
slope per channel. A negative slope indicates myoelectric fatigue
(spectral compression toward lower frequencies).

## Usage

``` r
emgFatigueSlope(
  x,
  feature = c("mdf", "mnf"),
  normalize = TRUE,
  window_sec = 1,
  overlap = 0.5,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- feature:

  Metric to regress: "mdf" (median frequency) or "mnf" (mean frequency).

- normalize:

  If TRUE, also report the slope as a percentage of the initial
  (intercept) value per minute.

- window_sec:

  Analysis window in seconds (default: 1.0).

- overlap:

  Overlap fraction between windows (default: 0.5).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per channel, containing columns:

- channel:

  Integer channel index.

- slope_hz_per_min:

  Regression slope in Hz per minute (negative under fatigue).

- norm_slope_pct_per_min:

  Slope as a percentage of the initial value per minute (`NA` if
  `normalize = FALSE`).

- intercept_hz:

  Fitted value at time zero (initial frequency), Hz.

- r_squared:

  Coefficient of determination of the fit.

- p_value:

  Two-sided p-value for the slope.

## References

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for the underlying windowed frequency estimates,
[`emgDimitrovIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDimitrovIndex.md)
for the spectral-moment fatigue index,
[`emgFatigueIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueIndex.md)
for the initial/final ratio metric

## Examples

``` r
# decreasing-frequency (fatiguing) signal
sr <- 1000; n <- 8000; t <- seq_len(n) / sr
f <- 90 - 3 * t                      # 90 Hz falling to ~66 Hz
sig <- sin(2 * pi * cumsum(f) / sr) + rnorm(n, sd = 0.1)
pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                       samplingRate = sr)
emgFatigueSlope(pe, feature = "mdf")
#>   channel slope_hz_per_min norm_slope_pct_per_min intercept_hz r_squared
#> 1       1        -179.5714              -202.8101     88.54167 0.9969317
#>        p_value
#> 1 1.004725e-17
```
