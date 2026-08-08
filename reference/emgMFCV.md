# Muscle-Fiber Conduction Velocity (MFCV)

Estimates muscle-fiber conduction velocity from the propagation delay of
the EMG signal between pairs of electrodes aligned along the fiber
direction. The delay between the two channels of each pair is estimated
either by cross-correlation (with parabolic sub-sample interpolation) or
from the slope of the cross-spectrum phase, and the velocity is the
inter-electrode distance divided by the delay.

## Usage

``` r
emgMFCV(
  x,
  electrode_pairs,
  ied_mm,
  method = c("xcorr", "phase"),
  band = c(20, 250),
  max_lag_ms = 25,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- electrode_pairs:

  Channel pairs aligned along the muscle fibres, given as a list of
  length-2 integer vectors `list(c(prox, dist), ...)` or as a two-column
  matrix (one pair per row; column 1 proximal, column 2 distal).

- ied_mm:

  Inter-electrode distance in millimetres. A single value applied to all
  pairs, or one value per pair.

- method:

  Delay-estimation method: "xcorr" (cross-correlation, default) or
  "phase" (cross-spectrum phase slope).

- band:

  Frequency band `c(fmin, fmax)` in Hz used by the "phase" method
  (default `c(20, 250)`).

- max_lag_ms:

  Maximum delay searched by the "xcorr" method, in milliseconds
  (default: 25).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per electrode pair, containing columns:

- pair:

  Integer pair index.

- ch1, ch2:

  Proximal and distal channel indices.

- delay_ms:

  Estimated propagation delay in milliseconds.

- velocity_m_s:

  Conduction velocity in metres per second (`ied_mm` / delay).

- quality:

  Peak cross-correlation ("xcorr") or phase-fit R^2 ("phase").

## References

Farina, D. & Merletti, R. (2000). "Comparison of algorithms for
estimation of EMG variables during voluntary isometric contractions."
Journal of Electromyography and Kinesiology, 10(5), 337-349.
doi:10.1016/S1050-6411(00)00025-0

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigueSlope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueSlope.md)
and
[`emgDimitrovIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDimitrovIndex.md)
for spectral fatigue metrics,
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for median/mean frequency tracking

## Examples

``` r
# two electrodes 10 mm apart; distal is the proximal signal delayed 3 samples
sr <- 2000
set.seed(1)
prox <- as.numeric(stats::filter(rnorm(4000), rep(1, 5), sides = 2))
prox[is.na(prox)] <- 0
dist <- c(rep(0, 3), prox[seq_len(length(prox) - 3)])
pe <- PhysioExperiment(assays = list(raw = cbind(prox, dist)),
                       samplingRate = sr)
emgMFCV(pe, electrode_pairs = list(c(1, 2)), ied_mm = 10)
#>   pair ch1 ch2 delay_ms velocity_m_s   quality
#> 1    1   1   2 1.500021     6.666574 0.9999945
```
