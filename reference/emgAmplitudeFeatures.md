# EMG Amplitude Features (ARV, MAV, iEMG, RMS)

Computes time-domain amplitude features of an EMG signal over a sliding
window. The signal is divided into overlapping windows and, for each
window and channel, the requested features are computed:

- `arv`:

  Average rectified value, \\\frac{1}{N}\sum \|x_i\|\\.

- `mav`:

  Mean absolute value, \\\frac{1}{N}\sum \|x_i\|\\ (identical to ARV for
  an unweighted window).

- `iemg`:

  Integrated EMG, \\\sum \|x_i\|\\\Delta t\\ with \\\Delta t = 1 /
  f_s\\; the time integral of the rectified signal, in signal-units
  times seconds.

- `rms`:

  Root mean square, \\\sqrt{\frac{1}{N}\sum x_i^2}\\.

The window geometry (window and step sizes in samples) matches
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
and
[`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
for identical `window_sec` and `overlap`, so amplitude and spectral
features align window-for-window.

## Usage

``` r
emgAmplitudeFeatures(
  x,
  features = c("arv", "mav", "iemg", "rms"),
  window_sec = 1,
  overlap = 0.5,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- features:

  Character vector selecting which features to return; any of `"arv"`,
  `"mav"`, `"iemg"`, `"rms"` (default: all four).

- window_sec:

  Analysis window in seconds (default: 1.0).

- overlap:

  Overlap fraction between consecutive windows, in \\\[0, 1)\\ (default:
  0.5).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per analysis window per channel, with columns:

- channel:

  Integer channel index.

- window:

  Integer window index within the channel (1-based).

- time_sec:

  Window start time in seconds.

- arv, mav, iemg, rms:

  The requested amplitude features (only the columns named in `features`
  are present).

If no full window fits in the signal, a 0-row data.frame with these
columns is returned.

## References

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

Hermens, H.J. et al. (2000). "Development of recommendations for SEMG
sensors and sensor placement procedures (SENIAM)." Journal of
Electromyography and Kinesiology, 10(5), 361-374.
doi:10.1016/S1050-6411(00)00027-4

## See also

[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for a per-sample amplitude envelope,
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for windowed spectral fatigue features,
[`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
for spectral moments on the same window grid

## Examples

``` r
# 2 s of two-channel EMG at 1000 Hz
set.seed(1)
m <- matrix(rnorm(2000 * 2, sd = 0.3), nrow = 2000, ncol = 2)
pe <- PhysioExperiment(
  assays = list(raw = m),
  colData = S4Vectors::DataFrame(label = c("EMG1", "EMG2"),
                                 type = c("EMG", "EMG")),
  samplingRate = 1000)
feats <- emgAmplitudeFeatures(pe, window_sec = 0.5, overlap = 0.5)
head(feats)
#>   channel window time_sec       arv       mav      iemg       rms
#> 1       1      1     0.00 0.2395552 0.2395552 0.1197776 0.3033508
#> 2       1      2     0.25 0.2542875 0.2542875 0.1271437 0.3167806
#> 3       1      3     0.50 0.2545804 0.2545804 0.1272902 0.3171735
#> 4       1      4     0.75 0.2496032 0.2496032 0.1248016 0.3129471
#> 5       1      5     1.00 0.2424656 0.2424656 0.1212328 0.3025909
#> 6       1      6     1.25 0.2467168 0.2467168 0.1233584 0.3094300
```
