# EMG Signal Quality Check

Computes per-channel surface-EMG quality metrics: signal-to-noise ratio,
resting baseline noise, power-line interference ratio at 50 and 60 Hz,
clipping/saturation percentage, and an ECG-contamination flag based on
the periodicity of the QRS-band envelope. When the PhysioECG package is
available the R-peak count it detects is reported as an additional
diagnostic.

## Usage

``` r
emgQualityCheck(
  x,
  ecg_score_threshold = 0.5,
  rms_window_ms = 100,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- ecg_score_threshold:

  QRS-band periodicity score above which a channel is flagged as
  ECG-contaminated (default: 0.5).

- rms_window_ms:

  Window (ms) for the moving-RMS used to estimate signal and rest levels
  (default: 100).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per channel:

- channel:

  Integer channel index.

- snr_db:

  Signal-to-noise ratio in dB (active vs rest RMS).

- baseline_noise:

  Resting baseline RMS (noise floor).

- powerline_50, powerline_60:

  Fraction of total power within +/-1 Hz of 50 and 60 Hz.

- powerline_ratio:

  The larger of `powerline_50`/`powerline_60`.

- clipping_pct:

  Percentage of samples at the signal extremes (ADC saturation).

- ecg_score:

  QRS-band periodicity score from 0 (none) to 1 (strong).

- ecg_contamination:

  TRUE if `ecg_score` exceeds the threshold.

## References

McManus, L., De Vito, G. & Lowery, M.M. (2020). "Analysis and biophysics
of surface EMG for physiotherapists and kinesiologists." Frontiers in
Neurology, 11, 576729. doi:10.3389/fneur.2020.576729

## See also

[`emgRemoveECG()`](https://x-biosignal.github.io/PhysioEMG/reference/emgRemoveECG.md)
to remove ECG contamination,
[`emgQCgate()`](https://x-biosignal.github.io/PhysioEMG/reference/emgQCgate.md)
for pass/fail gating,
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)

## Examples

``` r
set.seed(1)
sig <- rnorm(4000, sd = 0.02)
sig[1500:2500] <- rnorm(1001, sd = 0.4)          # a contraction
pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                       samplingRate = 1000)
emgQualityCheck(pe)
#>   channel   snr_db baseline_noise powerline_50 powerline_60 powerline_ratio
#> 1       1 26.55082     0.01887598  0.006956421  0.007947266     0.007947266
#>   clipping_pct ecg_score ecg_contamination ecg_rpeaks
#> 1        0.025 0.2413754             FALSE          0
```
