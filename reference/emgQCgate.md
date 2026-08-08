# EMG Quality-Control Gate

Applies pass/fail thresholds to the
[`emgQualityCheck()`](https://x-biosignal.github.io/PhysioEMG/reference/emgQualityCheck.md)
metrics and reports, per channel, whether the signal is acceptable and,
if not, why.

## Usage

``` r
emgQCgate(
  x,
  min_snr = 6,
  max_powerline = 0.1,
  max_clipping = 1,
  reject_ecg = TRUE,
  assay_name = NULL,
  ...
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- min_snr:

  Minimum acceptable SNR in dB (default: 6).

- max_powerline:

  Maximum acceptable power-line ratio (default: 0.10).

- max_clipping:

  Maximum acceptable clipping percentage (default: 1).

- reject_ecg:

  If TRUE (default), ECG-contaminated channels fail the gate.

- assay_name:

  Input assay name (default: first assay).

- ...:

  Passed to
  [`emgQualityCheck()`](https://x-biosignal.github.io/PhysioEMG/reference/emgQualityCheck.md).

## Value

A list with `pass` (TRUE if every channel passes) and `channels`, a
data.frame with columns `channel`, `pass` and `reasons` (a
comma-separated string of failed criteria, or "").

## References

McManus, L., De Vito, G. & Lowery, M.M. (2020). "Analysis and biophysics
of surface EMG for physiotherapists and kinesiologists." Frontiers in
Neurology, 11, 576729. doi:10.3389/fneur.2020.576729

## See also

[`emgQualityCheck()`](https://x-biosignal.github.io/PhysioEMG/reference/emgQualityCheck.md),
[`emgRemoveECG()`](https://x-biosignal.github.io/PhysioEMG/reference/emgRemoveECG.md)

## Examples

``` r
set.seed(1)
sig <- rnorm(4000, sd = 0.02); sig[1500:2500] <- rnorm(1001, sd = 0.4)
pe <- PhysioExperiment(assays = list(raw = matrix(sig, ncol = 1)),
                       samplingRate = 1000)
emgQCgate(pe)$pass
#> [1] TRUE
```
