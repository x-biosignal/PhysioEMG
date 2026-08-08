# Remove ECG Contamination from EMG

Removes cardiac (ECG/QRS) contamination from surface EMG using a
high-pass filter, QRS template subtraction, or gating (blanking) of the
QRS complexes.

## Usage

``` r
emgRemoveECG(
  x,
  method = c("highpass", "template", "gating"),
  ecg_channel = NULL,
  hp_cutoff = 30,
  qrs_ms = 60,
  assay_name = NULL,
  output_assay = "ecg_removed"
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- method:

  Removal method: "highpass" (default; FFT high-pass at `hp_cutoff`),
  "template" (average-QRS template subtraction), or "gating"
  (linear-interpolation blanking of QRS windows).

- ecg_channel:

  Optional channel index of a dedicated ECG reference used for R-peak
  detection in the "template"/"gating" methods; if NULL, R-peaks are
  detected from each EMG channel itself.

- hp_cutoff:

  High-pass cutoff in Hz for the "highpass" method (default: 30).

- qrs_ms:

  Half-width in ms of the QRS window for "template"/"gating" (default:
  60).

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: "ecg_removed").

## Value

A PhysioExperiment object with an added assay `output_assay` containing
the cleaned signal.

## References

Drake, J.D.M. & Callaghan, J.P. (2006). "Elimination of
electrocardiogram contamination from electromyogram signals: An
evaluation of currently used removal techniques." Journal of
Electromyography and Kinesiology, 16(2), 175-187.
doi:10.1016/j.jelekin.2005.07.003

Willigenburg, N.W., Daffertshofer, A., Kingma, I. & van Dieen, J.H.
(2012). "Removing ECG contamination from EMG recordings." Journal of
Electromyography and Kinesiology, 22(3), 485-493.
doi:10.1016/j.jelekin.2012.01.001

## See also

[`emgQualityCheck()`](https://x-biosignal.github.io/PhysioEMG/reference/emgQualityCheck.md)
to detect contamination,
[`emgQCgate()`](https://x-biosignal.github.io/PhysioEMG/reference/emgQCgate.md)

## Examples

``` r
set.seed(1)
sr <- 1000; n <- 5000
emg <- rnorm(n, sd = 0.1)
qrs <- numeric(n)
for (p in seq(500, n - 200, by = 1000)) qrs[p + -15:15] <- 3 * exp(-(-15:15)^2 / 20)
pe <- PhysioExperiment(assays = list(raw = matrix(emg + qrs, ncol = 1)),
                       samplingRate = sr)
clean <- emgRemoveECG(pe, method = "highpass")
```
