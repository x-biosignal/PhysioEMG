# Normalize EMG Amplitude

Normalizes EMG amplitude data by a per-channel reference value. `"mvc"`
divides by the per-channel maximum of a maximum-voluntary- contraction
trial (percentage-of-MVC). `"peak"` divides by the within-trial peak so
each channel ranges from 0 to 1. `"rvc"` divides by the per-channel mean
amplitude of a sub-maximal reference voluntary contraction trial
(percentage-of-RVC). `"dynamic_peak"` divides by a centered moving
maximum, tracking a time-varying peak within the trial.

## Usage

``` r
emgAmplitudeNormalize(
  x,
  method = c("mvc", "peak", "rvc", "dynamic_peak"),
  mvc_data = NULL,
  rvc_data = NULL,
  rvc_window = NULL,
  assay_name = NULL,
  output_assay = "normalized"
)
```

## Arguments

- x:

  A PhysioExperiment object (amplitude data).

- method:

  Normalization method: "mvc" (maximum voluntary contraction), "peak"
  (within-trial peak), "rvc" (sub-maximal reference voluntary
  contraction), or "dynamic_peak" (centered moving maximum).

- mvc_data:

  A PhysioExperiment containing MVC trial data (required for "mvc").
  Must have the same number of channels as `x`.

- rvc_data:

  A PhysioExperiment containing the sub-maximal reference-task data
  (required for "rvc"). Must have the same number of channels as `x`.

- rvc_window:

  Optional window. For "rvc", a length-2 numeric `c(start, end)` in
  seconds selecting the portion of `rvc_data` over which the mean
  reference amplitude is computed (default: the whole reference trial).
  For "dynamic_peak", the moving-maximum window length in seconds
  (default: 0.5).

- assay_name:

  Assay to normalize (default: first assay).

- output_assay:

  Output assay name (default: "normalized").

## Value

A PhysioExperiment object with an additional assay named `output_assay`
containing normalized amplitude values. For "peak" and "dynamic_peak",
values are scaled to a peak. For "mvc"/"rvc", values are proportions of
the reference (multiply by 100 for percentage-of-MVC /
percentage-of-RVC).

## Details

Input is assumed to be amplitude data (a rectified signal or an
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
output). MVC and RVC normalization are scale-invariant to a gain applied
to both the signal and its reference trial.

## References

De Luca, C.J. (1997). "The use of surface electromyography in
biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
doi:10.1123/jab.13.2.135

Yang, J.F. & Winter, D.A. (1984). "Electromyographic amplitude
normalization methods: improving their sensitivity as diagnostic tools
in gait analysis." Archives of Physical Medicine and Rehabilitation,
65(9), 517-521.

Burden, A. & Bartlett, R. (1999). "Normalisation of EMG amplitude: an
evaluation and comparison of old and new methods." Medical Engineering &
Physics, 21(4), 247-257. doi:10.1016/S1350-4533(99)00054-5

## See also

[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for computing amplitude envelopes prior to normalization,
[`emgAmplitudeFeatures()`](https://x-biosignal.github.io/PhysioEMG/reference/emgAmplitudeFeatures.md)
for windowed amplitude features,
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for fatigue analysis,
[`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
for muscle activation onset detection

## Examples

``` r
set.seed(1)
amp <- matrix(abs(rnorm(1000 * 2, sd = 0.3)), nrow = 1000, ncol = 2)
pe <- PhysioExperiment(
  assays = list(raw = amp),
  colData = S4Vectors::DataFrame(label = c("EMG1", "EMG2"),
                                 type = c("EMG", "EMG")),
  samplingRate = 1000)
# percentage-of-RVC using a sub-maximal reference (here, the trial itself)
pe_rvc <- emgAmplitudeNormalize(pe, method = "rvc", rvc_data = pe)
colMeans(SummarizedExperiment::assay(pe_rvc, "normalized"))  # ~ 1
#> [1] 1 1
```
