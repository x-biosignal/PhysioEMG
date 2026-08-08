# EMG Envelope Extraction

Extracts the amplitude envelope from EMG signals using various methods.
The RMS method computes root mean square over a sliding window. The
Hilbert method uses the analytic signal via the Hilbert transform. The
lowpass method rectifies the signal and applies a moving-average lowpass
filter.

## Usage

``` r
emgEnvelope(
  x,
  method = c("rms", "hilbert", "lowpass"),
  window_ms = 50,
  cutoff = 6,
  assay_name = NULL,
  output_assay = "envelope"
)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- method:

  Envelope method: "rms" (root mean square), "hilbert" (Hilbert
  transform), or "lowpass" (rectification + lowpass filter).

- window_ms:

  Window size in milliseconds for RMS method (default: 50).

- cutoff:

  Cutoff frequency in Hz for lowpass method (default: 6).

- assay_name:

  Input assay name (default: first assay).

- output_assay:

  Output assay name (default: "envelope").

## Value

A PhysioExperiment object with an additional assay named `output_assay`
containing the amplitude envelope. The envelope matrix has the same
dimensions as the input (time x channels) with non-negative values
representing instantaneous signal amplitude.

## References

De Luca, C.J. (1997). "The use of surface electromyography in
biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
doi:10.1123/jab.13.2.135

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgAmplitudeNormalize()`](https://x-biosignal.github.io/PhysioEMG/reference/emgAmplitudeNormalize.md)
for normalizing envelope values,
[`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
for onset detection from envelope data,
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for fatigue analysis using spectral features
