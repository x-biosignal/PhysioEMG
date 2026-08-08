# Create EMG PhysioExperiment with Known Contraction

Generates a synthetic single-channel EMG `PhysioExperiment` with a
clearly defined contraction period. The contraction region has
substantially higher amplitude than the baseline, making it suitable for
testing onset detection algorithms.

## Usage

``` r
make_emg_contraction(
  n_time = 5000,
  sr = 1000,
  contraction_start = 0.3,
  contraction_end = 0.7,
  baseline_sd = 0.01,
  contraction_sd = 0.5
)
```

## Arguments

- n_time:

  Number of time points (default: 5000).

- sr:

  Sampling rate in Hz (default: 1000).

- contraction_start:

  Proportion of signal where contraction begins (default: 0.3).

- contraction_end:

  Proportion of signal where contraction ends (default: 0.7).

- baseline_sd:

  Standard deviation of baseline noise (default: 0.01).

- contraction_sd:

  Standard deviation of contraction activity (default: 0.5).

## Value

A `PhysioExperiment` object with a single `"raw"` assay containing a
one-channel EMG signal with a known contraction window. The contraction
region is defined by `contraction_start` and `contraction_end` as
proportions of total signal length.

## References

Hodges, P.W. & Bui, B.H. (1996). "A comparison of computer-based methods
for the determination of onset of muscle contraction using
electromyography." Electroencephalography and Clinical Neurophysiology,
101(6), 511-519. doi:10.1016/S0921-884X(96)95190-5

## See also

[`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
for detecting the contraction onset,
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for extracting the amplitude envelope,
[`make_emg()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg.md)
for basic multi-channel EMG,
[`make_emg_fatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_fatigue.md)
for EMG with fatigue progression

## Examples

``` r
pe <- make_emg_contraction()
pe
#> class: PhysioExperiment
#> dim: 5000 x 1 
#> assays(1): raw
#> samplingRate: 1000 Hz
#> channels(1): EMG1
#> colData names(2): label, type
onset <- emgOnsetDetect(pe)
onset$onsets
#>   channel sample time_sec
#> 1       1   1489    1.488
```
