# Create a Basic EMG PhysioExperiment

Generates a synthetic multi-channel EMG `PhysioExperiment` object with
simulated muscle activity. Each channel contains baseline noise with a
burst of higher-amplitude activity in the middle 40 percent of the
signal (from 30 to 70 percent), mimicking a typical voluntary
contraction.

## Usage

``` r
make_emg(n_time = 2000, n_channels = 4, sr = 1000)
```

## Arguments

- n_time:

  Number of time points (default: 2000).

- n_channels:

  Number of EMG channels (default: 4).

- sr:

  Sampling rate in Hz (default: 1000).

## Value

A `PhysioExperiment` object with a single `"raw"` assay containing
simulated EMG data (time x channels matrix), channel metadata in
`colData`, and the specified sampling rate.

## References

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`make_emg_contraction()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_contraction.md)
for EMG with a defined contraction window,
[`make_emg_fatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_fatigue.md)
for EMG with fatigue progression,
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for extracting amplitude envelopes

## Examples

``` r
pe <- make_emg()
pe
#> class: PhysioExperiment
#> dim: 2000 x 4 
#> assays(1): raw
#> samplingRate: 1000 Hz
#> channels(4): EMG1, EMG2, EMG3, EMG4
#> colData names(2): label, type
dim(SummarizedExperiment::assay(pe, "raw"))
#> [1] 2000    4
```
