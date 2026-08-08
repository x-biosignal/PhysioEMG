# Create EMG PhysioExperiment with Fatigue Progression

Generates a synthetic single-channel EMG `PhysioExperiment` that
simulates a fatiguing isometric contraction. The median frequency of the
signal decreases progressively over time (from 80 Hz to approximately 44
Hz across 10 segments), mimicking the spectral compression
characteristic of muscle fatigue.

## Usage

``` r
make_emg_fatigue(n_time = 10000, sr = 1000)
```

## Arguments

- n_time:

  Number of time points (default: 10000).

- sr:

  Sampling rate in Hz (default: 1000).

## Value

A `PhysioExperiment` object with a single `"raw"` assay containing a
one-channel EMG signal exhibiting progressive median frequency decrease
across 10 equal-length segments.

## References

De Luca, C.J. (1997). "The use of surface electromyography in
biomechanics." Journal of Applied Biomechanics, 13(2), 135-163.
doi:10.1123/jab.13.2.135

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for tracking median frequency over time,
[`emgFatigueIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueIndex.md)
for computing a summary fatigue metric,
[`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
for spectral moment analysis,
[`make_emg()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg.md)
for basic multi-channel EMG

## Examples

``` r
pe <- make_emg_fatigue()
pe
#> class: PhysioExperiment
#> dim: 10000 x 1 
#> assays(1): raw
#> samplingRate: 1000 Hz
#> channels(1): EMG1
#> colData names(2): label, type
fatigue <- emgFatigue(pe)
head(fatigue)
#>   channel window time_sec median_freq mean_freq rms_amplitude
#> 1       1      1      0.0          80  156.5962     0.2946107
#> 2       1      2      0.5          80  154.2770     0.2881832
#> 3       1      3      1.0          76  155.5747     0.2885875
#> 4       1      4      1.5          76  155.9015     0.2935949
#> 5       1      5      2.0          72  153.4184     0.3045947
#> 6       1      6      2.5          72  152.0734     0.3074131
```
