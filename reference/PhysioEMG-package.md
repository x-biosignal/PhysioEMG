# PhysioEMG: EMG Analysis Functions for PhysioExperiment Objects

Provides electromyography (EMG) analysis functions for PhysioExperiment
objects. Includes envelope extraction (RMS, Hilbert), muscle activation
onset detection (Hodges-Bui, Teager-Kaiser), fatigue analysis (median
frequency shift), and muscle synergy decomposition (NMF, PCA, ICA).

## See also

Useful links:

- <https://github.com/x-biosignal/PhysioEMG>

- <https://x-biosignal.r-universe.dev/PhysioEMG>

- <https://x-biosignal.github.io/PhysioEMG/>

- Report bugs at <https://github.com/x-biosignal/PhysioEMG/issues>

## Author

**Maintainer**: Yusuke Matsui <mail.to.matsui@gmail.com>

## Examples

``` r
# Simulate EMG, extract an RMS envelope, and decompose into synergies
pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
pe <- emgEnvelope(pe, method = "rms")
muscleSynergy(pe, n_synergies = 2, method = "nmf", seed = 1,
              assay_name = "envelope")$vaf
#> [1] 0.9852343
```
