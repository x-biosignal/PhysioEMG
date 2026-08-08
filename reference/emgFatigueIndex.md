# EMG Fatigue Index

Computes a fatigue index as the ratio of final to initial median
frequency. Values less than 1 indicate fatigue (frequency decrease). The
signal is internally analyzed with
[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
using 0.5-second windows at 50 percent overlap, then the initial and
final portions are compared.

## Usage

``` r
emgFatigueIndex(x, initial_pct = 0.2, final_pct = 0.2, assay_name = NULL)
```

## Arguments

- x:

  A PhysioExperiment object with EMG data.

- initial_pct:

  Percentage of signal used for initial estimate (default: 0.2).

- final_pct:

  Percentage of signal used for final estimate (default: 0.2).

- assay_name:

  Input assay name (default: first assay).

## Value

A data.frame with one row per channel, containing columns:

- channel:

  Integer channel index.

- fatigue_index:

  Ratio of final to initial median frequency. Values less than 1
  indicate fatigue.

- initial_mdf:

  Mean median frequency (Hz) in the initial portion.

- final_mdf:

  Mean median frequency (Hz) in the final portion.

## References

De Luca, C.J. (1984). "Myoelectrical manifestations of localized
muscular fatigue in humans." Critical Reviews in Biomedical Engineering,
11(4), 251-279.

Merletti, R. & Parker, P.A. (2004). "Electromyography: Physiology,
Engineering, and Non-Invasive Applications." Wiley-IEEE Press.
doi:10.1002/0471678384

## See also

[`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
for detailed windowed fatigue tracking,
[`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
for spectral moment analysis,
[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for amplitude envelope extraction
