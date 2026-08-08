# PhysioEMG ![PhysioEMG logo](reference/figures/logo.png)

**EMG Analysis Functions for PhysioExperiment Objects**

PhysioEMG provides 21 exported functions for electromyography (EMG)
analysis, built on top of PhysioCore. It covers the complete EMG
analysis pipeline from signal conditioning through clinical
interpretation: envelope extraction and amplitude normalization,
spectral analysis, muscle activation onset detection, fatigue
monitoring, muscle synergy decomposition, and inter-muscular
connectivity network analysis – all operating directly on
`PhysioExperiment` objects.

## Installation

You can install PhysioEMG from
[r-universe](https://x-biosignal.r-universe.dev):

``` r

install.packages("PhysioEMG",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioEMG")
```

## Quick Start

``` r

library(PhysioEMG)

# Generate simulated EMG with muscle bursts
pe <- make_emg(n_time = 5000, n_channels = 4, sr = 1000)

# Detect muscle activation onsets
onsets <- emgOnsetDetect(pe, method = "hodges_bui")

# Extract RMS envelope and normalize
pe_env <- emgEnvelope(pe, method = "rms", window_ms = 50)
pe_norm <- emgAmplitudeNormalize(pe_env, method = "peak")

# Decompose into muscle synergies
syn <- muscleSynergy(pe, n_synergies = 3, method = "nmf")

# Analyze inter-muscular coordination network
net <- emgCoherenceNetwork(pe, freq_band = c(10, 50))
coord <- emgCoordinationStructure(net)
```

## Features

### Envelope Extraction and Amplitude Normalization

Signal conditioning for amplitude analysis:

- [`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
  – extract signal envelope using RMS (sliding window), Hilbert
  transform, or lowpass rectification
- [`emgAmplitudeNormalize()`](https://x-biosignal.github.io/PhysioEMG/reference/emgAmplitudeNormalize.md)
  – normalize amplitude to maximum voluntary contraction (MVC) or peak
  value

### Spectral Analysis

Frequency-domain characterization of EMG signals:

- [`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
  – compute spectral moments (mean frequency, median frequency,
  bandwidth) over sliding windows

### Onset Detection

Automatic identification of muscle activation timing:

- [`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
  – detect muscle activation onsets with two algorithms:
  - **Hodges-Bui:** threshold-based detection on the rectified/smoothed
    signal
  - **Teager-Kaiser:** energy operator for improved sensitivity to rapid
    onsets

### Fatigue Analysis

Monitor neuromuscular fatigue during sustained or repeated contractions:

- [`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
  – track median frequency shift over time in sliding windows
  (progressive decrease indicates fatigue)
- [`emgFatigueIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueIndex.md)
  – compute fatigue index by comparing spectral properties between
  initial and final contraction segments

### Muscle Synergy Decomposition

Extract coordinated muscle activation patterns underlying motor control:

- [`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
  – decompose multi-channel EMG into synergies using:
  - **NMF:** non-negative matrix factorization (physiologically
    interpretable, non-negative weights)
  - **PCA:** principal component analysis (orthogonal decomposition)
  - **ICA:** independent component analysis (statistically independent
    sources)
- [`synergyReconstruct()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyReconstruct.md)
  – reconstruct EMG signals from a reduced set of synergies (assess
  reconstruction quality)
- [`synergyCompare()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyCompare.md)
  – compare synergy structures between conditions, sessions, or subjects
  using similarity metrics

### Inter-Muscular Network Analysis

Characterize functional connectivity and coordination between muscles:

- [`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md)
  – magnitude-squared coherence networks within specified frequency
  bands
- [`emgWPLINetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgWPLINetwork.md)
  – weighted phase lag index networks (robust to volume conduction
  artifacts)
- [`emgPartialCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgPartialCoherenceNetwork.md)
  – partial coherence networks controlling for common input effects
- [`emgDirectedGCNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDirectedGCNetwork.md)
  – directed Granger causality networks for causal inter-muscular
  coupling
- [`emgDynamicWaveletNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDynamicWaveletNetwork.md)
  – time-varying connectivity using wavelet coherence (track
  coordination changes during movement)
- [`emgCoordinationStructure()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoordinationStructure.md)
  – extract network topology metrics (modularity, hub muscles,
  clustering coefficient)
- [`emgInterpretNetworkKG()`](https://x-biosignal.github.io/PhysioEMG/reference/emgInterpretNetworkKG.md)
  – interpret network results using anatomical and functional knowledge
  graphs

### Simulated Data Generators

Ready-to-use data for testing, demonstration, and teaching:

- [`make_emg()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg.md)
  – multi-channel EMG with realistic burst patterns
- [`make_emg_contraction()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_contraction.md)
  – EMG with controlled contraction-relaxation cycles
- [`make_emg_fatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_fatigue.md)
  – EMG with progressive fatigue characteristics (spectral shift)

## Dependencies

- **R** (\>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **stats**

## PhysioExperiment Ecosystem

PhysioEMG is the EMG analysis layer of the PhysioExperiment ecosystem, a
suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|----|----|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| [PhysioEEG](https://github.com/x-biosignal/PhysioEEG) | EEG analysis (ICA, ERP, source, BCI, sleep) |
| **PhysioEMG** | EMG analysis (synergy, fatigue, onset) |
| [PhysioECG](https://github.com/x-biosignal/PhysioECG) | ECG and HRV analysis |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to
browse all available packages.

## License

MIT License. See
[LICENSE](https://x-biosignal.github.io/PhysioEMG/LICENSE) for details.

## Author

Yusuke Matsui

## Governance & support

Part of the [Physio ecosystem](https://x-biosignal.r-universe.dev).
Community and policy documents live in the umbrella repository:

- [Code of
  Conduct](https://github.com/x-biosignal/PhysioExperiment/blob/main/CODE_OF_CONDUCT.md)
- [Contributing](https://github.com/x-biosignal/PhysioExperiment/blob/main/CONTRIBUTING.md)
- [Governance](https://github.com/x-biosignal/PhysioExperiment/blob/main/GOVERNANCE.md)
- [Support](https://github.com/x-biosignal/PhysioExperiment/blob/main/SUPPORT.md)
- [Security
  policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/SECURITY.md)
- [Deprecation & lifecycle
  policy](https://github.com/x-biosignal/PhysioExperiment/blob/main/DEPRECATION.md)
