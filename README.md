# PhysioEMG <img src="man/figures/logo.png" align="right" height="139" alt="PhysioEMG logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/x-biosignal/PhysioEMG/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/x-biosignal/PhysioEMG/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/PhysioEMG)](https://CRAN.R-project.org/package=PhysioEMG)
[![r-universe](https://x-biosignal.r-universe.dev/badges/PhysioEMG)](https://x-biosignal.r-universe.dev/PhysioEMG)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**EMG Analysis Functions for PhysioExperiment Objects**

PhysioEMG provides 21 exported functions for electromyography (EMG) analysis, built on top of PhysioCore. It covers the complete EMG analysis pipeline from signal conditioning through clinical interpretation: envelope extraction and amplitude normalization, spectral analysis, muscle activation onset detection, fatigue monitoring, muscle synergy decomposition, and inter-muscular connectivity network analysis -- all operating directly on `PhysioExperiment` objects.

## Installation

You can install PhysioEMG from [r-universe](https://x-biosignal.r-universe.dev):

```r
install.packages("PhysioEMG",
  repos = c("https://x-biosignal.r-universe.dev", "https://cloud.r-project.org"))
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("x-biosignal/PhysioEMG")
```

## Quick Start

```r
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

- `emgEnvelope()` -- extract signal envelope using RMS (sliding window), Hilbert transform, or lowpass rectification
- `emgAmplitudeNormalize()` -- normalize amplitude to maximum voluntary contraction (MVC) or peak value

### Spectral Analysis

Frequency-domain characterization of EMG signals:

- `emgSpectralMoments()` -- compute spectral moments (mean frequency, median frequency, bandwidth) over sliding windows

### Onset Detection

Automatic identification of muscle activation timing:

- `emgOnsetDetect()` -- detect muscle activation onsets with two algorithms:
  - **Hodges-Bui:** threshold-based detection on the rectified/smoothed signal
  - **Teager-Kaiser:** energy operator for improved sensitivity to rapid onsets

### Fatigue Analysis

Monitor neuromuscular fatigue during sustained or repeated contractions:

- `emgFatigue()` -- track median frequency shift over time in sliding windows (progressive decrease indicates fatigue)
- `emgFatigueIndex()` -- compute fatigue index by comparing spectral properties between initial and final contraction segments

### Muscle Synergy Decomposition

Extract coordinated muscle activation patterns underlying motor control:

- `muscleSynergy()` -- decompose multi-channel EMG into synergies using:
  - **NMF:** non-negative matrix factorization (physiologically interpretable, non-negative weights)
  - **PCA:** principal component analysis (orthogonal decomposition)
  - **ICA:** independent component analysis (statistically independent sources)
- `synergyReconstruct()` -- reconstruct EMG signals from a reduced set of synergies (assess reconstruction quality)
- `synergyCompare()` -- compare synergy structures between conditions, sessions, or subjects using similarity metrics

### Inter-Muscular Network Analysis

Characterize functional connectivity and coordination between muscles:

- `emgCoherenceNetwork()` -- magnitude-squared coherence networks within specified frequency bands
- `emgWPLINetwork()` -- weighted phase lag index networks (robust to volume conduction artifacts)
- `emgPartialCoherenceNetwork()` -- partial coherence networks controlling for common input effects
- `emgDirectedGCNetwork()` -- directed Granger causality networks for causal inter-muscular coupling
- `emgDynamicWaveletNetwork()` -- time-varying connectivity using wavelet coherence (track coordination changes during movement)
- `emgCoordinationStructure()` -- extract network topology metrics (modularity, hub muscles, clustering coefficient)
- `emgInterpretNetworkKG()` -- interpret network results using anatomical and functional knowledge graphs

### Simulated Data Generators

Ready-to-use data for testing, demonstration, and teaching:

- `make_emg()` -- multi-channel EMG with realistic burst patterns
- `make_emg_contraction()` -- EMG with controlled contraction-relaxation cycles
- `make_emg_fatigue()` -- EMG with progressive fatigue characteristics (spectral shift)

## Dependencies

- **R** (>= 4.2)
- **[PhysioCore](https://github.com/x-biosignal/PhysioCore)**
- **SummarizedExperiment**
- **stats**

## PhysioExperiment Ecosystem

PhysioEMG is the EMG analysis layer of the PhysioExperiment ecosystem, a suite of R packages for multi-modal physiological signal analysis:

| Package | Description |
|---------|-------------|
| [PhysioCore](https://github.com/x-biosignal/PhysioCore) | Core data structures and accessors |
| [PhysioIO](https://github.com/x-biosignal/PhysioIO) | File I/O (EDF, HDF5, BIDS, CSV, MAT) |
| [PhysioPreprocess](https://github.com/x-biosignal/PhysioPreprocess) | Preprocessing (filters, ICA, resampling) |
| [PhysioAnalysis](https://github.com/x-biosignal/PhysioAnalysis) | Analysis and visualization |
| [PhysioEEG](https://github.com/x-biosignal/PhysioEEG) | EEG analysis (ICA, ERP, source, BCI, sleep) |
| **PhysioEMG** | EMG analysis (synergy, fatigue, onset) |
| [PhysioECG](https://github.com/x-biosignal/PhysioECG) | ECG and HRV analysis |

Visit the [r-universe page](https://x-biosignal.r-universe.dev) to browse all available packages.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Author

Yusuke Matsui
