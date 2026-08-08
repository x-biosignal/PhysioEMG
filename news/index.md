# Changelog

## PhysioEMG 0.2.1

### Bug fixes

- [`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
  now applies the rectify-and-smooth conditioning of the published
  methods before thresholding: the rectified signal (“hodges_bui”) or
  rectified Teager-Kaiser trace (“teager_kaiser”) is smoothed with a
  centered moving average controlled by the new `smooth_ms` argument
  (default 25 ms, in the 10–50 ms range examined by Hodges & Bui 1996;
  `smooth_ms = 0` restores the previous behavior). Previously the
  threshold was applied to the unsmoothed rectified trace, which dips
  below threshold at every carrier zero-crossing, so contiguous
  suprathreshold runs almost never reached `min_duration_ms` and both
  threshold methods failed to detect genuine band-limited bursts (zero
  detections at SNR 20 dB); the statistical “bonato” and “aglr”
  detectors were unaffected.

## PhysioEMG 0.2.0

Initial release as a standalone package in the x-biosignal ecosystem,
split from the PhysioExperiment monorepo. PhysioEMG provides
electromyography (EMG) analysis for `PhysioExperiment` objects, spanning
amplitude, timing, spectral fatigue, muscle synergy, and inter-muscular
network methods.

### New Features

- Amplitude analysis:
  - [`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
    extracts the amplitude envelope via sliding-window RMS (O(1)
    prefix-sum implementation), the Hilbert analytic signal, or
    rectification with a moving-average lowpass filter.
  - [`emgAmplitudeNormalize()`](https://x-biosignal.github.io/PhysioEMG/reference/emgAmplitudeNormalize.md)
    normalizes amplitude to maximum voluntary contraction (MVC) or
    within-trial peak.
- Muscle activation timing:
  - [`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
    detects activation onsets and offsets using the Hodges-Bui
    baseline-SD threshold or the Teager-Kaiser energy operator, with a
    minimum-duration gate to reject spurious bursts.
- Spectral fatigue analysis:
  - [`emgFatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigue.md)
    tracks windowed median and mean frequency plus RMS amplitude to
    characterize myoelectric fatigue.
  - [`emgFatigueIndex()`](https://x-biosignal.github.io/PhysioEMG/reference/emgFatigueIndex.md)
    summarizes fatigue as the final-to-initial median frequency ratio.
  - [`emgSpectralMoments()`](https://x-biosignal.github.io/PhysioEMG/reference/emgSpectralMoments.md)
    computes spectral moments (M0, M1, M2) over sliding windows.
- Muscle synergy decomposition:
  - [`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
    factorizes multi-channel EMG into synergies via NMF, PCA, or ICA,
    reporting weight and activation matrices with variance-accounted-for
    (VAF).
  - [`synergyReconstruct()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyReconstruct.md)
    rebuilds signals from a chosen number of synergies, and
    [`synergyCompare()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyCompare.md)
    best-match pairs and correlates two synergy solutions.
- Inter-muscular network analysis:
  - Static connectivity via
    [`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md)
    (Welch magnitude-squared coherence),
    [`emgPartialCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgPartialCoherenceNetwork.md)
    (partial coherence from the inverse cross-spectral matrix), and
    [`emgWPLINetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgWPLINetwork.md)
    (weighted phase-lag index, with an optional debiased estimator).
  - Directed connectivity via
    [`emgDirectedGCNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDirectedGCNetwork.md)
    using pairwise Granger causality (F-statistic or delta-R-squared
    scoring, with p-values).
  - Time-varying connectivity via
    [`emgDynamicWaveletNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDynamicWaveletNetwork.md)
    using Morlet wavelet coherence in sliding windows with
    cone-of-influence masking.
  - [`emgCoordinationStructure()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoordinationStructure.md)
    summarizes network topology (module detection, weighted modularity,
    efficiency, clustering, participation and within-module z-scores).
  - [`emgInterpretNetworkKG()`](https://x-biosignal.github.io/PhysioEMG/reference/emgInterpretNetworkKG.md)
    ranks and annotates high-weight edges with optional node metadata
    and knowledge-graph links.
- Synthetic data generators for examples and testing:
  [`make_emg()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg.md),
  [`make_emg_contraction()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_contraction.md),
  and
  [`make_emg_fatigue()`](https://x-biosignal.github.io/PhysioEMG/reference/make_emg_fatigue.md).

### Bug Fixes

- Corrected the RMS envelope computation in
  [`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
  for accurate windowed amplitude estimates.
