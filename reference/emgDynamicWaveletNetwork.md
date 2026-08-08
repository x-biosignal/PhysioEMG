# Dynamic Wavelet Coherence Network for EMG

Builds a time-varying coordination network by computing wavelet
coherence between channel pairs and aggregating coherence within sliding
windows.

## Usage

``` r
emgDynamicWaveletNetwork(
  x,
  frequencies = seq(5, 120, by = 5),
  freq_band = NULL,
  channels = NULL,
  window_sec = 0.5,
  step_sec = 0.1,
  n_cycles = 7,
  smoothing_cycles = 3,
  assay_name = NULL,
  aggregate = c("mean", "max", "median"),
  threshold = NULL,
  respect_coi = TRUE
)
```

## Arguments

- x:

  A PhysioExperiment object.

- frequencies:

  Numeric vector of wavelet center frequencies (Hz).

- freq_band:

  Optional numeric vector `c(low, high)` for aggregation.

- channels:

  Integer vector of channel indices to include. If NULL, uses all.

- window_sec:

  Sliding window length in seconds.

- step_sec:

  Sliding window step in seconds.

- n_cycles:

  Number of Morlet cycles (default: 7).

- smoothing_cycles:

  Smoothing width in cycles (default: 3).

- assay_name:

  Input assay name. If NULL, uses default assay.

- aggregate:

  Aggregation across time-frequency bins: "mean", "max", or "median".

- threshold:

  Optional threshold for binary adjacency network per window.

- respect_coi:

  Logical; if TRUE, masks frequencies below COI before aggregation.

## Value

A list with:

- network:

  3D array (window x channel x channel).

- adjacency:

  Logical 3D array thresholded from `network`, or NULL.

- window_times:

  Window center times in seconds.

- static_summary:

  Mean network across windows.

- frequencies:

  Frequency vector used for wavelet transform.

- coi:

  Cone-of-influence frequency at each time sample.

## See also

[`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md)
for static spectral network,
[`emgInterpretNetworkKG()`](https://x-biosignal.github.io/PhysioEMG/reference/emgInterpretNetworkKG.md)
for annotation-aware interpretation.
