# EMG Coherence Network Analysis

Builds a static muscle coordination network from pairwise
magnitude-squared coherence between EMG channels.

## Usage

``` r
emgCoherenceNetwork(
  x,
  freq_band = NULL,
  channels = NULL,
  nperseg = 256L,
  noverlap = NULL,
  assay_name = NULL,
  aggregate = c("mean", "max", "median"),
  threshold = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object.

- freq_band:

  Optional numeric vector `c(low, high)` in Hz. If NULL, uses all
  frequencies.

- channels:

  Integer vector of channel indices to include. If NULL, uses all.

- nperseg:

  Segment length for Welch estimation (default: 256).

- noverlap:

  Overlap length (default: `floor(nperseg / 2)`).

- assay_name:

  Input assay name. If NULL, uses default assay.

- aggregate:

  Aggregation across frequency bins: "mean", "max", or "median".

- threshold:

  Optional threshold for binary adjacency matrix.

## Value

A list with:

- network:

  Numeric matrix (channel x channel) of coherence strength.

- adjacency:

  Logical matrix after thresholding, or NULL.

- coherence:

  3D array (freq x channel x channel).

- frequencies:

  Frequency vector (Hz) corresponding to `coherence`.

- channel_names:

  Channel labels used in the network.

## See also

[`emgDynamicWaveletNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDynamicWaveletNetwork.md)
for time-varying networks,
[`emgInterpretNetworkKG()`](https://x-biosignal.github.io/PhysioEMG/reference/emgInterpretNetworkKG.md)
for annotation-aware interpretation.
