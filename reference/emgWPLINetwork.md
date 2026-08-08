# EMG Weighted Phase-Lag Index (wPLI) Network

Builds a static network from pairwise weighted phase-lag index (wPLI),
which is less sensitive to zero-lag coupling artifacts.

## Usage

``` r
emgWPLINetwork(
  x,
  freq_band = NULL,
  channels = NULL,
  nperseg = 256L,
  noverlap = NULL,
  assay_name = NULL,
  aggregate = c("mean", "max", "median"),
  threshold = NULL,
  debiased = FALSE
)
```

## Arguments

- x:

  A PhysioExperiment object.

- freq_band:

  Optional numeric vector `c(low, high)` in Hz.

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

- debiased:

  Logical; if TRUE, uses debiased wPLI estimator.

## Value

A list with:

- network:

  Numeric matrix (channel x channel) of wPLI values.

- adjacency:

  Logical matrix after thresholding, or NULL.

- wpli:

  3D array (freq x channel x channel).

- frequencies:

  Frequency vector (Hz).

- channel_names:

  Channel labels used in the network.

## See also

[`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md),
[`emgPartialCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgPartialCoherenceNetwork.md)
