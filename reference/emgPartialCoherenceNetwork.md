# EMG Partial Coherence Network Analysis

Builds a static network from pairwise partial coherence, estimated from
the inverse cross-spectral density matrix at each frequency.

## Usage

``` r
emgPartialCoherenceNetwork(
  x,
  freq_band = NULL,
  channels = NULL,
  nperseg = 256L,
  noverlap = NULL,
  assay_name = NULL,
  aggregate = c("mean", "max", "median"),
  threshold = NULL,
  ridge = 1e-06
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

- ridge:

  Ridge regularization added to spectral matrix inversion.

## Value

A list with:

- network:

  Numeric matrix (channel x channel) of partial coherence.

- adjacency:

  Logical matrix after thresholding, or NULL.

- partial_coherence:

  3D array (freq x channel x channel).

- frequencies:

  Frequency vector (Hz).

- channel_names:

  Channel labels used in the network.

## See also

[`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md),
[`emgWPLINetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgWPLINetwork.md)
