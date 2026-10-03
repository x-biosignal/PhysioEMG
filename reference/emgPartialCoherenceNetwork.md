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

## Examples

``` r
pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
res <- emgPartialCoherenceNetwork(pe, freq_band = c(20, 150))
res$network
#>           EMG1      EMG2      EMG3      EMG4
#> EMG1 1.0000000 0.2264374 0.2042632 0.2311192
#> EMG2 0.2264374 1.0000000 0.1870111 0.1448351
#> EMG3 0.2042632 0.1870111 1.0000000 0.2229306
#> EMG4 0.2311192 0.1448351 0.2229306 1.0000000
```
