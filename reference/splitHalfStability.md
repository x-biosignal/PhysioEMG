# Split-half stability of inter-muscle coordination networks

Tests whether an inter-muscle coordination network is real signal or a
split-specific artefact: the recording is split into two independent
halves, the same network is built on each with
[`coordinationNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/coordinationNetwork.md),
and the edge weights (upper triangle) and node strengths (row sums) are
correlated between the halves. A correlation near 1 means each half
rebuilds the same network. Different estimators trade sensitivity for
specificity differently, and a sparse graph's node strengths are more
split-sensitive than its edges, so build several and trust what
reproduces in BOTH edges and topology.

## Usage

``` r
splitHalfStability(
  x,
  methods = c("coherence", "partial_coherence", "wpli", "directed_gc"),
  freq_band = NULL,
  channels = NULL,
  nperseg = 256L,
  noverlap = NULL,
  max_lag = 10L,
  debiased = FALSE,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` holding a multi-channel EMG assay.

- methods:

  Estimators to test (default: all four; see
  [`coordinationNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/coordinationNetwork.md)).

- freq_band, channels, nperseg, noverlap, max_lag, debiased, assay_name:

  Passed to
  [`coordinationNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/coordinationNetwork.md)
  for each half.

## Value

A data frame (class `split_half_stability`), one row per method, with
`method`, `split_half_edge_correlation`,
`split_half_node_strength_correlation`, `first_half_mean_edge`, and
`second_half_mean_edge`.

## See also

[`coordinationNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/coordinationNetwork.md),
[`emgCoordinationStructure()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoordinationStructure.md)

## Examples

``` r
pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
splitHalfStability(pe, methods = "coherence", freq_band = c(20, 150))
#>      method split_half_edge_correlation split_half_node_strength_correlation
#> 1 coherence                   0.4699161                            0.5115727
#>   first_half_mean_edge second_half_mean_edge
#> 1            0.4870532             0.3375243
```
