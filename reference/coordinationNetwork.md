# Build an inter-muscle coordination network with a chosen estimator

A single entry point over the four inter-muscle connectivity estimators
([`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md),
[`emgPartialCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgPartialCoherenceNetwork.md),
[`emgWPLINetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgWPLINetwork.md),
[`emgDirectedGCNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDirectedGCNetwork.md)):
builds the coordination network of a multi-muscle EMG recording with the
chosen `method` and routes the relevant parameters to the underlying
estimator. Use it (with
[`splitHalfStability()`](https://x-biosignal.github.io/PhysioEMG/reference/splitHalfStability.md))
to answer "which muscles co-activate, and is that structure real?".

## Usage

``` r
coordinationNetwork(
  x,
  method = c("coherence", "partial_coherence", "wpli", "directed_gc"),
  freq_band = NULL,
  channels = NULL,
  nperseg = 256L,
  noverlap = NULL,
  max_lag = 10L,
  debiased = FALSE,
  assay_name = NULL,
  ...
)
```

## Arguments

- x:

  A `PhysioExperiment` holding a multi-channel EMG assay.

- method:

  Estimator: `"coherence"` (magnitude-squared coherence),
  `"partial_coherence"` (coherence with shared influence removed),
  `"wpli"` (weighted phase-lag index), or `"directed_gc"` (directed
  Granger causality).

- freq_band:

  Optional `c(low, high)` Hz band for the spectral estimators.

- channels:

  Optional integer channel indices to include (default all).

- nperseg, noverlap:

  Welch segment length / overlap for the spectral estimators (coherence
  / partial coherence / wPLI).

- max_lag:

  Model order for directed Granger causality.

- debiased:

  Use the debiased wPLI estimator.

- assay_name:

  Assay to build the network from (default: the default assay).

- ...:

  Passed through to the underlying estimator.

## Value

The estimator's network object – a list whose `$network` is the
inter-muscle adjacency matrix.

## See also

[`splitHalfStability()`](https://x-biosignal.github.io/PhysioEMG/reference/splitHalfStability.md),
[`emgCoordinationStructure()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoordinationStructure.md),
[`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md)
