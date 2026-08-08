# Directed EMG Network via Pairwise Granger Causality

Estimates a directed muscle coordination network using pairwise Granger
causality in the time domain.

## Usage

``` r
emgDirectedGCNetwork(
  x,
  channels = NULL,
  assay_name = NULL,
  max_lag = 10L,
  score = c("f_stat", "delta_r2"),
  threshold = NULL,
  p_value_cutoff = NULL,
  standardize = TRUE
)
```

## Arguments

- x:

  A PhysioExperiment object.

- channels:

  Integer vector of channel indices to include. If NULL, uses all.

- assay_name:

  Input assay name. If NULL, uses default assay.

- max_lag:

  Lag order (in samples) for autoregressive modeling.

- score:

  Directed edge metric: "f_stat" or "delta_r2".

- threshold:

  Optional threshold for adjacency based on selected score.

- p_value_cutoff:

  Optional p-value threshold for adjacency.

- standardize:

  Logical; if TRUE, z-score each channel before GC.

## Value

A list with:

- network:

  Directed numeric matrix (source x target).

- p_values:

  Directed matrix of GC p-values.

- adjacency:

  Logical directed matrix, or NULL.

- channel_names:

  Channel labels used in the network.

- lag:

  Lag order used for modeling.
