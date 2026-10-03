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

## Examples

``` r
pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
res <- emgDirectedGCNetwork(pe, max_lag = 5)
res$network
#>          EMG1      EMG2      EMG3     EMG4
#> EMG1 0.000000 3.3421045 4.9050805 2.508087
#> EMG2 5.342390 0.0000000 0.8378886 2.237689
#> EMG3 2.355595 1.7247195 0.0000000 1.201503
#> EMG4 4.071673 0.9506827 3.2564382 0.000000
```
