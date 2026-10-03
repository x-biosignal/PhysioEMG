# EMG Coordination Structure Summary from Network Topology

Quantifies higher-order muscle coordination structure from a weighted
network matrix, including module structure, efficiency, and node roles.

## Usage

``` r
emgCoordinationStructure(
  network,
  threshold = NULL,
  n_modules = NULL,
  max_modules = 6L,
  directed = FALSE,
  symmetrize = c("mean", "max", "min"),
  normalize = TRUE
)
```

## Arguments

- network:

  Numeric square matrix (channels x channels) or a list containing
  `$network`.

- threshold:

  Optional edge-weight threshold. Values below threshold are set to zero
  before topology estimation.

- n_modules:

  Optional number of modules. If NULL, chooses a value automatically by
  maximizing weighted modularity over candidates.

- max_modules:

  Maximum number of candidate modules for automatic search.

- directed:

  Logical; set TRUE if `network` is directed/asymmetric.

- symmetrize:

  Method to convert directed matrices to undirected form: "mean", "max",
  or "min".

- normalize:

  Logical; if TRUE, rescales weights to `[0, 1]`.

## Value

A list with:

- network:

  Processed undirected weighted network matrix.

- node_metrics:

  Data.frame of node-level topology features.

- modules:

  Named integer vector of module assignments.

- summary:

  Global network topology summary.

## See also

[`emgCoherenceNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoherenceNetwork.md),
[`emgDynamicWaveletNetwork()`](https://x-biosignal.github.io/PhysioEMG/reference/emgDynamicWaveletNetwork.md)

## Examples

``` r
pe <- make_emg(n_time = 2000, n_channels = 4, sr = 1000)
net <- emgCoherenceNetwork(pe, freq_band = c(20, 150))
cs <- emgCoordinationStructure(net)
cs$summary
#> $n_nodes
#> [1] 4
#> 
#> $n_edges
#> [1] 6
#> 
#> $density
#> [1] 1
#> 
#> $mean_edge_weight
#> [1] 0.779085
#> 
#> $global_efficiency
#> [1] 0.779085
#> 
#> $mean_clustering
#> [1] 0.7673964
#> 
#> $modularity
#> [1] -0.09834864
#> 
#> $n_modules
#> [1] 2
#> 
```
