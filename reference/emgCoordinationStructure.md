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
