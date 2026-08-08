# Interpret EMG Network with Knowledge-Graph Metadata

Adds metadata/KG context to high-weight network edges. This function
does not require a specific backend and works with user-provided tables
exported from PhysioAnnotationHub/physioKG workflows.

## Usage

``` r
emgInterpretNetworkKG(
  network,
  node_metadata = NULL,
  kg_edges = NULL,
  threshold = NULL,
  top_n = 20L,
  window = NULL
)
```

## Arguments

- network:

  A square matrix (channels x channels) or 3D array (window x channels x
  channels).

- node_metadata:

  Optional data.frame containing at least `channel`. Additional columns
  (e.g. `kg_node`, `muscle_name`, `muscle_group`) are carried into the
  edge table.

- kg_edges:

  Optional data.frame of KG links. Expected columns are `node_a`,
  `node_b`, and optional `relation`.

- threshold:

  Optional edge threshold. Default is 75th percentile of upper-triangle
  weights.

- top_n:

  Maximum number of edges returned after thresholding.

- window:

  Optional window index when `network` is a 3D array. If NULL, uses the
  mean across windows.

## Value

A list with:

- edge_table:

  Ranked edge table with optional metadata/KG annotations.

- threshold:

  Applied threshold value.

- network_matrix:

  Matrix used for interpretation.

- summary:

  List of summary statistics.

- kg_relation_summary:

  Relation counts from matched KG links, or NULL.
