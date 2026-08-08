# Muscle Synergy Decomposition

Decomposes multi-channel EMG into muscle synergies using matrix
factorization.

## Usage

``` r
muscleSynergy(
  x,
  n_synergies,
  method = c("nmf", "pca", "ica"),
  n_restarts = 1L,
  max_iter = 200L,
  tol = 1e-04,
  seed = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with multi-channel EMG.

- n_synergies:

  Number of synergies to extract.

- method:

  Decomposition method: "nmf" (non-negative matrix factorization), "pca"
  (principal component analysis), or "ica" (independent component
  analysis).

- n_restarts:

  Number of random NMF restarts; the fit with the lowest reconstruction
  error (highest VAF) is returned (default: 1). Ignored for "pca".

- max_iter:

  Maximum iterations for NMF (default: 200).

- tol:

  Convergence tolerance for NMF (default: 1e-4).

- seed:

  Optional integer seed for reproducible random initialization (NMF
  restarts and ICA). If set, the returned best-of-restarts solution is
  deterministic.

- assay_name:

  Input assay name (default: first assay).

## Value

A list with:

- `W`: Synergy weight matrix (n_synergies x channels)

- `H`: Activation pattern matrix (time x n_synergies)

- `vaf`: Variance accounted for (0-1)

- `method`: Method used

- `n_restarts`: Number of restarts used

- `all_vaf`: VAF of each NMF restart (NA for pca/ica)

- `convergence`: For NMF, a list with the best restart's `iterations`
  and `converged` flag (NULL otherwise)

- `original_data`: Original data matrix for reconstruction

## References

Lee, D.D. & Seung, H.S. (1999). "Learning the parts of objects by
non-negative matrix factorization." Nature, 401(6755), 788-791.
doi:10.1038/44565

Tresch, M.C., Cheung, V.C.K. & d'Avella, A. (2006). "Matrix
factorization algorithms for the identification of muscle synergies."
Journal of Neurophysiology, 95(4), 2199-2212. doi:10.1152/jn.00222.2005

## See also

[`muscleSynergyOrder()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyOrder.md)
for model-order (synergy count) selection,
[`synergyReconstruct()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyReconstruct.md)
for reconstructing data from synergies,
[`synergyCompare()`](https://x-biosignal.github.io/PhysioEMG/reference/synergyCompare.md)
for comparing synergy solutions
