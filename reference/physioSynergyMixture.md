# Mixture-based muscle synergy analysis on PhysioExperiment objects

Thin entry point that runs
[`muscleSynergyMixture()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyMixture.md)
directly on ecosystem data: it reads each subject's `time x muscle`
activation from a PhysioExperiment (or accepts matrices) and fits the
subgroup mixture.

## Usage

``` r
physioSynergyMixture(
  pe_list,
  n_clusters,
  n_factors,
  method = c("mfa", "mpca"),
  assay_name = NULL,
  ...
)
```

## Arguments

- pe_list:

  A list with one element per subject, each a `PhysioExperiment` (its
  assay is read) or a `time x muscle` matrix.

- n_clusters, n_factors:

  Number of subgroups and synergies (factors).

- method:

  `"mfa"` or `"mpca"` (see
  [`muscleSynergyMixture()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyMixture.md)).

- assay_name:

  Assay to read from each PhysioExperiment (default: its default assay).

- ...:

  Further arguments passed to
  [`muscleSynergyMixture()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergyMixture.md).

## Value

A `"synergy_mixture"` object.

## Examples

``` r
set.seed(1)
mk <- function(L) {
  pe <- PhysioCore::PhysioExperiment(
    assays = list(raw = t(L %*% matrix(rnorm(2 * 50), 2, 50) +
                          matrix(rnorm(6 * 50, 0, .3), 6, 50))),
    samplingRate = 100)
  pe
}
LA <- matrix(abs(rnorm(12)), 6, 2); LB <- matrix(abs(rnorm(12)), 6, 2)
pes <- c(replicate(5, mk(LA), simplify = FALSE),
         replicate(5, mk(LB), simplify = FALSE))
fit <- physioSynergyMixture(pes, n_clusters = 2, n_factors = 2, seed = 1)
fit$cluster
#>  [1] 1 1 1 1 1 2 2 2 2 2
```
