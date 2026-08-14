# Mixture-based muscle synergy analysis

Fits a mixture of factor analysers to per-subject muscle activation data
to identify subgroups of subjects with distinct muscle-synergy
structure, rather than forcing a single global synergy model on a
heterogeneous population. Each cluster's factor loading matrix is that
subgroup's synergy set.

## Usage

``` r
muscleSynergyMixture(
  list_of_data,
  n_clusters,
  n_factors,
  method = c("mfa", "mpca"),
  max_iter = 100,
  tol = 1e-04,
  n_init = 5,
  seed = NULL
)
```

## Arguments

- list_of_data:

  A list of non-empty numeric matrices, one per subject, each
  `time x muscle` (rows = samples/time, columns = muscles; all subjects
  share the same set of muscles/columns).

- n_clusters:

  Number of subgroups (mixture components).

- n_factors:

  Number of synergies (factors) per subgroup.

- method:

  `"mfa"` (mixture of factor analysers, per-muscle diagonal noise;
  default) or `"mpca"` (mixture of probabilistic PCA, isotropic noise).

- max_iter:

  Maximum EM iterations (default 100).

- tol:

  Log-likelihood convergence tolerance (default 1e-4).

- n_init:

  Random/k-means restarts; the best log-likelihood is kept (default 5).

- seed:

  Optional integer seed.

## Value

An object of class `"synergy_mixture"`: a list with `cluster` (subgroup
label per subject), `synergies` (length-`n_clusters` list of
`muscle x factor` loading matrices), `uniqueness` (per-cluster diagonal
noise), `mean` (per-cluster mean), `proportions` (mixing weights),
`responsibilities` (`subject x cluster`), `loglik`, and `bic`.

## References

Ghahramani Z, Hinton GE (1996). The EM algorithm for mixtures of factor
analyzers. Technical Report CRG-TR-96-1, University of Toronto. Matsui
Y. synergyMixR: Mixture-Based Muscle Synergy Analysis.
<https://github.com/matsui-lab/synergyMixR>

## Examples

``` r
set.seed(1)
M <- 6; r <- 2
LA <- matrix(abs(rnorm(M * r)), M, r)          # subgroup A synergies
LB <- matrix(abs(rnorm(M * r)), M, r)          # subgroup B synergies
gen <- function(L, n = 60) t(L %*% matrix(rnorm(r * n), r, n) +
                             matrix(rnorm(M * n, 0, 0.3), M, n))
data <- c(replicate(6, gen(LA), simplify = FALSE),
          replicate(6, gen(LB), simplify = FALSE))
fit <- muscleSynergyMixture(data, n_clusters = 2, n_factors = 2, seed = 1)
fit$cluster
#>  [1] 1 1 1 1 1 1 2 2 2 2 2 2
```
