# Select Muscle Synergy Model Order

Sweeps the number of synergies from 1 to `max_synergies`, fits a
multi-restart NMF at each, and selects the model order from the
variance-accounted-for (VAF) curve. The selected order is the smallest
number of synergies whose global VAF reaches `vaf_threshold`; if the
threshold is never reached, the knee (elbow) of the VAF curve is used
instead.

## Usage

``` r
muscleSynergyOrder(
  x,
  max_synergies,
  vaf_threshold = 0.9,
  per_synergy_vaf = 0.75,
  n_restarts = 10L,
  max_iter = 200L,
  tol = 1e-04,
  seed = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object with multi-channel EMG.

- max_synergies:

  Maximum number of synergies to test (capped at the number of
  channels).

- vaf_threshold:

  Global VAF threshold for selection (default: 0.90).

- per_synergy_vaf:

  Per-channel VAF each muscle should reach at the selected order;
  reported as a quality check (default: 0.75).

- n_restarts:

  Number of NMF restarts per synergy count (default: 10).

- max_iter:

  Maximum NMF iterations (default: 200).

- tol:

  NMF convergence tolerance (default: 1e-4).

- seed:

  Optional integer seed making the whole sweep reproducible.

- assay_name:

  Input assay name (default: first assay).

## Value

A list with:

- vaf_curve:

  A data.frame with columns `n_synergies` and `vaf`.

- selected:

  The selected number of synergies.

- selection_rule:

  "vaf_threshold" or "knee".

- vaf_threshold, per_synergy_vaf:

  The thresholds used.

- per_channel_vaf:

  Per-channel VAF at the selected order.

- per_synergy_vaf_met:

  TRUE if every channel's VAF meets `per_synergy_vaf`.

- fit:

  The
  [`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
  result at the selected order.

## References

Cheung, V.C.K. et al. (2005). "Central and sensory contributions to the
activation and organization of muscle synergies during natural motor
behaviors." Journal of Neuroscience, 25(27), 6419-6434.
doi:10.1523/JNEUROSCI.4904-04.2005

Tresch, M.C., Cheung, V.C.K. & d'Avella, A. (2006). "Matrix
factorization algorithms for the identification of muscle synergies."
Journal of Neurophysiology, 95(4), 2199-2212. doi:10.1152/jn.00222.2005

## See also

[`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
for the underlying decomposition

## Examples

``` r
# three synergies driving disjoint muscle pairs
set.seed(1)
tt <- seq_len(500)
H <- sapply(c(150, 250, 350), function(c0) exp(-((tt - c0) / 40)^2))
W <- rbind(c(1, 1, 0, 0, 0, 0), c(0, 0, 1, 1, 0, 0), c(0, 0, 0, 0, 1, 1))
data <- H %*% W + matrix(abs(rnorm(500 * 6, sd = 0.02)), 500, 6)
pe <- PhysioExperiment(assays = list(env = data), samplingRate = 1000)
ord <- muscleSynergyOrder(pe, max_synergies = 5, seed = 1)
ord$selected
#> [1] 3
ord$vaf_curve
#>   n_synergies       vaf
#> 1           1 0.1909256
#> 2           2 0.6017851
#> 3           3 0.9989987
#> 4           4 0.9991099
#> 5           5 0.9995725
```
