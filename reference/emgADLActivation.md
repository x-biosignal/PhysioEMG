# EMG muscle-activation summary of an ADL task

Summarises how an agonist-antagonist muscle pair drives one ADL task
from an EMG envelope: the co-contraction index (via
[`emgCoContraction()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoContraction.md)),
each muscle's peak and mean activation (as percent of a maximal
reference when `mvc_data` is given, else raw envelope units), and the
fraction of the task each muscle is active. Higher antagonist
co-contraction means less efficient, more guarded movement.

## Usage

``` r
emgADLActivation(
  x,
  agonist,
  antagonist,
  mvc_data = NULL,
  task = c("reaching", "drinking", "feeding", "dressing", "grooming"),
  cc_method = c("falconer_winter", "rudolph", "frost"),
  active_frac = 0.2,
  assay_name = NULL
)
```

## Arguments

- x:

  A `PhysioExperiment` of EMG envelope data (rectified / smoothed).

- agonist, antagonist:

  Channel index or `colData$label` of the pair.

- mvc_data:

  Optional `PhysioExperiment` of maximal-effort (MVC) data with the same
  channels; enables percent-MVC amplitudes.

- task:

  ADL task realised: `"reaching"` (d445), `"drinking"` (d560),
  `"feeding"` (d550), `"dressing"` (d540) or `"grooming"` (d520).

- cc_method:

  Co-contraction index method (see
  [`emgCoContraction()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoContraction.md)).

- active_frac:

  Activation threshold as a fraction of a muscle's own peak (default
  0.2) for the active-time summary.

- assay_name:

  EMG assay to use (default: the first assay).

## Value

an `emg_adl_activation` list: `task`, `icf_code`, `cocontraction_index`,
`unit` (`"%MVC"` or `"envelope"`), per-muscle
`agonist_peak`/`agonist_mean`/`antagonist_peak`/`antagonist_mean` and
`agonist_active_frac`/`antagonist_active_frac`.

## See also

[`emgCoContraction()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoContraction.md),
[`emgAmplitudeNormalize()`](https://x-biosignal.github.io/PhysioEMG/reference/emgAmplitudeNormalize.md)

## Examples

``` r
set.seed(1)
n <- 1000; t <- seq_len(n)
ag <- exp(-((t - 500) / 150)^2) + rnorm(n, 0, 0.02)     # agonist burst
an <- 0.2 * exp(-((t - 500) / 150)^2) + rnorm(n, 0, 0.02) # low antagonist
pe <- PhysioCore::PhysioExperiment(
  assays = S4Vectors::SimpleList(envelope = cbind(BIC = ag, TRI = an)),
  colData = S4Vectors::DataFrame(label = c("BIC", "TRI"), type = "EMG"),
  samplingRate = 1000)
emgADLActivation(pe, "BIC", "TRI", task = "drinking")$cocontraction_index
#> [1] 86.30546
```
