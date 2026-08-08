# Plot Agonist-Antagonist Co-Activation

Overlays the agonist and antagonist envelopes and shades their common
(co-contraction) area, annotated with the co-contraction index.

## Usage

``` r
plotCoactivation(
  x,
  agonist,
  antagonist,
  method = c("falconer_winter", "rudolph", "frost"),
  normalize = TRUE,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object of EMG envelopes.

- agonist, antagonist:

  Channel index or label of the two muscles.

- method:

  Co-contraction index for the annotation (see
  [`emgCoContraction()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoContraction.md)).

- normalize:

  Peak-normalize each channel before plotting (default: TRUE).

- assay_name:

  Input assay name (default: first assay).

## Value

A ggplot object.

## References

Falconer, K. & Winter, D.A. (1985). "Quantitative assessment of
co-contraction at the ankle joint in walking." Electromyography and
Clinical Neurophysiology, 25(2-3), 135-149.

## See also

[`emgCoContraction()`](https://x-biosignal.github.io/PhysioEMG/reference/emgCoContraction.md)

## Examples

``` r
t <- seq_len(1000)
a <- exp(-((t - 400) / 120)^2); b <- exp(-((t - 600) / 120)^2)
pe <- PhysioExperiment(assays = list(env = cbind(a, b)),
  colData = S4Vectors::DataFrame(label = c("TA", "GAS")), samplingRate = 1000)
plotCoactivation(pe, "TA", "GAS")
```
