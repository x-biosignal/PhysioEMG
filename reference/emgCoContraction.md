# EMG Co-Contraction Index

Quantifies simultaneous activation of an agonist-antagonist muscle pair
from their amplitude envelopes, using one of three established indices.
Input channels are treated as (non-negative) amplitude envelopes; with
`normalize = TRUE` each is first scaled to its own peak so the two
muscles are amplitude-comparable.

## Usage

``` r
emgCoContraction(
  x,
  agonist,
  antagonist,
  method = c("falconer_winter", "rudolph", "frost"),
  normalize = TRUE,
  window_sec = NULL,
  assay_name = NULL
)
```

## Arguments

- x:

  A PhysioExperiment object of EMG amplitude envelopes.

- agonist, antagonist:

  Channel index (integer) or channel label (character) of the two
  muscles. The indices are interchangeable for all three methods (the
  indices are symmetric in the pair).

- method:

  Co-contraction index: "falconer_winter" (default), "rudolph", or
  "frost".

- normalize:

  If TRUE (default), peak-normalize each channel to \[0, 1\] before
  computing the index.

- window_sec:

  Optional non-overlapping window length in seconds for a time-resolved
  index. If NULL (default) only the whole-signal summary is returned.

- assay_name:

  Input assay name (default: first assay).

## Value

A list with:

- method:

  The method used.

- summary:

  The whole-signal co-contraction index (a scalar).

- timeseries:

  If `window_sec` is given, a data.frame with columns `window`,
  `time_sec` and `cci` (one row per window); otherwise `NULL`.

## Details

- `"falconer_winter"`:

  Falconer & Winter (1985) percent co-contraction: \\2\sum\min(a,b) /
  \sum(a+b) \times 100\\. 100 when the two envelopes are identical, 0
  when they never overlap.

- `"rudolph"`:

  Rudolph et al. (2000) co-contraction index: the mean over time of
  \\(\ell/h)(\ell+h)\\ with \\\ell=\min(a,b)\\, \\h=\max(a,b)\\.
  Symmetric in the two muscles and non-negative.

- `"frost"`:

  Frost et al. (1997) variant: the time-average of the instantaneous
  co-activation ratio \\2\min(a,b)/(a+b)\\, as a percentage.

## References

Falconer, K. & Winter, D.A. (1985). "Quantitative assessment of
co-contraction at the ankle joint in walking." Electromyography and
Clinical Neurophysiology, 25(2-3), 135-149.

Rudolph, K.S., Axe, M.J. & Snyder-Mackler, L. (2000). "Dynamic stability
after ACL injury: who can hop?" Knee Surgery, Sports Traumatology,
Arthroscopy, 8(5), 262-269. doi:10.1007/s001670000130

Frost, G., Dowling, J., Dyson, K. & Bar-Or, O. (1997). "Cocontraction in
three age groups of children during treadmill locomotion." Journal of
Electromyography and Kinesiology, 7(3), 179-186.
doi:10.1016/S1050-6411(97)84626-3

## See also

[`emgEnvelope()`](https://x-biosignal.github.io/PhysioEMG/reference/emgEnvelope.md)
for computing amplitude envelopes,
[`emgAmplitudeNormalize()`](https://x-biosignal.github.io/PhysioEMG/reference/emgAmplitudeNormalize.md)
for reference normalization,
[`emgOnsetDetect()`](https://x-biosignal.github.io/PhysioEMG/reference/emgOnsetDetect.md)
for activation timing

## Examples

``` r
set.seed(1)
t <- seq_len(1000)
a <- exp(-((t - 400) / 120)^2)            # agonist envelope
b <- exp(-((t - 600) / 120)^2)            # antagonist envelope
pe <- PhysioExperiment(
  assays = list(env = cbind(a, b)),
  colData = S4Vectors::DataFrame(label = c("TA", "GAS"),
                                 type = c("EMG", "EMG")),
  samplingRate = 1000)
emgCoContraction(pe, agonist = "TA", antagonist = "GAS")$summary
#> [1] 23.85986
```
