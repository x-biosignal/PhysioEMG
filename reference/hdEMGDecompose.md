# Decompose HD-sEMG into motor units (PhysioHDEMG bridge)

Convenience bridge that forwards to
[`PhysioHDEMG::hdEMGDecompose()`](https://x-biosignal.r-universe.dev/PhysioHDEMG/reference/hdEMGDecompose.html)
for high-density surface EMG motor-unit decomposition by convolutive
blind source separation. The PhysioHDEMG package must be installed. See
[`PhysioHDEMG::hdEMGDecompose()`](https://x-biosignal.r-universe.dev/PhysioHDEMG/reference/hdEMGDecompose.html)
for the full argument list and return value.

## Usage

``` r
hdEMGDecompose(x, ...)
```

## Arguments

- x:

  HD-sEMG data: an `n_time x n_channels` matrix, a `PhysioExperiment`,
  or an `hdemg_sim` object.

- ...:

  Further arguments passed to
  [`PhysioHDEMG::hdEMGDecompose()`](https://x-biosignal.r-universe.dev/PhysioHDEMG/reference/hdEMGDecompose.html).

## Value

An `hdemg_decomposition` object (see PhysioHDEMG).

## See also

[`muscleSynergy()`](https://x-biosignal.github.io/PhysioEMG/reference/muscleSynergy.md)
for whole-muscle synergy decomposition.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- PhysioHDEMG::make_hdemg_sim(n_units = 2, duration_sec = 3)
dec <- hdEMGDecompose(sim, n_units = 4)
} # }
```
