# Determines optimal scaling constant for a list of CIs

Determines optimal scaling constant for a list of CIs

## Usage

``` r
autoscale(cis, save_as_pngs = TRUE, targetpath)
```

## Arguments

- cis:

  List of classification images, each a list of pixel matrices holding
  at least the noise (`$ci`), plus the base image (`$base`) if PNGs are
  to be written.

- save_as_pngs:

  Boolean: combine each autoscaled noise pattern with its base image and
  save it as a PNG, named after its key in `cis`.

- targetpath:

  Directory to save PNGs to. Required when `save_as_pngs = TRUE`; there
  is no default. The directory is created if it does not exist; to just
  try the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

## Value

The input `cis` list, with each element's `$scaled` matrix replaced by
its autoscaled version. The scaling constant is printed to the console,
not returned.

**Look at `$scaled`, not `$combined`.** `$combined` is returned exactly
as it was passed in, on purpose, so that existing scripts that plot it
keep producing the same image. Only `$scaled` holds the autoscaled
result. To show the autoscaled noise over the base image, compute
`(ci$scaled + ci$base) / 2`: that is exactly what `save_as_pngs = TRUE`
writes to disk.

This catches people out after
[`batchGenerateCI`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
or
[`batchGenerateCI2IFC`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md).
Both scale with `'none'` before calling this function, so their
`$combined` overlays the *unscaled* noise and looks almost blank, while
`$scaled` is the image you want.

## Examples

``` r
cis <- list(
  participant1 = list(ci = matrix(runif(64, -0.2, 0.2), 8, 8), base = matrix(0.5, 8, 8)),
  participant2 = list(ci = matrix(runif(64, -0.3, 0.3), 8, 8), base = matrix(0.5, 8, 8))
)
scaled_cis <- autoscale(cis, save_as_pngs = FALSE)
#> Using scaling factor constant:0.298241481184959
```
