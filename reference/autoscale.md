# Determines optimal scaling constant for a list of CIs

Determines optimal scaling constant for a list of CIs

## Usage

``` r
autoscale(cis, save_as_pngs = TRUE, targetpath)
```

## Arguments

- cis:

  List of cis, each of which are a list containing the pixel matrices of
  at least the noise pattern (`$ci`) and if the noise patterns need to
  be written to PNGs, also the base image (`$base`).

- save_as_pngs:

  Boolean, when set to true, the autoscaled noise patterns will be
  combined with their respective base images and saved as PNGs (using
  the key of the list as name).

- targetpath:

  String specifying the directory to save PNGs to. Required when
  `save_as_pngs = TRUE`; there is no default path. It is created if it
  does not exist. Use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) if you only want
  to try the function out.

## Value

The input `cis` list, with the `$scaled` pixel matrix of each element
replaced by its autoscaled version. The scaling constant that was
determined is printed to the console, not returned.

**Look at `$scaled`, not `$combined`.** `$combined` is left exactly as
it was passed in: autoscaling deliberately does not disturb it, so a
combination made before this call survives unchanged. Only `$scaled`
reflects the autoscaled result. If you want the autoscaled noise shown
over the base image, build it yourself with `(ci$scaled + ci$base) / 2`
— that is exactly what `save_as_pngs = TRUE` writes to disk.

This matters most after
[`batchGenerateCI`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
or
[`batchGenerateCI2IFC`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md),
which scale with `'none'` before handing over to this function: their
`$combined` is therefore an overlay of the *unscaled* noise and will
look almost blank, while `$scaled` is the image you want.

## Examples

``` r
cis <- list(
  participant1 = list(ci = matrix(runif(64, -0.2, 0.2), 8, 8), base = matrix(0.5, 8, 8)),
  participant2 = list(ci = matrix(runif(64, -0.3, 0.3), 8, 8), base = matrix(0.5, 8, 8))
)
scaled_cis <- autoscale(cis, save_as_pngs = FALSE)
#> Using scaling factor constant:0.298241481184959
```
