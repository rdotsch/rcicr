# Plots a Z-map

Plots a Z-map given a matrix of z-scores that maps onto a specified base
image.

## Usage

``` r
plotZmap(
  zmap,
  bgimage = "",
  sigma,
  threshold = 3,
  mask = NULL,
  decoration = TRUE,
  targetpath,
  filename = "zmap",
  size = 512,
  ...
)
```

## Arguments

- zmap:

  A matrix containing z-scores that map onto a given base image. zmap
  and baseimage must have the same dimensions.

- bgimage:

  A matrix containing the grayscale image to use as a background. This
  should be either the base image or the final CI. If not this argument
  is not given, only the Z-map will be drawn.

- sigma:

  The sigma of the smoothing that was applied to the CI to create the
  Z-map.

- threshold:

  Integer specifying the threshold z-score (default: 3). Z-scores below
  the threshold will not be plotted on the z-map.

- mask:

  Optional. A binary matrix with the same dimensions as zmap: cells that
  are 0 (or FALSE) are masked, cells that are 1 (or TRUE) are kept. Can
  also be the filename (as a string) of a black and white PNG image, in
  which case black (0) is masked and white (1) is kept. This is the same
  convention
  [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
  uses for its own `mask` argument, so a mask can be used with both.
  Note that earlier versions of this documentation stated the opposite
  for the matrix form; the description here matches what the code does.

- decoration:

  Optional boolean specifying whether the Z-map should be plotted with
  margins, text (sigma, threshold) and a scale (default: TRUE).

- targetpath:

  String specifying the directory to save the Z-map PNG to. Required:
  this function exists to write a file, and has no default path. It is
  created if it does not exist. Use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) if you only want
  to try the function out.

- filename:

  Optional string to specify a file name for the Z-map PNG.

- size:

  Integer specifying the width and height of the PNG image (default:
  512).

- ...:

  Additional arguments to be passed to
  [`graphics::image`](https://rdrr.io/r/graphics/image.html). Only
  applied when decoration is TRUE.

## Value

Nothing. It writes a Z-map image.

## Details

This function takes in a matrix of z-scores (as returned by generateCI)
and an Rdata file containing a base image. It returns a Z-map image in
PNG format. Unlisted additional arguments will be passed to
[`graphics::image`](https://rdrr.io/r/graphics/image.html). For example,
a different color palette can be specified using the `col` argument. See
[`graphics::image`](https://rdrr.io/r/graphics/image.html) for details.
Versions up to and including 1.2.3 passed these to the raster package's
plot method instead; `col` works the same way in both, but an argument
specific to that method will no longer be understood.

## Reproducibility across platforms

The z-scores themselves are ordinary R arithmetic and do not depend on
your operating system, and neither do classification images, scaling or
informational value. The PNG this function writes is different: it is
drawn through a graphics device, and graphics devices differ between
platforms both in colour management and in whether they write an alpha
channel. The same z-map rendered on Linux and on macOS gives visibly
identical figures whose files are not byte-identical – macOS renders a
mid-grey background at roughly 0.573 where the cairo device gives 0.502.

So when you are checking that an analysis reproduces, compare the
numbers rather than the rendered figures. A z-map PNG that differs
pixel-for-pixel on a colleague's machine is not a different result, and
regenerating figures on another platform is safe.

This applies only to `plotZmap()`, which is the only function in the
package that opens a graphics device. Every other PNG written by `rcicr`
– stimuli, classification images, autoscaled classification images – is
written straight from the pixel array by
[`png::writePNG()`](https://rdrr.io/pkg/png/man/writePNG.html) and
carries no such dependence.

## Examples

``` r
set.seed(1)
zmap <- matrix(rnorm(64, sd = 5), 8, 8)
plotZmap(zmap, sigma = 3, threshold = 3, decoration = FALSE,
         targetpath = tempdir(), size = 200)
```
