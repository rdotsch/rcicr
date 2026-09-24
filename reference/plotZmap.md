# Plots a Z-map

Plots a z-map: a matrix of z-scores drawn over a background image.

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
  ...,
  pointsize = 12
)
```

## Arguments

- zmap:

  Matrix of z-scores. Must have the same dimensions as `bgimage`.

- bgimage:

  Greyscale image matrix to draw behind the z-map, normally the base
  image or the final CI. Without it, only the z-map is drawn.

- sigma:

  Sigma of the smoothing applied to the CI to create the z-map, shown in
  the decoration.

- threshold:

  Threshold z-score (default: 3). Z-scores below it are not drawn on the
  z-map.

- mask:

  Optional mask: a binary matrix the size of `zmap` (0 or `FALSE` =
  masked, 1 or `TRUE` = kept), or the path to a black-and-white PNG
  image (black = masked, white = kept). This is the same convention as
  [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)'s
  `mask`, so one mask works for both. Earlier versions of this
  documentation described the matrix the wrong way round.

- decoration:

  Boolean: draw the z-map with margins, a caption (sigma, threshold) and
  a scale (default: `TRUE`).

- targetpath:

  Directory to save the z-map PNG to. Required, since writing that file
  is this function's purpose; there is no default. The directory is
  created if it does not exist; to just try the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- filename:

  Optional file name for the z-map PNG.

- size:

  Width and height of the PNG image, in pixels (default: 512).

- ...:

  Further arguments passed to
  [`graphics::image`](https://rdrr.io/r/graphics/image.html). Used only
  when `decoration = TRUE`.

- pointsize:

  Text size of the decoration, in points (default: 12, the graphics
  device's own default). Margins are measured in lines of text, so this
  also sets how much of the image the decoration takes up. The minimum
  image size is fixed in inches, so in pixels it depends on the device's
  resolution: roughly `12.3 * pointsize` pixels at 72 ppi (Linux, macOS)
  and `16.4 * pointsize` at 96 ppi (Windows), about 160 and 200 pixels
  at the default. Below that, `plotZmap()` stops and names the minimum
  for your device. A lower `pointsize` fits a decorated z-map onto a
  small image, at the cost of a smaller map: the margins shrink but the
  labels still need room. Ignored when `decoration = FALSE`, which has
  no margins and works at any size.

## Value

Nothing; the z-map is written as a PNG.

## Details

Takes a matrix of z-scores, such as the `zmap` that
[`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
returns, and optionally a background image matrix, and writes the z-map
as a PNG. Further arguments are passed to
[`graphics::image`](https://rdrr.io/r/graphics/image.html); for example,
`col` sets a different colour palette. Versions up to and including
1.2.3 passed them to the raster package's plot method instead. `col`
works the same in both, but arguments specific to that method are no
longer understood.

## Reproducibility across platforms

The z-scores are ordinary R arithmetic and do not depend on your
operating system; nor do classification images, scaling or informational
value. The PNG this function writes does: it is drawn through a graphics
device, and devices differ between platforms in colour management and in
whether they write an alpha channel. The same z-map rendered on Linux
and on macOS gives figures that look identical but are not
byte-identical. macOS renders a mid-grey background at roughly 0.573,
where the cairo device gives 0.502.

So to check that an analysis reproduces, compare the numbers, not the
rendered figures. A z-map PNG that differs pixel for pixel on a
colleague's machine is not a different result, and regenerating figures
on another platform is safe.

This applies only to `plotZmap()`, the one function in the package that
opens a graphics device. Every other PNG `rcicr` writes (stimuli,
classification images, autoscaled classification images) comes straight
from the pixel array via
[`png::writePNG()`](https://rdrr.io/pkg/png/man/writePNG.html).

## Examples

``` r
set.seed(1)
zmap <- matrix(rnorm(64, sd = 5), 8, 8)
plotZmap(zmap, sigma = 3, threshold = 3, decoration = FALSE,
         targetpath = tempdir(), size = 200)
```
