# Generates 2IFC classification image

Generate a classification image for a two-image forced-choice reverse
correlation task. This function wraps
[`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
and is kept so that older scripts still run; new code can call
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
directly.

## Usage

``` r
generateCI2IFC(
  stimuli,
  responses,
  baseimage,
  rdata,
  save_as_png = TRUE,
  filename = "",
  targetpath,
  antiCI = FALSE,
  scaling = "independent",
  constant = 0.1
)
```

## Arguments

- stimuli:

  Numeric vector of stimulus numbers, one per response and in the same
  order. Each must be a positive whole number no larger than the number
  of trials saved for the selected base image. Numbers may repeat and
  need not be consecutive. Factors, characters and logicals are
  rejected: if your data hold stimulus labels, check them against the
  generated stimulus filenames before converting them to numbers.

- responses:

  Vector of responses in the same order as `stimuli`: 1 where the
  original stimulus was chosen, -1 where the inverted one was.

- baseimage:

  String naming the base image: not its file name, but its key in the
  `base_face_files` list passed to
  [`generateStimuli2IFC`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md).

- rdata:

  Path to the `.Rdata` file written when the stimuli were generated. It
  holds the contrast parameters of every stimulus.

- save_as_png:

  Boolean: also save the CI as a PNG image.

- filename:

  Optional file name for the PNG image.

- targetpath:

  Directory to save PNGs to. Required when `save_as_png = TRUE`; there
  is no default. The directory is created if it does not exist; to just
  try the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- antiCI:

  Boolean: compute the anti-CI, the classification image with its sign
  flipped, instead of the CI.

- scaling:

  Scaling method: `none`, `constant`, `matched` or `independent`
  (default).

- constant:

  Scaling constant for the noise. Used only when `scaling = 'constant'`.

## Value

List of pixel matrices: the raw classification noise (`ci`), the scaled
noise (`scaled`), the base image (`base`) and the two combined
(`combined`).

## Details

This function returns the classification image (CI) and, by default,
saves it as a PNG. How the CI is scaled for display decides what the
image looks like and whether two images can be compared. The default,
`'independent'`, picks the lowest scaling constant that avoids clipping
this particular CI (see `'constant'` below for the formula). Each CI
therefore gets its own constant, and CIs with different noise ranges
cannot be compared by eye.

`'matched'` scaling matches the range of the CI's pixel intensities to
the range of the base image's. This is nonlinear and depends on the
ranges of both. It also shifts the zero point of the noise: a pixel that
would not have changed the base image before scaling may change it
afterwards, and the other way round. Use it as a quick look at how the
noise affects the base image, not for reporting.

`'constant'` scaling does not depend on the base image or the noise
range, but the constant is yours to choose, with the `constant`
argument. The noise is scaled as
`scaled <- (ci + constant) / (2 * constant)`. Pixel intensities must lie
between 0 and 1; if the scaled noise falls outside that range you get a
warning and should pick a higher constant. The higher the constant, the
fainter the noise in the resulting image. Use the same constant for
every classification image you want to compare.

For several classification images, the lowest constant that works for
all of them is a good choice.
[`autoscale`](https://rdotsch.github.io/rcicr/reference/autoscale.md)
finds it for you.

## Examples

``` r
# a synthetic square grayscale image stands in for a real base face photo
base_face <- tempfile(fileext = ".png")
png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)

stimulus_path <- tempdir()
generateStimuli2IFC(
  base_face_files = list(face = base_face),
  n_trials = 6,
  img_size = 32,
  stimulus_path = stimulus_path,
  seed = 1,
  ncores = 1,
  nscales = 1,
  save_as_png = FALSE
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
rdata_file <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)[1]

responses <- sample(c(1, -1), 6, replace = TRUE)
ci <- generateCI2IFC(
  stimuli = 1:6, responses = responses, baseimage = "face",
  rdata = rdata_file, save_as_png = FALSE
)
```
