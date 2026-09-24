# Generates classification image

Generate a classification image for any reverse correlation task.

## Usage

``` r
generateCI(
  stimuli,
  responses,
  baseimage,
  rdata,
  participants = NA,
  save_individual_cis = FALSE,
  save_as_png = TRUE,
  filename = "",
  targetpath,
  antiCI = FALSE,
  scaling = "independent",
  scaling_constant = 0.1,
  individual_scaling = "independent",
  individual_scaling_constant = 0.1,
  zmap = FALSE,
  zmapmethod = "quick",
  zmapdecoration = TRUE,
  sigma = 3,
  threshold = 3,
  zmaptargetpath,
  n_cores = default_ncores(),
  mask = NA,
  zmappointsize = 12
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

- participants:

  Optional vector with one participant ID per trial, the same length as
  `stimuli` and `responses`. When given, the CI is computed in two
  steps: one CI per participant, then their average. When missing or all
  `NA`, one CI is computed from all trials together. Some, but not all,
  `NA` is an error: give every trial an ID, or remove the trials without
  one.

- save_individual_cis:

  Boolean: when `participants` is given, also save each participant's CI
  as a PNG image.

- save_as_png:

  Boolean: also save the CI as a PNG image.

- filename:

  Optional file name for the PNG image.

- targetpath:

  Directory to save PNGs to. Required when `save_as_png = TRUE` or
  `save_individual_cis = TRUE`; there is no default. The directory is
  created if it does not exist; to just try the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- antiCI:

  Boolean: compute the anti-CI, the classification image with its sign
  flipped, instead of the CI.

- scaling:

  Scaling method: `none`, `constant`, `matched` or `independent`
  (default). When both individual and group CIs are computed, this
  applies to the group CI.

- scaling_constant:

  Scaling constant for the noise, used only when `scaling = 'constant'`.
  When both individual and group CIs are computed, this applies to the
  group CI.

- individual_scaling:

  Scaling method for the individual CIs: `none`, `constant` or
  `independent` (default).

- individual_scaling_constant:

  Scaling constant for the individual CIs, used only when
  `individual_scaling = 'constant'`.

- zmap:

  Boolean: also create a z-map (default: `FALSE`).

- zmapmethod:

  Method for the z-map: `quick` (default) or `t.test`.

- zmapdecoration:

  Boolean: draw the z-map with margins, a caption (sigma, threshold) and
  a scale (default: `TRUE`).

- sigma:

  Amount of smoothing applied when creating the z-map (default: 3).

- threshold:

  Threshold z-score (default: 3). Z-scores below it are not drawn on the
  z-map.

- zmaptargetpath:

  Directory to save z-map PNGs to. Required when `zmap = TRUE`; there is
  no default. The directory is created if it does not exist; to just try
  the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- n_cores:

  Number of CPU cores used to create the z-map (default:
  `detectCores() - 1`; 2 under `R CMD check`, per CRAN policy).

- mask:

  Optional mask for the CI: a 2D matrix (0 = masked, 1 = kept) or the
  path to a greyscale PNG image (black = masked, white = kept). Default:
  `NA`, no mask. Documentation up to and including 1.1.0 described the
  matrix the wrong way round (1 = masked); the code has always masked
  where the matrix is 0.

- zmappointsize:

  Text size of the z-map decoration, in points (default: 12). Passed to
  [`plotZmap`](https://rdotsch.github.io/rcicr/reference/plotZmap.md),
  which makes the z-map image `img_size` pixels wide. The decoration
  needs roughly `12.3 * zmappointsize` pixels on a 72 ppi device and
  `16.4 * zmappointsize` on a 96 ppi one. At the default, a stimulus set
  smaller than about 160 to 200 pixels is too small for it, and
  `generateCI()` stops with an error naming the minimum for your device.
  Lower this value to fit the decoration on a small z-map, or set
  `zmapdecoration = FALSE`.

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
range, but the constant is yours to choose, with the `scaling_constant`
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

## Repeated presentations of the same stimulus

When `participants` is `NA` (the default), repeated presentations of the
same stimulus are averaged before the CI is built, so each unique
stimulus gets equal weight however often it was shown. If every stimulus
was shown equally often, this is the same as weighting each trial
equally. If not, it changes what the CI estimates: a stimulus shown
three times counts the same as one shown once, not three times as much.

So in an unbalanced design, where a participant saw some stimuli more
often than others (an adaptive procedure, a crashed session, or by
design), the CI is the average response per unique stimulus, not per
trial. The two can differ a lot: on an 8-trial set with counts 4/2/1/1
they correlate at 0.77.

[`computeCumulativeCICorrelation`](https://rdotsch.github.io/rcicr/reference/computeCumulativeCICorrelation.md)
does *not* average repeats; it weights each trial equally. With unequal
counts, the final CI it computes itself therefore differs from the one
this function returns. To compare against the CI you will report, pass
this function's output as its `targetci`.

## Examples

``` r
# a synthetic square grayscale image stands in for a real base face photo
base_face <- tempfile(fileext = ".png")
png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)

stimulus_path <- tempfile("stimuli")
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
ci <- generateCI(
  stimuli = 1:6, responses = responses, baseimage = "face",
  rdata = rdata_file, save_as_png = FALSE
)
```
