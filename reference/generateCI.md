# Generates classification image

Generate classification image for any reverse correlation task.

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
  mask = NA
)
```

## Arguments

- stimuli:

  Vector with stimulus numbers (should be numeric) that were presented
  in the order of the response vector. Stimulus numbers must match those
  in file name of the generated stimuli.

- responses:

  Vector specifying the responses in the same order of the stimuli
  vector, coded 1 for original stimulus selected and -1 for inverted
  stimulus selected.

- baseimage:

  String specifying which base image was used. Not the file name, but
  the key used in the list of base images at time of generating the
  stimuli.

- rdata:

  String pointing to .RData file that was created when stimuli were
  generated. This file contains the contrast parameters of all generated
  stimuli.

- participants:

  Optional vector specifying participant IDs. If specified, will compute
  the requested CIs in two steps: step 1, compute CI for each
  participant. Step 2, compute final CI by averaging participant CIs. If
  unspecified, the function defaults to averaging all data in the
  stimuli and responses vector.

- save_individual_cis:

  Optional boolean specifying whether individual CIs should be save as
  PNG images when the `participants` parameter is used.

- save_as_png:

  Optional boolean stating whether to additionally save the CI as PNG
  image.

- filename:

  Optional string to specify a file name for the PNG image.

- targetpath:

  String specifying the directory to save PNGs to. Required when
  `save_as_png = TRUE` or `save_individual_cis = TRUE`; there is no
  default path. It is created if it does not exist. Use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) if you only want
  to try the function out.

- antiCI:

  Optional boolean specifying whether antiCI instead of CI should be
  computed.

- scaling:

  Optional string specifying scaling method: `none`, `constant`,
  `matched`, or `independent` (default). This scaling applies to the
  group-level CIs if both individual-level and group-level CIs are being
  generated.

- scaling_constant:

  Optional number specifying the value used as constant scaling factor
  for the noise (only works for `scaling='constant'`). This scaling
  applies to the group-level CIs if both individual-level and
  group-level CIs are being generated.

- individual_scaling:

  Optional string specifying scaling method for individual CIs: `none`,
  `constant`, `independent` (default).

- individual_scaling_constant:

  Optional number specifying the value used as constant scaling factor
  for the noise of individual CIs (only works for
  `individual_scaling='constant'`).

- zmap:

  Boolean specifying whether a z-map should be created (default: FALSE).

- zmapmethod:

  String specifying the method to create the z-map. Can be: `quick`
  (default), `t.test`.

- zmapdecoration:

  Optional boolean specifying whether the Z-map should be plotted with
  margins, text (sigma, threshold) and a scale (default: TRUE).

- sigma:

  Integer specifying the amount of smoothing to apply when generating
  the z-maps (default: 3).

- threshold:

  Integer specifying the threshold z-score (default: 3). Z-scores below
  the threshold will not be plotted on the z-map.

- zmaptargetpath:

  String specifying the directory to save z-map PNGs to. Required when
  `zmap = TRUE`; there is no default path. It is created if it does not
  exist. Use [`tempdir()`](https://rdrr.io/r/base/tempfile.html) if you
  only want to try the function out.

- n_cores:

  Optional integer specifying the number of CPU cores to use to generate
  the z-map (default: `detectCores()-1`; 2 under `R CMD check`, per CRAN
  policy).

- mask:

  Optional 2D matrix that defines the mask to be applied to the CI (0 =
  masked, 1 = unmasked). May also be a string specifying the path to a
  grayscale PNG image (black = masked, white = unmasked). Default: NA.
  Note the matrix convention was documented the wrong way round (as 1 =
  masked) up to and including 1.1.0; the code has always masked where
  the matrix is 0, matching the PNG form, and that is what is described
  here.

## Value

List of pixel matrix of classification noise only, scaled classification
noise only, base image only and combined.

## Details

This function saves the classification image as PNG to a folder and
returns the CI. Your choice of scaling matters. The default,
`'independent'`, picks the lowest scaling constant that avoids clipping
this particular classification image (see `'constant'` scaling below for
the formula), so it is not comparable across classification images with
different noise ranges.

`'matched'` scaling will match the range of the intensity of the pixels
to the range of the base image pixels. This scaling is nonlinear and
depends on the range of both base image and noise pattern. It is truly
suboptimal, because it shifts the 0 point of the noise (that is, pixels
that would not have changed the base image at all before scaling may
change the base image after scaling and vice versa). It is however the
quick and dirty way to see how the CI noise affects the base image.

For more control, use `'constant'` scaling, where the scaling is
independent of the base image and noise range, but where the choice of
constant is arbitrary (provided by the user with the `constant`
parameter). The noise is then scale as follows:
`scaled <- (ci + constant) / (2*constant)`. Note that pixels can take
intensity values between 0 and 1. If your scaled noise exceeds those
values, a warning will be given. You should pick a higher constant (but
do so consistently for different classification images that you want to
compare). The higher the constant, the less visible the noise will be in
the resulting image.

When creating multiple classification images a good strategy is to find
the lowest constant that works for all classification images. This can
be automatized using the `autoscale` function.

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
ci <- generateCI(
  stimuli = 1:6, responses = responses, baseimage = "face",
  rdata = rdata_file, save_as_png = FALSE
)
```
