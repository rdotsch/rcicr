# Generates 2IFC classification image

Generate classification image for 2 images forced choice reverse
correlation task. This function exists for backwards compatibility. You
can also just use
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md),
which this function wraps.

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

- save_as_png:

  Boolean stating whether to additionally save the CI as PNG image.

- filename:

  Optional string to specify a file name for the PNG image.

- targetpath:

  String specifying the directory to save PNGs to. Required when
  `save_as_png = TRUE`; there is no default path. It is created if it
  does not exist. Use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) if you only want
  to try the function out.

- antiCI:

  Optional boolean specifying whether antiCI instead of CI should be
  computed.

- scaling:

  Optional string specifying scaling method: `none`, `constant`,
  `matched`, or `independent` (default).

- constant:

  Optional number specifying the value used as constant scaling factor
  for the noise (only works for `scaling='constant'`).

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
ci <- generateCI2IFC(
  stimuli = 1:6, responses = responses, baseimage = "face",
  rdata = rdata_file, save_as_png = FALSE
)
```
