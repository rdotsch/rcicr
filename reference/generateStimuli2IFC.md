# Generates 2IFC stimuli

Generate stimuli for 2 images forced choice reverse correlation task.

## Usage

``` r
generateStimuli2IFC(
  base_face_files,
  n_trials = 770,
  img_size = 512,
  stimulus_path,
  label = "rcic",
  use_same_parameters = TRUE,
  seed = 1,
  maximize_baseimage_contrast = TRUE,
  noise_type = "sinusoid",
  nscales = 5,
  sigma = 25,
  ncores = default_ncores(),
  return_as_dataframe = FALSE,
  save_as_png = TRUE,
  save_rdata = TRUE
)
```

## Arguments

- base_face_files:

  Named list of base face file names used as base images for stimuli,
  e.g. `list(aName = 'baseface.jpg')`. Accepts JPEG and PNG images,
  recognised by a `.png`, `.jpg` or `.jpeg` extension. Each name labels
  that base image's stimulus files and indexes the .Rdata file that
  [`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
  reads back, so every element must be named, and named uniquely. Each
  image must be square and exactly `img_size` by `img_size` pixels:
  rcicr does not resize base images. All of this is checked before any
  stimuli are generated, and the message names the offending entry.

- n_trials:

  Number specifying how many trials the task will have (function will
  generate two images for each trial per base image: original and
  inverted/negative noise).

- img_size:

  Number specifying the number of pixels that the stimulus image will
  span horizontally and vertically (will be square, so only one integer
  needed).

- stimulus_path:

  String specifying the directory to save the stimuli and the .Rdata
  file to. Required unless both `save_as_png` and `save_rdata` are
  FALSE; there is no default path. It is created if it does not exist.
  Use [`tempdir()`](https://rdrr.io/r/base/tempfile.html) if you only
  want to try the function out.

- label:

  Label to prepend to each file for your convenience.

- use_same_parameters:

  Boolean specifying whether for each base image the same set of
  parameters is used (TRUE) or a unique set is created for each base
  image (FALSE).

- seed:

  Integer seeding the random number generator (for reproducibility).

- maximize_baseimage_contrast:

  Boolean specifying whether the pixel values of the base image should
  be rescaled to maximize its contrast. A base image with no contrast at
  all — every pixel the same value — has nothing to rescale, and is
  rejected with an error rather than silently turned into an all-`NaN`
  base image. Such an image is still usable with
  `maximize_baseimage_contrast = FALSE`.

- noise_type:

  String specifying noise pattern type (defaults to `sinusoid`; other
  options: `gabor`).

- nscales:

  Integer specifying the number of incremental spatial scales. Defaults
  to 5. Higher numbers will add higher spatial frequency scales.

- sigma:

  Number specifying the sigma of the Gabor patch if noise_type is set to
  `gabor` (defaults to 25).

- ncores:

  Number of CPU cores to use (default: `detectCores()-1`; 2 under
  `R CMD check`, per CRAN policy).

- return_as_dataframe:

  Boolean specifying whether to return a data frame with the raw noise
  of the stimuli that were generated (default: FALSE). Data frame
  columns represent pixel values, data frame rows represent stimuli. The
  frame holds one noise image per trial, so it is meaningful only under
  the default `use_same_parameters = TRUE`, where every base image
  shares a single parameter set and one noise image per trial is all
  there is. With `use_same_parameters = FALSE` and more than one base
  image, only the first base image's noise is returned; the frame's
  shape cannot represent trial x base image. Stimuli are still written
  to disk for every base image either way, and `save_rdata = TRUE`
  records the full parameter set, so nothing is lost from the files
  themselves.

- save_as_png:

  Boolean specifying whether to write the stimuli as images to disk
  (default: TRUE).

- save_rdata:

  Boolean specifying whether .RData file with stimulus parameters will
  be saved (default: TRUE). Note: you always need to save the .RData
  file so that you can retrieve the stimulus parameters to compute
  classification images. This function argument exists primarily for
  internal rcicr use.

## Value

Nothing, everything is saved to files, unless return_as_dataframe is set
to TRUE.

## Details

Will save the stimuli as PNGs to a folder, including .Rdata file needed
for analysis of data after data collection. This .Rdata file contains
the parameters that were used to generate each stimulus.

## Examples

``` r
# a synthetic square grayscale image stands in for a real base face photo
base_face <- tempfile(fileext = ".png")
png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)

generateStimuli2IFC(
  base_face_files = list(face = base_face),
  n_trials = 4,
  img_size = 32,
  stimulus_path = tempdir(),
  seed = 1,
  ncores = 1,
  nscales = 1
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
```
