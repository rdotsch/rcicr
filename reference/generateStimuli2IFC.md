# Generates 2IFC stimuli

Generate stimuli for a two-image forced-choice reverse correlation task.

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

  Named list of base image files, e.g. `list(aName = 'baseface.jpg')`.
  JPEG and PNG images are accepted, recognised by a `.png`, `.jpg` or
  `.jpeg` extension. Each name labels that base image's stimulus files
  and is the key
  [`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
  uses to find it in the `.Rdata` file, so every element needs a unique
  name. Each image must be square and exactly `img_size` pixels wide:
  rcicr does not resize base images. All of this is checked before any
  stimuli are generated, and an error names the offending entry.

- n_trials:

  Number of trials. Each trial gets two images per base image: one with
  the noise added (original) and one with it subtracted (inverted).

- img_size:

  Width and height of the square stimulus images, in pixels.

- stimulus_path:

  Directory to save the stimuli and the `.Rdata` file to. Required
  unless both `save_as_png` and `save_rdata` are FALSE; there is no
  default. The directory is created if it does not exist; to just try
  the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- label:

  Label put at the start of each file name.

- use_same_parameters:

  Boolean: all base images share one set of noise parameters (`TRUE`) or
  each gets its own (`FALSE`).

- seed:

  Seed for the random number generator, for reproducibility. It is saved
  in the `.Rdata` file, where the default InfoVal reference replays it.
  With `seed = NULL` there is nothing to replay, so InfoVal references
  for that file need an explicit `response_seed`.

- maximize_baseimage_contrast:

  Boolean: rescale the base image's pixel values to maximize its
  contrast. A base image with no contrast at all, every pixel the same
  value, cannot be rescaled and is rejected with an error. It can still
  be used with `maximize_baseimage_contrast = FALSE`.

- noise_type:

  Noise pattern type: `sinusoid` (default) or `gabor`.

- nscales:

  Number of spatial scales (default: 5). Each additional scale adds a
  higher spatial frequency.

- sigma:

  Sigma of the Gabor patches when `noise_type = 'gabor'` (default: 25).

- ncores:

  Number of CPU cores to use (default: `detectCores() - 1`; 2 under
  `R CMD check`, per CRAN policy).

- return_as_dataframe:

  Boolean: return a data frame with the raw noise of the generated
  stimuli (default: `FALSE`), one row per pixel and one column per
  trial. With the default `use_same_parameters = TRUE` every base image
  shares the same noise, so that is all of it. With
  `use_same_parameters = FALSE` and more than one base image, only the
  first base image's noise is returned, because one column per trial
  cannot hold several. The stimuli are still written for every base
  image, and `save_rdata = TRUE` records every parameter set, so nothing
  is missing from the files.

- save_as_png:

  Boolean: write the stimuli to disk as PNG images (default: `TRUE`).

- save_rdata:

  Boolean: save the `.Rdata` file with the stimulus parameters (default:
  `TRUE`). Computing classification images needs that file, so keep this
  `TRUE`; the argument exists mainly for internal use.

## Value

Nothing: everything is saved to files. With
`return_as_dataframe = TRUE`, the data frame described there.

## Details

Saves the stimuli as PNGs, together with an `.Rdata` file holding the
parameters used to generate each stimulus. Analysing the responses later
requires that file.

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
