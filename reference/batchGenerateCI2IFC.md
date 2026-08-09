# Generates multiple 2IFC classification images by participant or condition

Generate classification image for 2 images forced choice reverse
correlation task.

## Usage

``` r
batchGenerateCI2IFC(
  data,
  by,
  stimuli,
  responses,
  baseimage,
  rdata,
  save_as_png = TRUE,
  targetpath,
  antiCI = FALSE,
  scaling = "autoscale",
  constant = 0.1,
  label = ""
)
```

## Arguments

- data:

  Data frame

- by:

  String specifying column name that specifies the smallest unit
  (participant, condition) to subset the data on and calculate CIs for.

- stimuli:

  String specifying column name in data frame that contains the stimulus
  numbers of the presented stimuli.

- responses:

  String specifying column name in data frame that contains the
  responses coded 1 for original stimulus selected and -1 for inverted
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
  `matched`, `independent`, or `autoscale` (default).

- constant:

  Optional number specifying the value used as constant scaling factor
  for the noise (only works for `scaling='constant'`).

- label:

  Optional string to insert in file names of PNGs to make them easier to
  identify.

## Value

List of classification image data structures (which are themselves lists
of pixel matrix of classification noise only, scaled classification
noise only, base image only and combined).

## Details

This function saves the classification images by participant or
condition as PNG to a folder and returns the CIs.

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

# two "participants", three trials each
data <- data.frame(
  participant = rep(c("p1", "p2"), each = 3),
  stimulus = 1:6,
  response = sample(c(1, -1), 6, replace = TRUE)
)

cis <- suppressWarnings(batchGenerateCI2IFC(
  data = data, by = "participant", stimuli = "stimulus", responses = "response",
  baseimage = "face", rdata = rdata_file, save_as_png = FALSE
))
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%Using scaling factor constant:0.182416456662573
#> 
```
