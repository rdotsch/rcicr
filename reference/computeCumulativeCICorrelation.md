# Computes cumulative trial CIs correlations with final/target CI

Computes cumulative trial CIs correlations with final/target CI.

## Usage

``` r
computeCumulativeCICorrelation(
  stimuli,
  responses,
  baseimage,
  rdata,
  targetci = list(),
  step = 1
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

- targetci:

  List Target CI object generated with rcicr functions to correlate
  cumulative CIs with.

- step:

  Step size in sequence of trials to compute correlations with.

## Value

Vector containing correlation between cumulative CI and final/target CI.

## Details

Use for instance for plotting curves of trial-final/target CI
correlations to estimate how many trials are necessary in your task

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
correlations <- suppressWarnings(computeCumulativeCICorrelation(
  stimuli = 1:6, responses = responses, baseimage = "face", rdata = rdata_file
))
#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%
```
