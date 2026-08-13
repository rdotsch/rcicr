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

## Repeated presentations of the same stimulus

This function walks trials in the order they were presented and does not
aggregate repeated presentations of a stimulus, unlike
[`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md),
which averages the responses to each unique stimulus before building its
classification image. That is deliberate: collapsing repeats would
discard the presentation order a cumulative curve is entirely about.

One consequence is worth knowing. With no `targetci`, the final CI
computed here is built from the same un-aggregated trials as the curve.
Where the evaluated trials reach the last one – always so at the default
`step = 1` – the curve's final point compares that CI with itself and is
exactly 1: self-consistency, not evidence of convergence. A larger
`step` can stop short, because trials are taken at
`seq(1, length(responses), step)`: with six responses and `step = 2` the
last one evaluated is the fifth, and the curve ends at whatever that
partial CI correlates to – 0.97 in one such set, not 1.

Both statements assume the CI being compared against varies at all.
Responses that cancel exactly – every presentation of a stimulus
answered both ways – average to a uniformly zero CI, and a correlation
against a constant is undefined, so **every** point on the curve is `NA`
rather than the last one being 1. Such a curve means the responses carry
no net signal, not that the call failed.

A `targetci` carrying masked pixels –
[`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
stores `NA` in every pixel a `mask` excludes – is handled by correlating
over the unmasked pixels only. If the mask covers *every* pixel, there
are no complete pairs and the curve is all-`NA`, same as the
zero-variance case above.

Where every stimulus was presented the same number of times, that final
CI is identical to the one `generateCI` returns. Where repeat counts
differ, the two weight the data differently – each trial equally here,
each unique stimulus equally there – and they diverge: on an 8-trial set
with counts 4/2/1/1 they correlate at 0.77.

So to see how the CI approaches the one you will actually report, pass
it as `targetci = generateCI(...)` rather than relying on the
self-computed default – built without a `mask`, per the note above.

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
