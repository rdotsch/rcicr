# Computes cumulative trial CIs correlations with final/target CI

Correlates the CI built from the first trials with the final or target
CI, adding trials one step at a time.

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

- targetci:

  Optional target CI to correlate the cumulative CIs with, as returned
  by
  [`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md).
  Without it, the final CI of these trials is used.

- step:

  Number of trials added between successive correlations.

## Value

Vector of correlations between each cumulative CI and the final or
target CI.

## Details

Plot the resulting curve to estimate how many trials your task needs.

## Repeated presentations of the same stimulus

This function takes the trials in the order they were presented and does
not average repeated presentations of a stimulus.
[`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
does average them, per unique stimulus, before building its
classification image. Averaging here would discard the presentation
order that a cumulative curve is about.

Without a `targetci`, the final CI is built here from the same trials as
the curve. When the evaluated trials reach the last one, as they always
do at the default `step = 1`, the curve's last point compares that CI
with itself and is exactly 1. That shows self-consistency, not
convergence. A larger `step` can stop short, because trials are taken at
`seq(1, length(responses), step)`: with six responses and `step = 2`,
the last trial evaluated is the fifth, and the curve ends at whatever
that partial CI correlates to (0.97 in one such set, not 1).

This assumes the CI compared against varies at all. Responses that
cancel exactly, with every presentation of a stimulus answered both
ways, average to a CI that is zero everywhere. A correlation with a
constant is undefined, so then **every** point on the curve is `NA`.
Such a curve means the responses carry no net signal, not that the call
failed.

A `targetci` with masked pixels
([`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
stores `NA` in every pixel a `mask` excludes) is correlated over the
unmasked pixels only. If the mask covers *every* pixel, no pairs remain
and the whole curve is `NA`, as above.

If every stimulus was shown equally often, the final CI computed here is
identical to the one `generateCI` returns. If not, the two weight the
data differently (each trial equally here, each unique stimulus equally
there) and they diverge: on an 8-trial set with counts 4/2/1/1 they
correlate at 0.77.

So to see how the CI approaches the one you will report, pass that CI as
`targetci = generateCI(...)` instead of relying on the default, which is
always built without a `mask`.

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
