# Generates reference distribution

Generates reference distribution of norms for a particular set of task
parameters.

## Usage

``` r
generateReferenceDistribution2IFC(
  rdata,
  iter = 10000,
  ncores = default_ncores(),
  response_seed = NULL,
  save_rdata = TRUE
)
```

## Arguments

- rdata:

  String pointing to .RData file that was created when stimuli were
  generated. This file contains the contrast parameters of all generated
  stimuli.

- iter:

  Number of iterations for the simulation (i.e., the number of norms
  generated with classification images based on random responding).

- ncores:

  Number of CPU cores to use when re-generating the stimuli (default:
  `detectCores()-1`; 2 under `R CMD check`, per CRAN policy).

- response_seed:

  Optional seed for the simulated random responses. The default (`NULL`)
  draws them from the state left by the stimulus re-generation, which is
  the reproducible behaviour described under Reproducibility. Supply a
  number to obtain an independent draw of the null from the same
  stimuli.

- save_rdata:

  Boolean specifying whether the reference distribution should be
  written back into the `rdata` file (default `TRUE`). Set to `FALSE` to
  compute a distribution without changing what later calls to
  [`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
  will use – worth doing whenever `response_seed` is set, so a one-off
  null does not become the file's permanent reference.

## Value

The reference distribution, invisibly, as a numeric vector of `iter`
norms. Unless `save_rdata = FALSE`, it is also added to the supplied
`rdata` file as `reference_norms` (alongside `reference_norms_seed`,
recording the `response_seed` it was generated with), so a later call to
[`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
using the same file can reuse it instead of re-simulating.

## Details

In order to compute the Informational Value metric. Saves its results in
the supplied rdata file for later reuse.

## Reproducibility

With the default `response_seed = NULL`, the reference distribution is
determined by the stimulus `.Rdata` file alone. It does not depend on
the ambient random number state of the calling session, and it does not
depend on `ncores`. Two researchers who compute InfoVal from the same
stimulus file therefore get the same number, and the same reference
distribution, on different machines and in different sessions.

This is a guarantee, not a coincidence, and it is relied upon: the
function re-generates the stimuli through
[`generateStimuli2IFC`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md),
whose internal [`set.seed()`](https://rdrr.io/r/base/Random.html) call
uses the seed stored in the `.Rdata` file and lands before the random
responses below are drawn.

Pass an explicit `response_seed` to draw a \*different\* null from the
same stimuli – for instance to check how much Monte Carlo error a given
`iter` leaves in your InfoVal. This changes only the simulated
responses; the stimuli themselves, and so the noise basis the null is
built on, are unaffected.

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

# iter is kept tiny here for a fast example; in practice use iter >= 10000.
suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3, ncores = 1))
#> Re-generating stimuli based on rdata file, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%Computing reference distribution, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> 
#> Saving simulated reference distribution to rdata file...
```
