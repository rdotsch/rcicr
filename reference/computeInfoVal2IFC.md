# Computes Informational Value

Computes Informational Value for a single CI in a 2IFC task.

## Usage

``` r
computeInfoVal2IFC(
  target_ci,
  rdata,
  iter = 10000,
  force_gen_ref_dist = FALSE,
  response_seed = NULL
)
```

## Arguments

- target_ci:

  A classification image object (list-type) as returned by generateCI

- rdata:

  String pointing to .RData file that was created when stimuli were
  generated. This file contains the contrast parameters of all generated
  stimuli and possibly its corresponding reference distribution
  generated with generateReferenceDistribution().

- iter:

  Number of iterations for the simulation of the reference distribution
  (only used if reference distribution is not already pre-generated and
  present in rdata file)

- force_gen_ref_dist:

  Boolean specifying whether to override the default behavior to use
  pre-computed values for the reference distribution for specific task
  parameters and instead force to recompute the reference distribution
  (default: FALSE).

- response_seed:

  Optional seed for the simulated random responses used to build the
  reference distribution. The default (`NULL`) uses the reference
  distribution stored in the `rdata` file, or generates the reproducible
  default one described under Reproducibility in
  [`generateReferenceDistribution2IFC`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md).
  Supplying a number forces a fresh reference distribution to be
  simulated from an independent draw, which is how you check how much
  Monte Carlo error `iter` leaves in this Informational Value. The
  result is deliberately *not* written back to the `rdata` file, so a
  one-off check cannot change the number every later analysis of that
  stimulus set reports.

## Value

Informational value (z-score)

## Details

The Informational Value metric can be considered as a z-score that
quantifies the signal present in a classification image. The higher the
Informational Value, the more signal. It is possible to use a cut-off
such as z = 1.96 to select classification images with significant signal
under alpha = 0.05.

Informational Value is computed by simulating random responding under
identical task parameters to an empirical dataset (called the reference
distribution). The metric quantifies how unlikely it is to observe these
data under the null-hypothesis that there is no signal (i.e., that there
is only random responding).

The simulation to compute the reference distribution takes a long time,
and is only run locally when pre-computed values for the reference
distribution matching the stimulus set in the .Rdata file have not been
supplied by the rcicr package.

For more information see Brinkman, Goffin, Aarts, van Haren, & Dotsch
(in prep).

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

# compute (and cache in rdata_file) a reference distribution; iter is kept
# tiny here for a fast example, in practice use iter >= 10000.
suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3, ncores = 1))
#> Re-generating stimuli based on rdata file, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%Computing reference distribution, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> 
#> Saving simulated reference distribution to rdata file...

responses <- sample(c(1, -1), 6, replace = TRUE)
target_ci <- generateCI(
  stimuli = 1:6, responses = responses, baseimage = "face",
  rdata = rdata_file, save_as_png = FALSE
)

computeInfoVal2IFC(target_ci = target_ci, rdata = rdata_file)
#> Using reference distribution found in rdata file.
#> Informational value: z = -1.01238122740691 (ci norm = 0.985406648738523; reference median = 1.50357804391794; MAD = 0.51183425882624; iterations = 3)
#> [1] -1.012381
```
