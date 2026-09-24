# Computes Informational Value

Computes the Informational Value of a single CI from a 2IFC task.

## Usage

``` r
computeInfoVal2IFC(
  target_ci,
  rdata,
  iter = 10000,
  force_gen_ref_dist = FALSE,
  response_seed = NULL,
  baseimage = NULL
)
```

## Arguments

- target_ci:

  A classification image, as the list returned by
  [`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md).

- rdata:

  Path to the `.Rdata` file written when the stimuli were generated. It
  holds the contrast parameters of every stimulus and, once computed,
  the reference distribution (see
  [`generateReferenceDistribution2IFC`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)).

- iter:

  Number of simulated classification images in the reference
  distribution. Used only when the reference distribution has to be
  simulated.

- force_gen_ref_dist:

  Boolean: simulate the reference distribution again even if the `rdata`
  file already holds one (default: `FALSE`).

- response_seed:

  Optional seed for the simulated random responses behind the reference
  distribution. The default, `NULL`, uses the reference distribution
  stored in the `rdata` file, or simulates the reproducible default one
  described under Reproducibility in
  [`generateReferenceDistribution2IFC`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md),
  which needs the stimulus seed saved in the file. For a file without
  one, pass a number. A number forces a fresh reference distribution
  from an independent draw; use it to check how much Monte Carlo error
  `iter` leaves in the Informational Value. That result is *not* written
  back to the `rdata` file, so a one-off check cannot change the number
  every later analysis of the stimulus set reports.

- baseimage:

  The base-image label `target_ci` was computed for, the same one passed
  to
  [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md).
  Required when the base images have different noise parameters; each
  base then gets its own stored reference distribution, and a shared one
  left by an older version is ignored. With a single base image, or base
  images sharing one parameter set, leave it at `NULL`.

## Value

The Informational Value, a z-score.

## Details

The Informational Value is a z-score for the signal in a classification
image: the higher it is, the more signal. A cut-off such as z = 1.96
selects classification images with significant signal at alpha = 0.05.

It is computed against a reference distribution: classification images
simulated from random responses under the same task parameters as the
real data. The Informational Value expresses how unlikely the observed
CI is under the null hypothesis that the responses were random.

Simulating the reference distribution takes a long time. It is simulated
whenever the `rdata` file does not already hold one, and then stored
there for reuse. If the file cannot be written, the simulated reference
is used for this call only: a note names the file (and, for independent
base images, the base) that could not be updated, and the next call
simulates it again. An archive on read-only media can therefore still be
scored, at the cost of simulating each time.

For the method, see Brinkman, L., Goffin, S., van de Schoot, R., van
Haren, N. E. M., Dotsch, R., & Aarts, H. (2019). Quantifying the
informational value of classification images. *Behavior Research
Methods*, *51*, 2059-2073.
[doi:10.3758/s13428-019-01232-2](https://doi.org/10.3758/s13428-019-01232-2)

## Examples

``` r
# a synthetic square grayscale image stands in for a real base face photo
base_face <- tempfile(fileext = ".png")
png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)

stimulus_path <- tempfile("stimuli")
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
#> Building the reference from the saved noise, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%
#> Computing reference distribution, please wait...
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
