# Generates reference distribution

Generates the reference distribution of norms for a stimulus set.

## Usage

``` r
generateReferenceDistribution2IFC(
  rdata,
  iter = 10000,
  ncores = default_ncores(),
  response_seed = NULL,
  save_rdata = TRUE,
  baseimage = NULL
)
```

## Arguments

- rdata:

  Path to the `.Rdata` file written when the stimuli were generated. It
  holds the contrast parameters of every stimulus.

- iter:

  Number of simulated classification images, each built from random
  responses; the distribution holds one norm per image.

- ncores:

  Number of CPU cores used to rebuild the saved noise (default:
  `detectCores() - 1`; 2 under `R CMD check`, per CRAN policy).

- response_seed:

  Optional seed for the simulated random responses. The default, `NULL`,
  continues from the state the stimulus generator left behind, as
  described under Reproducibility; it needs the stimulus seed saved in
  the file. A number gives an independent draw of the null from the same
  stimuli.

- save_rdata:

  Boolean: write the reference distribution into the `rdata` file
  (default `TRUE`). With `FALSE`, later calls to
  [`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
  keep using what the file already holds. Set it to `FALSE` whenever you
  set `response_seed`, so a one-off null does not become the file's
  permanent reference. The exception is a file without a stimulus seed:
  there the seeded reference is the one to keep, as the file's
  reproducible reference.

- baseimage:

  Base-image label, the same key passed to
  [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md).
  Required when the base images have different noise parameters. With a
  single base image, or base images sharing one parameter set, leave it
  at `NULL`.

## Value

The reference distribution, invisibly, as a numeric vector of `iter`
norms. Unless `save_rdata = FALSE`, it is also added to the `rdata` file
as `reference_norms`, with `reference_norms_seed` recording the
`response_seed` it was drawn with. A later
[`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
call on the same file then reuses it instead of simulating again. For
independent base images it goes in `reference_norms_by_base` instead, as
described above.

## Details

[`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
scores a classification image against this distribution. By default the
result is saved in the `rdata` file for later reuse.

## Reproducibility

With the default `response_seed = NULL`, the reference distribution
depends only on the stimulus `.Rdata` file and the session's
[`RNGkind`](https://rdrr.io/r/base/Random.html): not on the current
random state, and not on `ncores`. Two researchers computing InfoVal
from the same stimulus file get the same reference distribution, and the
same number, on any machine and in any session, provided both use the
same RNG kind.

This needs the stimulus seed saved in the file. A file without one (made
with `generateStimuli2IFC(seed = NULL)`, or with the field removed) has
no stream to continue, so the default stops with an error; pass a
`response_seed` instead.

The RNG kind is the one gap.
[`set.seed()`](https://rdrr.io/r/base/Random.html) keeps whatever kind
the session already has, and no stimulus file records which kind was in
use. A session with a different
[`RNGkind()`](https://rdrr.io/r/base/Random.html) therefore draws
different simulated responses, and so a different null; the saved noise
basis and stimulus parameters stay the same. To reproduce an earlier
reference, use the RNG kind that built it; see
<https://github.com/rdotsch/rcicr/issues/315>.

The noise is rebuilt from the basis and parameters saved in the stimulus
file, so the base images are never reopened and a moved or archived
experiment can still be scored. The simulated responses continue the
random stream from the saved stimulus seed, after the draws for one
parameter matrix, which reproduces the historical default reference.

Pass a `response_seed` to draw a *different* null from the same stimuli,
for instance to check how much Monte Carlo error a given `iter` leaves
in your InfoVal. This changes only the simulated responses, not the
stimuli or the noise basis the null is built on.

## Independent base images

When the base images have different parameter matrices, `baseimage` says
whose noise to use. The distributions are then stored in
`reference_norms_by_base`, one entry per base label, each holding
`norms` and `response_seed`. A shared `reference_norms` in such a file
is neither used nor overwritten. Files whose base images share one
parameter matrix keep using `reference_norms`.

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
#> Building the reference from the saved noise, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%
#> Computing reference distribution, please wait...
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> 
#> Saving simulated reference distribution to rdata file...
```
