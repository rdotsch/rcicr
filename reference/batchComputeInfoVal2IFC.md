# Computes Informational Values for several classification images

Computes the Informational Value of every classification image in a
list, from one 2IFC stimulus file, as
[`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
would for each in turn.

## Usage

``` r
batchComputeInfoVal2IFC(
  target_cis,
  rdata,
  iter = 10000,
  force_gen_ref_dist = FALSE,
  response_seed = NULL,
  baseimage = NULL,
  reference_stimuli = NULL,
  reference_method = c("gram", "images")
)
```

## Arguments

- target_cis:

  A list of classification images, each as returned by
  [`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md);
  for example the result of
  [`batchGenerateCI2IFC`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md).

- rdata:

  Path to the `.Rdata` file written when the stimuli were generated.
  Every classification image in `target_cis` must come from it.

- iter, force_gen_ref_dist, response_seed, reference_method:

  As in
  [`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md),
  applied to every classification image.

- baseimage:

  As in
  [`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md):
  the one base-image label every classification image in `target_cis`
  was computed for. Score CIs for different base images in separate
  calls.

- reference_stimuli:

  The stimulus numbers each classification image was built from, when
  that is not every saved stimulus. `NULL` (the default) scores every
  image against the full set; a vector scores every image against that
  subset; a list with one element per image (`NULL` for the full set)
  gives each its own. To use what
  [`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
  recorded, pass
  `lapply(target_cis, function(ci) attr(ci, "trial_design")$stimuli)`.

## Value

A numeric vector of Informational Values, one per classification image,
named as `target_cis` is.

## Details

Each distinct reference distribution is found or simulated once and used
for every classification image scored against it. The stimulus file is
read at most three times to find stored references, however many images
there are, and again for each reference that has to be simulated. A loop
over
[`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
reads the file twice per image and, when the reference cannot be stored
(a read-only file, a `response_seed`, or `force_gen_ref_dist = TRUE`),
simulates it again for every image. The values are identical to that
loop's.

Messages about the trial design (see "Matching the reference to the
classification image" in
[`computeInfoVal2IFC`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md))
are collected into one message per kind, naming the classification
images concerned, rather than one per image.

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
#>   |                                                                              |                                                                      |   0%  |                                                                              |============                                                          |  17%  |                                                                              |=======================                                               |  33%  |                                                                              |===================================                                   |  50%  |                                                                              |===============================================                       |  67%  |                                                                              |==========================================================            |  83%  |                                                                              |======================================================================| 100%
rdata_file <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)[1]

# two "participants", three trials each
data <- data.frame(
  participant = rep(c("p1", "p2"), each = 3),
  stimulus = 1:6,
  response = sample(c(1, -1), 6, replace = TRUE)
)
cis <- suppressWarnings(batchGenerateCI2IFC(
  data = data, by = "participant", stimuli = "stimulus", responses = "response",
  baseimage = "face", rdata = rdata_file, save_as_png = FALSE, scaling = "none"
))
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%

# iter is kept tiny here for a fast example, in practice use iter >= 10000.
# Each participant saw three of the six stimuli, so each is scored over those.
suppressWarnings(batchComputeInfoVal2IFC(
  cis, rdata_file, iter = 20,
  reference_stimuli = list(1:3, 4:6)
))
#> Building the reference from the saved noise, please wait...
#> Computing reference distribution, please wait...
#> InfoVal reference computed with reference_method = "gram". References from rcicr 1.5.0 and earlier used "images"; pass reference_method = "images" to reproduce them bit for bit.
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> The reference distribution has been saved to the .Rdata file for reuse.
#> Reference for 1 of 2 classification images (reference over 3 of 6 stimuli; reference median = 2.16646130846797; MAD = 0.50885026846389; iterations = 20)
#> Building the reference from the saved noise, please wait...
#> Computing reference distribution, please wait...
#> InfoVal reference computed with reference_method = "gram". References from rcicr 1.5.0 and earlier used "images"; pass reference_method = "images" to reproduce them bit for bit.
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> The reference distribution has been saved to the .Rdata file for reuse.
#> Reference for 1 of 2 classification images (reference over 3 of 6 stimuli; reference median = 1.90258177977704; MAD = 0.489742072371927; iterations = 20)
#> face_participant_p1 face_participant_p2 
#>       -2.618192e-15       -1.287036e+00 
```
