# Generates multiple classification images by participant or condition

Generate one classification image per participant or condition, for any
reverse correlation task.

## Usage

``` r
batchGenerateCI(
  data,
  by,
  stimuli,
  responses,
  baseimage,
  rdata,
  save_as_png = TRUE,
  targetpath,
  label = "",
  antiCI = FALSE,
  scaling = "autoscale",
  constant = 0.1
)
```

## Arguments

- data:

  Data frame with one row per trial.

- by:

  Name of the column that splits the data into units, such as
  participants or conditions. One CI is computed per unit.

- stimuli:

  Name of the column holding the stimulus numbers of the presented
  stimuli.

- responses:

  Name of the column holding the responses: 1 where the original
  stimulus was chosen, -1 where the inverted one was.

- baseimage:

  String naming the base image: not its file name, but its key in the
  `base_face_files` list passed to
  [`generateStimuli2IFC`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md).

- rdata:

  Path to the `.Rdata` file written when the stimuli were generated. It
  holds the contrast parameters of every stimulus.

- save_as_png:

  Boolean: also save the CI as a PNG image.

- targetpath:

  Directory to save PNGs to. Required when `save_as_png = TRUE`; there
  is no default. The directory is created if it does not exist; to just
  try the function out, use
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html).

- label:

  Optional string added to the PNG file names to make them easier to
  identify.

- antiCI:

  Boolean: compute the anti-CI, the classification image with its sign
  flipped, instead of the CI.

- scaling:

  Scaling method: `none`, `constant`, `matched`, `independent` or
  `autoscale` (default). `autoscale` computes the CIs unscaled and then
  puts them on one scale with
  [`autoscale`](https://rdotsch.github.io/rcicr/reference/autoscale.md).

- constant:

  Scaling constant for the noise. Used only when `scaling = 'constant'`.

## Value

Named list with one classification image per unit. Each is itself a list
of pixel matrices: the raw noise (`ci`), the scaled noise (`scaled`),
the base image (`base`) and the two combined (`combined`).

## Details

Splits `data` by the `by` column, calls
[`generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
for each part, and returns the CIs. By default each CI is also saved as
a PNG.

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

cis <- suppressWarnings(batchGenerateCI(
  data = data, by = "participant", stimuli = "stimulus", responses = "response",
  baseimage = "face", rdata = rdata_file, save_as_png = FALSE
))
#>   |                                                                              |                                                                      |   0%  |                                                                              |===================================                                   |  50%  |                                                                              |======================================================================| 100%Using scaling factor constant:0.182416456662573
#> 
```
