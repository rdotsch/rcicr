# rcicr

[![CRAN status](https://www.r-pkg.org/badges/version/rcicr)](https://CRAN.R-project.org/package=rcicr)
[![R-CMD-check](https://github.com/rdotsch/rcicr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rdotsch/rcicr/actions/workflows/R-CMD-check.yaml)

`rcicr` implements **reverse correlation image classification**, a technique from psychophysics for visualizing internal mental representations (for example, of faces). It generates noise-based stimuli for two-image-forced-choice (2IFC) perceptual tasks, and computes "classification images" from participants' responses that reveal what visual features drove their choices.

## Installation

Install from GitHub:

``` r
install.packages('remotes')
remotes::install_github('rdotsch/rcicr')
```

> **`install.packages('rcicr')` does not currently work.** The package was archived on
> CRAN on 2021-06-08 because email to the maintainer had become undeliverable — an old
> university address that stopped working. The code was never the problem, and the
> maintainer address has since been updated. Returning to CRAN is being worked on; until
> then, GitHub is the only source. The last CRAN release was 0.3.4.1, which is several
> years behind.

## Quick example

A minimal 2IFC workflow: generate stimuli from a base face, then turn collected responses into a classification image.

``` r
library(rcicr)

# 1. Generate stimuli: writes an original + inverted noise-blended PNG per
#    trial to stimulus_path, plus an .Rdata file that later analysis needs.
generateStimuli2IFC(
  base_face_files = list(face = "path/to/base_face.jpg"),
  n_trials        = 770,
  img_size        = 512,
  stimulus_path   = "./stimuli",
  seed            = 1
)

# 2. After running the task and collecting responses (1 = original chosen,
#    -1 = inverted chosen), compute the classification image:
generateCI(
  stimuli   = 1:770,                # stimulus numbers, in presentation order
  responses = my_responses,         # 1 / -1 vector, same order as `stimuli`
  baseimage = "face",               # key used in base_face_files above
  rdata     = "./stimuli/rcic_seed_1_time_....Rdata"
)
```

## Documentation

Two vignettes ship with the package:

``` r
vignette("getting-started", package = "rcicr")  # shortest working example
vignette("reverse-correlation-walkthrough", package = "rcicr")  # the full method
```

The walkthrough covers designing a study, generating stimuli, computing classification images for several participants, choosing a scaling method, and telling signal from noise. Its code runs when the package is built, so it cannot drift out of date.

For example datasets and analysis scripts, see [rcicr_examples](https://github.com/rdotsch/rcicr_examples/). There is also an older [Medium post](https://medium.com/@rondotsch/reverse-correlation-image-classification-using-r-a0701648fb0/) covering similar ground; the vignette above supersedes it and is the version kept current with the code.

## Development

``` r
devtools::load_all()   # load the package for interactive development
devtools::test()       # run the test suite
devtools::check()      # full CRAN-style check
```

## Contributing

Contributions, thoughts, and criticisms are very welcome — please [open an issue](https://github.com/rdotsch/rcicr/issues/).

## License

GPL-2
