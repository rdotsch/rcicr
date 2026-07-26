# rcicr

[![CRAN status](https://www.r-pkg.org/badges/version/rcicr)](https://CRAN.R-project.org/package=rcicr)
[![R-CMD-check](https://github.com/rdotsch/rcicr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rdotsch/rcicr/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/rdotsch/rcicr/branch/main/graph/badge.svg)](https://codecov.io/gh/rdotsch/rcicr)

`rcicr` implements **reverse correlation image classification**, a technique from psychophysics for visualizing internal mental representations (for example, of faces). It generates noise-based stimuli for two-image-forced-choice (2IFC) perceptual tasks, and computes "classification images" from participants' responses that reveal what visual features drove their choices.

## Installation

Install the latest stable release from CRAN:

``` r
install.packages('rcicr')
```

Or install the development version from GitHub:

``` r
install.packages('devtools')
devtools::install_github('rdotsch/rcicr')
```

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

See [this Medium post](https://medium.com/@rondotsch/reverse-correlation-image-classification-using-r-a0701648fb0/) for a full walkthrough, and [rcicr_examples](https://github.com/rdotsch/rcicr_examples/) for example datasets and scripts.

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
