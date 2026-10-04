# Check a stimulus file against stimulus PNGs

Compares the noise a stimulus `.Rdata` file records with the stimulus
PNGs in a folder, trial by trial. An original and an inverted stimulus
hold the same base image, so their difference is the trial's noise
alone: `ori - inv` equals the noise divided by 0.6 wherever neither
image is clipped at black or white. The base images are therefore not
needed.

## Usage

``` r
checkStimulusPNGs2IFC(rdata, png_dir, label = NULL, seed = NULL)
```

## Arguments

- rdata:

  Path to the `.Rdata` file to check, as written by
  [`generateStimuli2IFC`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md).

- png_dir:

  Directory holding the stimulus PNGs. Required; there is no default.

- label, seed:

  The `label` and `seed` the PNGs were generated with, which their names
  carry. When omitted: the values stored in `rdata`. Pass the archive's
  own when `rdata` is a candidate regenerated with other settings, so
  the PNGs can be found; `seed = NULL` for a set generated with
  `seed = NULL`.

## Value

A data frame with one row per base label and trial: `base`, `trial`,
`share` (the share of compared pixels whose `ori - inv` agrees with the
file's noise divided by 0.6 to within 1/255), `compared` (the number of
pixels neither image clips) and `missing` (`TRUE` when either PNG is
absent, with `share` `NA`). Its `unchecked` attribute lists PNGs in
`png_dir` that carry `label` and `seed` but belong to no row, such as
trials beyond the file's `n_trials`.

No verdict is returned. In the configurations measured in
<https://github.com/rdotsch/rcicr/blob/main/analyses/stimulus-png-residuals.md>,
a wrong `nscales`, seed, noise type, base order or `use_same_parameters`
agreed on 2.1% to 7.4% of pixels per trial, but a Gabor `sigma` of 24
where the PNGs used 25 on 95% to 99.6%: compare candidates against the
same PNGs. Check every base label, since with several bases a wrong
`use_same_parameters` shows only after the first.

Warns about missing and unchecked PNGs, and stops if no PNG named for
`label` and `seed` exists.

## Details

Use it to confirm that a regenerated `.Rdata` file matches an archive of
stimulus PNGs whose original file was lost;
[`vignette("recipes", package = "rcicr")`](https://rdotsch.github.io/rcicr/articles/recipes.md),
"When the `.Rdata` file is lost", shows how.

## See also

[`vignette("recipes", package = "rcicr")`](https://rdotsch.github.io/rcicr/articles/recipes.md),
"When the `.Rdata` file is lost".

## Examples

``` r
base_face <- tempfile(fileext = ".png")
png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
stimulus_path <- tempfile("stimuli")
generateStimuli2IFC(list(face = base_face), n_trials = 3, img_size = 32,
                    stimulus_path = stimulus_path, seed = 1, ncores = 1, nscales = 2)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
rdata <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)
checkStimulusPNGs2IFC(rdata, stimulus_path)
#>   base trial share compared missing
#> 1 face     1     1     1024   FALSE
#> 2 face     2     1     1024   FALSE
#> 3 face     3     1     1024   FALSE
```
