# What rcicr stores

Two objects outlive the R session that made them: the `.Rdata` file that
[`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
writes, and the classification image that
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
returns, which you may save yourself. This vignette lists everything in
both. Its tables are checked against freshly generated objects every
time the package is built, so they describe the version you have
installed.

``` r

library(rcicr)
```

## The stimulus `.Rdata` file

The file is the only link between generating stimuli and analysing
responses. Without it, a stimulus set can only be regenerated from the
seed together with every generation setting;
[`vignette("recipes", package = "rcicr")`](https://rdotsch.github.io/rcicr/articles/recipes.md)
shows how, and how to check the result against the stimulus PNGs. Keep
it with your response data and with anything you publish.

[`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
writes one file per call, named
`<label>_seed_<seed>_time_<timestamp>.Rdata`. A small stimulus set shows
what is inside. The base image is synthetic, and everything is shrunk so
this vignette builds in seconds.

``` r

base_face <- tempfile(fileext = ".png")
png::writePNG(outer(seq(0, 1, length.out = 64), seq(1, 0, length.out = 64)), base_face)

stimulus_path <- tempfile("stimuli")
generateStimuli2IFC(list(face = base_face), n_trials = 20, img_size = 64,
                    stimulus_path = stimulus_path, seed = 1, ncores = 1, save_as_png = FALSE)
rdata_file <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)
```

``` r

basename(rdata_file)
#> [1] "rcic_seed_1_time_Oct_04_2026_07_21.Rdata"
stored <- new.env()
load(rdata_file, envir = stored)
ls(stored)
#>  [1] "base_face_files"     "base_faces"          "generator_version"  
#>  [4] "img_size"            "label"               "n_trials"           
#>  [7] "noise_type"          "nscales"             "p"                  
#> [10] "rng_kind"            "seed"                "sigma"              
#> [13] "stimuli_params"      "stimulus_path"       "use_same_parameters"
```

### Written when the stimuli are generated

| Field | What it is |
|:---|:---|
| `p` | The noise basis, described in the next table. This is the expensive part and the reason the file exists. |
| `stimuli_params` | Named list, one entry per base image, each an `n_trials` by `nparams` matrix of contrast weights in \[-1, 1\]. **Row *i* is the noise of stimulus *i***: this is what [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md) looks up and weights by the responses. |
| `base_faces` | Named list of the base images as greyscale matrices, after contrast maximization. The pixels themselves, not paths, so the file is self-contained. |
| `base_face_files` | The paths the base images were read from, for reference. |
| `img_size` | Width and height of the images, in pixels. |
| `n_trials` | Number of stimuli per base image. The InfoVal reference uses it to select trial rows. |
| `nscales` | Number of spatial scales in the basis. Added in 1.1.0. |
| `sigma` | Width of the Gabor patches, used when `noise_type` is `"gabor"`. Added in 1.1.0. |
| `noise_type` | `"sinusoid"` or `"gabor"`. |
| `seed` | The RNG seed. Regenerating the stimuli from it also needs the same generation settings and the same [`RNGkind()`](https://rdrr.io/r/base/Random.html). |
| `rng_kind` | [`RNGkind()`](https://rdrr.io/r/base/Random.html) when the stimuli were drawn. InfoVal references replay the seed’s stream under it, so they do not depend on the kind of the session that computes them. Files without it replay under the session’s kind. Added in the development version. |
| `use_same_parameters` | Whether every base image shared one parameter matrix (`TRUE`) or each got its own. |
| `label` | The label the stimulus files were named with. |
| `stimulus_path` | The directory the files were written to. |
| `generator_version` | The rcicr version that wrote the file, as a `package_version`. Unreliable in older files; see below. |

`p` is a list:

| Field | What it is |
|:---|:---|
| `patches` | An `img_size` by `img_size` by `12 * nscales` array of sinusoid or Gabor layers. |
| `patchIdx` | Which parameter drives each pixel of each layer. |
| `noise_type` | As above. |
| `generator_version` | The rcicr version that built the basis. Unlike the top-level field, this one has always been correct. |

The InfoVal reference reads the saved `p` (or `s` in old files) and
`stimuli_params` directly; it never rebuilds the basis from `img_size`,
`nscales`, `sigma` or `noise_type`, so it works when those are missing.

### Added when an informational value is computed

The first time you compute an informational value,
[`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
and
[`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
**add** fields to the same file. Simulating the reference distribution
is slow, so it is stored and reused.

``` r

ci <- generateCI(stimuli = 1:20, responses = rep(c(1, -1), 10), baseimage = "face",
                 rdata = rdata_file, save_as_png = FALSE)
computeInfoVal2IFC(ci, rdata_file, iter = 1000)
```

``` r

load(rdata_file, envir = stored)
grep("^reference", ls(stored), value = TRUE)
#> [1] "reference_norms"             "reference_norms_fingerprint"
#> [3] "reference_norms_method"      "reference_norms_seed"       
#> [5] "reference_norms_source"
```

Which fields appear depends on whether the base images share one
parameter matrix, and on which stimuli the classification image used:

| Field | What it is |
|:---|:---|
| `reference_norms` | The simulated null distribution: the norms of `iter` classification images built from random responses. Written when the base images share one parameter matrix. |
| `reference_norms_seed` | The `response_seed` those norms were drawn with, or `NULL` for the default stream. Added in 1.2.0. |
| `reference_norms_source` | What `reference_norms` was built from: `"saved_noise"` means the file’s own saved parameters and basis. A default-stream reference without this marker and a matching `reference_norms_fingerprint` is rebuilt and, if the file is writable, saved. A reference drawn with a `response_seed` is kept; if it was built from an incorrect reconstruction, regenerate it yourself before recomputing InfoVal. Added in 1.4.0. |
| `reference_norms_fingerprint` | A full copy of the norms, as `list(norms = ...)`, compared with [`identical()`](https://rdrr.io/r/base/identical.html). An older rcicr can keep the marker while replacing the norms; if the copy no longer matches, a default-stream reference is rebuilt. About 80 KB for 10,000 norms before compression. Added in 1.4.0. |
| `reference_norms_method` | How the norms were computed: `"gram"` (the default) or `"images"`. References from 1.5.0 and earlier, which lack this field, were computed with `"images"`; pass `reference_method = "images"` to reproduce them. Added in the development version. |
| `reference_norms_by_base` | Written instead of the fields above when the base images have *different* parameter matrices, because each base then needs a null built from its own saved noise. A list named by base label, one entry each; the entries are described below. A `reference_norms` already in such a file is left in place and ignored. Added in 1.4.0. |
| `reference_norms_by_stimuli` | References for a classification image that did not use every saved stimulus (`reference_stimuli`) or that is masked. A list with one entry per stimulus set and mask, described below. The default reference is never read or written for these. Added in the development version. |

Each entry of `reference_norms_by_base` and `reference_norms_by_stimuli`
holds:

| Field | What it is |
|:---|:---|
| `norms` | The simulated null distribution, as in `reference_norms`. |
| `response_seed` | As `reference_norms_seed`. |
| `source` | As `reference_norms_source`. |
| `fingerprint` | As `reference_norms_fingerprint`. |
| `method` | As `reference_norms_method`. Added in the development version. |
| `reference_stimuli` | `reference_norms_by_stimuli` only: the stimulus numbers the reference was built over. |
| `baseimage` | `reference_norms_by_stimuli` only: the base label, or `NULL` when the bases share one parameter matrix. |
| `mask` | `reference_norms_by_stimuli` only: the run-length encoding of the mask the classification image was scored under, or `NULL` when it was not masked. Added in the development version. |

### Attributes

Apart from `names` and `dim`, the only attributes in the file are
classes:

| Field | Attribute | What it is |
|:---|:---|:---|
| `generator_version` | `class` | `package_version`, so versions compare with `<` and `>=`. Files from before 1.2.0 hold a character string instead. |
| `p$generator_version` | `class` | `package_version`, in files from every version back to 0.3.3, the oldest in the repository’s history. |

### Before you write code against the file

- **Fields are only ever added.** They are never renamed or given a new
  meaning, so newer rcicr reads older files.
- **`generator_version` is unreliable in older files.** It was hardcoded
  as `'0.4.0'` until 1.2.0, so every file written by 0.4.0 through 1.1.0
  claims to be 0.4.0. `p$generator_version` has always held the real
  value. Compare versions with
  [`numeric_version()`](https://rdrr.io/r/base/numeric_version.html),
  never as text.

## The classification image

[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
returns a list of pixel matrices:

``` r

names(ci)
#> [1] "ci"       "scaled"   "base"     "combined"
```

| Element | What it is |
|:---|:---|
| `ci` | The raw classification noise: the stimuli’s noise, weighted by the responses and averaged. |
| `scaled` | The noise after scaling, as recorded in the `scaling` attribute. |
| `base` | The base image. |
| `combined` | The scaled noise over the base image: what a CI PNG shows. |
| `zmap` | The z-map. Only with `zmap = TRUE`. |

Two attributes record how it was made. They live on the object, not in
the stimulus file, so they are kept only if you save the whole result,
for example with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

``` r

str(attr(ci, "trial_design"))
#> List of 3
#>  $ stimuli       : int [1:20] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ repeated      : logi FALSE
#>  $ n_participants: int 1
str(attr(ci, "scaling"))
#> List of 2
#>  $ method  : chr "independent"
#>  $ constant: num 0.0461
```

`trial_design` lets
[`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
check that its reference matches the classification image:

| Field | What it is |
|:---|:---|
| `stimuli` | The saved stimuli the classification image was built from, sorted and unique. Pass it as `reference_stimuli` to score a CI that did not use every stimulus. |
| `repeated` | Whether any stimulus was presented more than once (per participant, when `participants` was given). |
| `n_participants` | How many participants contributed; 1 when `participants` was not given. |

`scaling` records how `scaled` was made:

| Field | What it is |
|:---|:---|
| `method` | The scaling method applied. An unrecognised one is recorded as the `"none"` used instead; after [`autoscale()`](https://rdotsch.github.io/rcicr/reference/autoscale.md), `"autoscale"`. |
| `constant` | The constant used: `NA` for `"none"` and `"matched"`, the one computed from this CI for `"independent"`, and the shared constant after [`autoscale()`](https://rdotsch.github.io/rcicr/reference/autoscale.md). |
| `individual` | With `participants`: the same `method` and `constant` for the individual CIs, with one constant per participant, named by ID, under `"independent"`. These are the CIs `save_individual_cis` writes, and the record is there whether or not they were written. |
| `combined` | After [`autoscale()`](https://rdotsch.github.io/rcicr/reference/autoscale.md): the `method` and `constant` that still describe `combined`, which [`autoscale()`](https://rdotsch.github.io/rcicr/reference/autoscale.md) leaves as it was. `NULL` when there was no earlier record. |

[`batchGenerateCI()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
and
[`batchGenerateCI2IFC()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md)
return a named list of these classification images, autoscaled by
default.
