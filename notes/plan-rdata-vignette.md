# Plan: what rcicr stores becomes a checked vignette

## What changes

A new shipped vignette, `vignette("stored-data")` ("What rcicr stores"), documents the two objects that outlive a session:

1. **The stimulus `.Rdata` file.** Its anatomy moves here from the README (745 of the README's ~2000 words). The README keeps its "How it works" paragraph on the file being the only link between the two halves, plus a link.
2. **The classification image `generateCI()` returns**, which researchers save themselves (`saveRDS()`). Its list elements and its `trial_design` and `scaling` attributes are documented today only in `?generateCI` and `?autoscale`.

A **vignette**, not a pkgdown article: `vignettes/articles/` is `.Rbuildignore`d, so an article exists only on the website and only for `main`. A vignette installs with each version and reads offline, which is what someone reopening an old analysis has.

## Verified layouts

All at 64px with 20 trials; each layout below was printed by generating it and walking `names()` and `attributes()` recursively, `class` included.

**`.Rdata` file.** After `generateStimuli2IFC()`: `base_face_files base_faces generator_version img_size label n_trials noise_type nscales p rng_kind seed sigma stimuli_params stimulus_path use_same_parameters`; `p` holds `patches patchIdx noise_type generator_version`. After `computeInfoVal2IFC()` the file also holds `reference_norms reference_norms_fingerprint reference_norms_method reference_norms_seed reference_norms_source`. A masked CI's reference goes to `reference_norms_by_stimuli`, whose entries hold `reference_stimuli baseimage mask norms response_seed source fingerprint method`. The only attribute beyond `names` and `dim` is `class`: `generator_version` and `p$generator_version` are `package_version` objects. Files from before 1.2.0 hold a character string there instead (`R/rdata.R:34-36`), so the vignette gives the type and its history.

Three of these fields are written but missing from the README table:

| field | written by |
|---|---|
| `reference_norms_method` | `R/generateReferenceDistribution.R:199` (#358) |
| `method` in each `reference_norms_by_base` and `reference_norms_by_stimuli` entry | `R/reference-base.R:357`, `R/reference-stimuli.R:138` (#358) |
| `mask` in each `reference_norms_by_stimuli` entry | `R/reference-stimuli.R:136` (#380) |

**Classification image.**

| call | list elements | `trial_design` | `scaling` |
|---|---|---|---|
| `generateCI()` | `ci scaled base combined` | `stimuli repeated n_participants` | `method constant` |
| with `participants` | same | same | adds `individual` (`method constant`) |
| with `zmap = TRUE` | adds `zmap` | same | same |
| after `autoscale()` | same | same | `method constant combined`, plus `individual` if present |
| `batchGenerateCI2IFC()` | a named list of such CIs | same | as after `autoscale()` |

## The tables check themselves

The tables are R data frames rendered with `knitr::kable()`, so the documentation is data. A hidden chunk builds every layout above and `stop()`s when:

- a field, list element or attribute in a generated object, at any depth listed above, has no row; attributes other than `names` and `dim` count, `class` included; or
- a row names one that no generated object contains.

The build then fails in `R CMD check`, pkgdown and CI on the PR that adds a field, while its reasoning is to hand. On CRAN it can only fail if a release ships with drift that CI already reported.

Everything runs at 64px with `iter = 1000` and an undecorated z-map (a decorated one needs at least 148px at the default pointsize). Before the PR leaves draft, the complete vignette is timed with `tools::buildVignettes()` and compared with the other three vignettes on the same machine; it must not be the slowest of them.

## Other edits in the same PR

- `?generateCI` and `?autoscale`: `@return` keeps one sentence on what each attribute is for and points to `vignette("stored-data")` for its fields, so the field lists exist once. `?generateStimuli2IFC`'s `@return` ("Nothing: everything is saved to files") points there too. `man/` regenerated with `roxygen2::roxygenise()`.
- `DECISIONS.md` → "The `.Rdata` anatomy belongs in `README.md`": rewritten in place for the vignette, with why it is a vignette and why its tables check themselves.
- `AGENTS.md:112` and `CONTRIBUTING.md:12` point to the vignette. Both files are near their budgets (99 and 172 words spare); the edits replace words, not add them.
- `_pkgdown.yml`: the vignette goes under "Get started", after `recipes`.
- `NEWS.md` → Documentation: one bullet.
- README → Documentation lists four vignettes.

## Most likely to fail

`kable()` cells containing Markdown (backticks, bold, links) rendering in both `html_vignette` and the pkgdown theme. If one of them escapes it, the fallback is writing the tables as Markdown and parsing their first column in the check chunk.
