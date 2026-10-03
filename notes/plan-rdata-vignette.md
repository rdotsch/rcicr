# Plan: the `.Rdata` anatomy becomes a vignette

## What changes

The README's "Anatomy of the `.Rdata` file" (745 of its ~2000 words) moves to a new shipped vignette, `vignette("rdata-file")`. The README keeps its "How it works" paragraph on the file being the only link between the two halves, plus one link to the vignette.

A **vignette**, not a pkgdown article: `vignettes/articles/` is `.Rbuildignore`d, so an article exists only on the website and only for `main`. A vignette installs with each version and reads offline, which is what someone reopening an old analysis has.

## The README table is already wrong

Generating a 64px, 20-trial file and computing an InfoVal on it (`generateStimuli2IFC()`, `generateCI()`, `computeInfoVal2IFC(iter = 1000)`, then `ls()` of the loaded file) gives:

```
after generate: base_face_files base_faces generator_version img_size label n_trials noise_type
                nscales p rng_kind seed sigma stimuli_params stimulus_path use_same_parameters
after infoval:  ... reference_norms reference_norms_fingerprint reference_norms_method
                reference_norms_seed reference_norms_source ...
```

Three fields are written but not documented:

| field | written by |
|---|---|
| `reference_norms_method` | `R/generateReferenceDistribution.R:199` (#358) |
| `method` in each `reference_norms_by_base` and `reference_norms_by_stimuli` entry | `R/reference-base.R:357`, `R/reference-stimuli.R:138` (#358) |
| `mask` in each `reference_norms_by_stimuli` entry | `R/reference-stimuli.R:136` (#380) |

The vignette documents all three.

A masked CI (`generateCI(mask = )`, then `computeInfoVal2IFC()`) writes its reference to `reference_norms_by_stimuli`; that entry's names are `reference_stimuli baseimage mask norms response_seed source fingerprint method`.

## Decision for review: show the fields, or check them

**(a) Show.** The vignette generates a small file and prints `ls()` next to the tables. A reader sees the real file, but nothing fails when a field goes undocumented: the three above went unnoticed through two PRs.

**(b) Check (recommended).** The tables are R data frames rendered with `knitr::kable()`, so the documentation *is* data. A hidden chunk builds the four file shapes the package writes (shared parameters; per-base parameters, `use_same_parameters = FALSE` with two bases; a `reference_stimuli` subset; a masked CI) and `stop()`s when:

- a field in a generated file, top-level or inside a `reference_norms_by_*` entry, has no row; or
- a row names a field that no generated file contains.

The build then fails in `R CMD check`, pkgdown and CI on the PR that adds a field, which is when the reasoning is to hand. On CRAN it can only fail if a release ships with drift that CI already reported.

Everything runs at 64px with `iter = 1000`. Before the PR leaves draft, the complete vignette is timed with `tools::buildVignettes()` and compared with the other three vignettes on the same machine; it must not be the slowest of them.

## Other edits in the same PR

- `DECISIONS.md` → "The `.Rdata` anatomy belongs in `README.md`": rewritten in place for the vignette, with the reason above and (if (b)) why the table checks itself.
- `AGENTS.md:112` and `CONTRIBUTING.md:12` point to the vignette instead of the README section. Both files are near their budgets (99 and 172 words spare); the edits replace words, not add them.
- `?generateStimuli2IFC`: `@return` says "Nothing: everything is saved to files"; it gains a pointer to `vignette("rdata-file")`. Regenerate `man/` with `roxygen2::roxygenise()`.
- `_pkgdown.yml`: the vignette goes under "Get started", after `recipes`.
- `NEWS.md` → Documentation: one bullet.
- README → Documentation lists four vignettes.

## Most likely to fail

`kable()` cells containing Markdown (backticks, bold, links) rendering in both `html_vignette` and the pkgdown theme. If one of them escapes it, the fallback is writing the tables as Markdown and parsing their first column in the check chunk.
