# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`rcicr` is an R package (CRAN-style) implementing the **reverse correlation image classification** technique from psychophysics: generating noise-based stimuli for 2-image-forced-choice (2IFC) perceptual tasks and computing "classification images" (CIs) from participant response data to visualize internal mental representations (e.g., of faces).

## Common commands

This is a standard R package (built with roxygen2 docs, with a testthat suite under `tests/testthat/` and GitHub Actions CI). Typical workflows, run from an R console with working directory set to the package root:

```r
# Regenerate NAMESPACE and man/*.Rd from roxygen comments in R/*.R
roxygen2::roxygenise()

# Load all package functions into an interactive session without installing
devtools::load_all()

# Run the testthat suite
devtools::test()

# Full package check (as CRAN would run it)
devtools::check()
# equivalent from shell:
# R CMD build . && R CMD check --as-cran rcicr_<version>.tar.gz

# Install local copy
devtools::install()
```

Note: `generateStimuli2IFC()` spawns parallel workers via `parallel::makeCluster()` that each `library(rcicr)` themselves (`.packages='rcicr'` in its `foreach` call) — any test or script that calls it needs the package actually **installed** (`devtools::install()`), not just `load_all()`-loaded, or workers fail with "there is no package called 'rcicr'".

Version/date/dependency metadata lives in `DESCRIPTION`; user-facing changes should be noted in `ChangeLog` following its existing per-version format.

## Testing and CI

- `tests/testthat/` has unit tests for every pure/deterministic function (`deg2rad`, `generateSinusoid`, `generateGabor`, `generateNoisePattern`, `generateNoiseImage`, `generateCINoise`), light/targeted tests for the I/O-heavy ones (`generateCI`, `generateCI2IFC`, `autoscale`, `plotZmap`, `computeInfoVal2IFC`, `batchGenerateCI`, `batchGenerateCI2IFC`, `computeCumulativeCICorrelation`, `generateReferenceDistribution2IFC`, `simulateNoiseIntensities`), and an end-to-end smoke test (`test-smoke-pipeline.R`: `generateStimuli2IFC()` → `generateCI()` at tiny image size).
- `tests/testthat/helper-fixtures.R` provides shared fixtures: `make_square_png()` (synthetic base face — never a real photo), `make_fixture_rdata()` (runs a tiny `generateStimuli2IFC()` and returns the `.Rdata` path), `seed_reference_norms()` (pre-seeds a `reference_norms` vector so `computeInfoVal2IFC()` skips the expensive/interactive reference-distribution path).
- `tests/testthat/test-fixed-bugs.R` holds regression tests for the P0 bugs listed in `BACKLOG.md`, all now fixed. They assert *intended* behaviour — never the buggy output they replaced. If one fails, that is a regression; do **not** "fix" it by asserting the broken result, which locks the bug back in. (The suite briefly shipped tests written that way; it was reversed deliberately.)
- `tests/testthat/test-regression-baseline.R` is a golden master pinning the numeric output of the default pipeline — noise basis, classification image, scaling, infoVal — to the values produced *before* the P0 fixes. It is the evidence that those fixes did not change results under normal usage. **If a change turns it red, that change alters researchers' results** and must be documented in `NEWS.md` under "Reproducibility impact" before merging. It is not a test to casually update.
- `.github/workflows/R-CMD-check.yaml` runs `R CMD check` on `ubuntu-latest` for R `release` and `devel` on every push/PR to `main`. `.github/workflows/test-coverage.yaml` runs `covr::package_coverage()` and uploads to Codecov (needs a `CODECOV_TOKEN` repo secret). `codecov.yml` sets lenient coverage thresholds since coverage is deliberately partial (I/O-heavy functions get lighter tests, not full coverage).
- **Any new top-level file must be added to `.Rbuildignore`** unless it genuinely belongs in the built package. `R CMD check` fails a "non-standard file/directory found at top level" NOTE otherwise (and a separate "hidden files" NOTE for dotfiles). This has already caught `codecov.yml`, `.pre-commit-config.yaml`, and `BACKLOG.md` — if you add a config, doc, or scratch file at the repo root, update `.Rbuildignore` in the same commit. Note the file is a set of *regexes anchored with `^`*, not globs (e.g. `^BACKLOG\.md$`).
  - It currently excludes `*.Rproj`, `.Rproj.user`, `.github`, `.claude`, `CLAUDE.md`, `.session-log.md`, `.pre-commit-config.yaml`, `codecov.yml`, `BACKLOG.md`.
  - `.Rbuildignore` itself was listed in `.gitignore` from 2016 until recently (likely an unintentional RStudio-template leftover), so it never actually shipped to CI or other contributors — don't re-add it there.
- `.pre-commit-config.yaml` drives the pre-commit.ci GitHub App, which runs on every PR. It uses only minimal language-agnostic hooks (trailing whitespace, end-of-file, YAML/merge-conflict checks). R-specific hooks (`styler`/`lintr`/`roxygenize`) were deliberately left out — `styler` would reformat nearly every file in one sweep and destroy `git blame`.

## Git and merge strategy

- **Merge pull requests to `main` with squash merges** (`gh pr merge <n> --squash`, or the
  "Squash and merge" button). Decided 2026-07-27; PRs up to and including #132 were
  ordinary merge commits, so `main`'s history changes shape from #131 onward. One commit
  per PR on `main` keeps the history readable and makes `git revert` of a whole change
  straightforward, which matters here because a single PR is usually one self-contained
  fix plus its test and its `NEWS.md` entry.
- The squash commit message is the place to preserve *why* a change was made — the
  per-commit detail on the branch disappears, so anything that a future reader needs
  (measurements, rejected alternatives, reproducibility impact) belongs in the squash
  message or in `NEWS.md`, not only in the branch commits.
- Note this is repo convention, not enforced by GitHub settings — merge commits and
  rebase merges are still permitted in the repository settings, so it is on whoever
  merges to pick the right one.

## Backlog

`BACKLOG.md` is the prioritized modernization backlog — read it before starting substantial work, and update it when you close something out. It was compiled from the open GitHub issues, the published literature, and a direct source review, and it records:

- The **guiding constraint**: researchers re-run old analysis scripts years later, so never change existing call syntax, argument meanings, or the numeric output of a function silently. Deprecate rather than delete; treat the `.Rdata` contract as append-only.
- Several **verified** bugs (reproduced, not inferred) — including that the package is archived on CRAN, that `mask` is unusable on R >= 4.2, that base images must already equal `img_size` (root cause of issue #124), and that `nscales`/`sigma` are missing from the saved `.Rdata` in a way that silently corrupts InfoVal.
- What is **already correct** and should not be "re-fixed" — notably the infoVal formula, which already matches the published Schmitz et al. erratum.

## Architecture

The package has two pipeline halves that share state only through an `.RData` file written at stimulus-generation time:

### 1. Stimulus generation (`generateStimuli2IFC.R`, `generateCI2IFC.R`, `generateStimuli.R`-family)

- `generateNoisePattern()` builds the base noise basis: a stack of sinusoid (or Gabor) patches at multiple orientations, phases, and spatial scales (`nscales`), returned as `p` (a list with `patches`, `patchIdx`, `noise_type`). This is generated **once** per stimulus set and reused for every trial.
- `generateNoiseImage()` takes a per-trial parameter vector (random contrast weights, one per patch/scale index) and combines it with `p` to produce a single noise image.
- `generateStimuli2IFC()` orchestrates the full generation loop: loads base face image(s) (PNG/JPEG, must be square), generates one random parameter vector per trial (`stimuli_params`), and for each trial writes two PNGs per base face — the original stimulus (noise + base face) and its inverted/negative counterpart — using `parallel`/`doParallel`/`foreach` for the per-trial loop.
- Critically, at the end it saves an `.RData` file containing `base_faces`, `stimuli_params`, `p`, `img_size`, `seed`, `noise_type`, etc. **This file is the only link between stimulus generation and later analysis** — every CI-computation function requires it via the `rdata` argument.

### 2. Analysis / classification image computation (`generateCI.R`, `generateCI2IFC.R`, `batchGenerateCI.R`, `batchGenerateCI2IFC.R`, `computeInfoVal2IFC.R`, `computeCumulativeCICorrelation.R`, `plotZmap.R`)

- `generateCI()` / `generateCI2IFC()` load the `rdata` file, look up the stimulus parameters (`stimuli_params[[baseimage]]`) that correspond to the stimulus numbers a participant actually saw, weight them by that participant's responses (1 = original chosen, -1 = inverted chosen), and sum/average them into a single noise image (the classification image). Scaling of the resulting pixel intensities (`none`, `constant`, `matched`, `independent`) is a key user-facing decision documented at length in `generateCI.R`'s roxygen header — read it before changing scaling logic.
- `batchGenerateCI()` / `batchGenerateCI2IFC()` wrap the single-CI functions to compute one CI per participant/condition (grouped by a `by` column in a data frame), optionally followed by `autoscale()` to rescale a whole batch of CIs consistently so they remain visually comparable.
- `autoscale()` finds the single scaling constant that works across a list of CIs without clipping any of them, then rewrites each `$scaled` field and optionally re-saves PNGs.
- `computeInfoVal2IFC()` computes an "informational value" (a z-score-like signal-strength metric) for a CI by comparing it to a simulated null/reference distribution (`generateReferenceDistribution.R`, `simulateNoiseIntensities.R`); a small lookup table of precomputed reference norms avoids expensive resimulation for common parameter combinations (seed/img_size/n_trials/iter) — see the `ref_lookup` tibble in `computeInfoVal2IFC.R`.
- `plotZmap()` computes/plots a spatial z-map of a CI using `spatstat.explore::blur` for smoothing, to highlight image regions with statistically reliable signal.

### Data flow summary

```
base face image(s) ─┐
                     ├─> generateStimuli2IFC() ─> PNGs on disk + <label>_seed_<seed>_time_<ts>.Rdata
random noise params ─┘                                  │
                                                         ▼
participant responses (stimulus #, 1/-1) ──────> generateCI()/generateCI2IFC() ──> classification image
                                                         │
                                            ┌────────────┼───────────────┐
                                            ▼            ▼               ▼
                                       autoscale()  computeInfoVal2IFC()  plotZmap()
```

### Notes on conventions in this codebase

- Functions generally use `save_as_png=TRUE` / `save_rdata=TRUE`-style side-effecting defaults — most analysis functions write PNGs to disk (`targetpath`/`stimulus_path` args) in addition to returning data structures.
- `pre_0.3.0` / `generator_version` fields exist to keep backward compatibility with `.Rdata` files produced by older versions of the package (index-counter starts at 0 vs 1) — do not remove without understanding this.
- `R/zzz.R` declares `utils::globalVariables()` for names that only exist at runtime after `load()`-ing an `.Rdata` file (e.g. `base_faces`, `stimuli_params`, `p`, `seed`) or inside `foreach` loops (`obs`) — when adding new variables loaded this way, add them here to avoid `R CMD check` NOTEs.
- Parallelism is via base `parallel` + `doParallel`/`foreach`, not newer alternatives (e.g. `future`) — match this pattern for new parallel code, and remember worker processes need `.packages='rcicr'` set on `foreach` calls.
