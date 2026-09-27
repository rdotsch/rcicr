# CRAN comments

## Submission type

Update of rcicr 1.4.1, which CRAN's auto-check service accepted on 2026-09-20. The maintainer address `rdotsch@gmail.com` is unchanged.

This update follows 1.4.1 closely because it fixes incorrect results: `batchGenerateCI()` and `batchGenerateCI2IFC()` returned a wrong classification image for grouped tibbles in versions 1.1.0 to 1.4.1.

## What changed

`NEWS.md` carries the full list. Its "Reproducibility impact" section describes changes to existing results, who is affected and what to do.

**`batchGenerateCI()` and `batchGenerateCI2IFC()` group a tibble correctly.** Called on a tibble with more than one group, versions 1.1.0 to 1.4.1 returned a single classification image built from alternating trials of different groups. A `data.frame`, and a tibble with one group, are unaffected.

**The default InfoVal reference needs the stimulus seed saved in the file.** For a stimulus file without one, the reference was seeded from the clock and could not be reproduced. Such a file now stops with an error that gives the call storing a reproducible, seeded reference. Files with a seed are unaffected.

**Stimulus `.Rdata` files are protected.** Saving a reference distribution into one keeps a verified backup until the save completes, so an interrupted save no longer leaves it unloadable. `generateStimuli2IFC()` never overwrites an existing one, and reserves its name with a lock beside it before generating anything.

**Clearer errors instead of unusable results.** `generateCI()` stops on partly missing participant IDs, which returned an all-`NA` image. An `img_size` that the noise scales cannot tile now gets an error naming the constraint.

## Test environments

The package checks below used commit `662cce416aa9dc96d57301ed67331e2acc6b4edc`, except the local check, which used `2ef6c7683df0321310de12429e04205ce45b72fa`; the only package change between them is the date in the `NEWS.md` release heading.

* Local: R 4.3.3, Ubuntu 24.04.3 LTS, x86_64.
* GitHub Actions (run 36100586167): R 4.6.1 and R-devel on Ubuntu 24.04.5 LTS, R 4.6.1 on Windows and on macOS.
* win-builder, Windows Server 2022 x64 (build 20348):
  * R 4.6.1 (2026-06-24 ucrt).
  * R-devel (2026-09-21 r90579 ucrt).
* R-hub (run 36103799444): R-devel on Linux, Windows and macOS.

## R CMD check results

* win-builder: 0 errors, 0 warnings, 1 NOTE on each environment:
  `Days since last update: 5`. CRAN incoming feasibility and both manuals
  passed. [R 4.6.1](https://win-builder.r-project.org/8h08Tr92o1Jx/),
  [R-devel](https://win-builder.r-project.org/PtX771ROH14R/).
* GitHub Actions: all four jobs passed; `--no-manual`.
* R-hub: all three jobs passed; R-hub uses `--no-manual` and skips the incoming
  feasibility checks.
* Local `R CMD check --as-cran`: 0 errors, 0 warnings, 2 NOTEs. This
  environment could not verify the current time and retained the installation
  lock directory `00LOCK-rcicr`; the package checks, tests, examples,
  vignettes, and PDF and HTML manuals completed successfully.

The full reproducibility gate (run 36100586164) also passed. Against the
published v1.0.1 baseline, 211 checks were identical within tolerance, 22
documented expected deviations fired, and there were 0 unexpected deviations.
Against v1.4.1, 243 checks were identical within tolerance, with 0 expected and
0 unexpected deviations.

## Downstream dependencies

There are no reverse dependencies in the current CRAN package index, checked
with `tools::package_dependencies(..., reverse = TRUE)` against
`https://cloud.r-project.org`.

## Notes

* `parallel`/`doSNOW` respect `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2 when it is
  set, so no example, test or vignette uses more than two cores under check.
* Ten test files call `skip_on_cran()`. All are development guards rather than checks that
  validate an installation: a golden-master regression baseline and the cross-platform z-map
  literals, a pipeline smoke test, a signal-recovery test, four serial-versus-parallel
  agreement checks, a progress-reporting check, and the slower path-handling and
  trial-validation cases. They run in CI. Where a file skips only part of itself,
  the single-core and fast cases run everywhere.
