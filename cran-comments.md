# CRAN comments

## Submission type

Update of rcicr 1.4.1, which CRAN's auto-check service accepted on 2026-09-20. The maintainer address `rdotsch@gmail.com` is unchanged.

## What changed

`NEWS.md` carries the full list. Its "Reproducibility impact" section describes changes to existing results, who is affected and what to do.

**`batchGenerateCI()` and `batchGenerateCI2IFC()` group a tibble correctly.** Called on a tibble with more than one group, versions 1.1.0 to 1.4.1 returned a single classification image built from alternating trials of different groups. A `data.frame`, and a tibble with one group, are unaffected.

**The default InfoVal reference needs the stimulus seed saved in the file.** For a stimulus file without one, the reference was seeded from the clock and could not be reproduced. Such a file now stops with an error that gives the call storing a reproducible, seeded reference. Files with a seed are unaffected.

**Stimulus `.Rdata` files are protected.** Saving a reference distribution into one keeps a verified backup until the save completes, so an interrupted save no longer leaves it unloadable. `generateStimuli2IFC()` never overwrites an existing one, and reserves its name with a lock beside it before generating anything.

**Clearer errors instead of unusable results.** `generateCI()` stops on partly missing participant IDs, which returned an all-`NA` image. An `img_size` that the noise scales cannot tile now gets an error naming the constraint.

## Test environments

Pending: filled in from the checks in `RELEASING.md` step 2, on the release branch.

## R CMD check results

Pending: as above.

## Downstream dependencies

Pending: to be checked against the CRAN package page before submission.

## Notes

* `parallel`/`doSNOW` respect `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2 when it is
  set, so no example, test or vignette uses more than two cores under check.
* Ten test files call `skip_on_cran()`. All are development guards rather than checks that
  validate an installation: a golden-master regression baseline and the cross-platform z-map
  literals, a pipeline smoke test, a signal-recovery test, four serial-versus-parallel
  agreement checks, a progress-reporting check, and the slower path-handling and
  trial-validation cases. They run in CI. Where a file skips only part of itself,
  the single-core and fast cases run everywhere.
