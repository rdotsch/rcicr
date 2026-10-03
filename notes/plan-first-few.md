# Plan: one helper for the "first five, ..." lists (#389)

Part of #384. Six error messages build `paste(utils::head(x, 5), collapse = ', ')` plus `', ...'` when longer: `R/ci-inputs.R` (twice), `R/batch-ci.R` (twice), `R/batchComputeInfoVal2IFC.R`, `R/reference-stimuli.R`.

## The change

`firstFew(x, n = 5)` in `R/ci-inputs.R`, returning the same string; each site calls it with the vector it passes today (`sort(repeated)` stays sorted at its call site).

## Proving nothing changed

- Before the change, a test per site pins the complete message for a 3-element and a 7-element case (the `...` branch), unless one exists. These are written first and must pass on `main`.
- **Full** release gate against `main` (`Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`): `0 expected deviations, 0 unexpected`.
- Full `testthat::test_local()` passing with **no existing test modified**.
- `lintr::lint_package()`: 0 lints.
