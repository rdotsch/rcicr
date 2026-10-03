# Plan: computeCumulativeCICorrelation() through the shared loaders (#387)

Part of #384.

## The change

Replace the function's own `load()` into its frame (with `captureArgs()`/`list2env()` to undo what `load()` overwrites), its three existence checks, the `s` → `p` conversion, the base lookup, the parameter selection and the 4096 → 4092 truncation with:

```r
loaded <- loadStimulusParams(rdata, require_img_size = FALSE)
base <- selectBaseImage(loaded$base_faces, baseimage)
params <- matrix(selectStimulusParams(loaded$stimuli_params, baseimage, stimuli), nrow = length(stimuli))
```

- `loadStimulusParams()` gains `require_img_size = TRUE`: this function never needed `img_size`, and making it required here would turn a call that works today into an error.
- `matrix(..., nrow = length(stimuli))` keeps the one-row matrix the cumulative loop's `params[1:trial, ]` relies on, where `selectStimulusParams()` returns a vector for a single stimulus.
- Loading into an environment of its own removes the argument-overwrite hazard instead of repairing it, so `captureArgs()` is no longer needed here.

## Messages

The base-label message is already identical in both. The three "did not contain" messages say "rdata argument" here and "rdata" in `loadStimulusParams()`; they become the latter. Every test matches on "did not contain X", which both satisfy. No NEWS entry: the wording loses one word and names the same field.

## Proving nothing changed

- **Full** release gate against `main` (`Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`): `0 expected deviations, 0 unexpected`.
- Full `testthat::test_local()` passing with **no existing test modified**.
- `lintr::lint_package()`: 0 lints. The gate's `cumulative` extra compares the curve.

## Step most likely to fail

A single-stimulus call: covered by an existing test (`drop = FALSE` regression) that must pass unchanged.
