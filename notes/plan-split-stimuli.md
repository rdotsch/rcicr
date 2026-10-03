# Plan: split generateStimuli2IFC() into its stages (#385)

Part of #384. `generateStimuli2IFC()` is 309 lines, cyclomatic complexity 45.

## The change

1. **`readBaseFace(label, filename, img_size, maximize_contrast)`**: the body of the base-image loop (read PNG/JPEG, square check, greyscale, size check, contrast). Same checks in the same order, same messages; the loop becomes `base_faces[[label]] <- readBaseFace(...)`.
2. **`drawStimulusParams(n_trials, nparams, labels, same)`**: one `matrix((runif(n_trials * nparams) * 2) - 1, n_trials, nparams, byrow = TRUE)` per parameter set, the shared set assigned to every label. Verified on #385 to give `identical()` parameters *and* an identical next random number.
3. **`renderStimulus(noise, base)`**: `((noise + 0.3) / 0.6 + base) / 2`, called with `trial_noise` and `-trial_noise` in the `foreach` body.
4. `1:n_trials` loops become `seq_len()`.

## Constraints kept

- `saveStimulusFile()` saves frame variables by name; every saved name stays bound in this frame (none of the loop temporaries moved into helpers is saved: checked against the list on the `saveStimulusFile()` line).
- The draw order, and so `seedResponseStream()`'s replay, is unchanged (point 2).
- Helpers called inside `%dopar%` are package functions, reachable through `.packages = 'rcicr'`.

## The one behaviour change, stated

`n_trials` is not validated. `0` or a negative number fails today with "subscript out of bounds", and `2.5` writes 2 trials into a file whose `n_trials` of 2.5 every InfoVal then rejects. `seq_len()` would change the first case silently, so the function checks `n_trials` up front (`validTrialCount()`, already used by the reference code) and stops with a clear message before anything is written. This is a NEWS bug-fix entry; every valid call is unchanged.

## Proving nothing changed

- **Full** release gate against `main` (`Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`): `0 expected deviations, 0 unexpected`.
- Full `testthat::test_local()` passing with **no existing test modified**.
- `lintr::lint_package()`: 0 lints. The gate compares `stimuli_params` exactly and every stimulus PNG by MD5, which is what the parameter and rendering helpers could break.

## Step most likely to fail

Point 2 inside the `use_same_parameters = FALSE` branch: each base draws its own block, in `names(base_faces)` order. The gate's `twobase-indep` configs pin that order.
