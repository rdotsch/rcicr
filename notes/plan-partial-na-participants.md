# Plan: stop on partly missing participant IDs (#337)

## Problem

`generateCI(participants = ...)` with some, but not all, IDs `NA` returns a group CI in which every pixel is `NA`. `computeParticipantCIs()` codes IDs with `as.numeric(factor(participants))`, which keeps `NA`; `pids == obs` is then `NA` for those trials, `params[pid.rows, ]` includes them as rows of `NA`, and every participant's CI, and so their average, becomes `NA`. #337 has the reproduction.

## Change

1. **`coerceTrialVectors()` (`R/ci-inputs.R`) stops** when `participants` is given (not all `NA`) and any entry is `NA`, after the existing length check. The message names how many trials lack an ID and which (the first few positions), and says to give every trial an ID or remove those trials from `stimuli`, `responses` and `participants` together. It runs before the `.Rdata` file is read, so nothing is computed or written.
2. **The all-`NA` convention is unchanged**: `NA`, or a vector of `NA`s, still means "no grouping".
3. **`?generateCI`**: the `participants` parameter says a partly missing vector is an error.
4. **`NEWS.md`, under "Bug fixes"**, after the #333 entry: what used to happen (a CI of `NA`s, a black PNG, `min`/`max` warnings), what happens now, and that no valid result changes.
5. **`DECISIONS.md`**, "Trial alignment errors stop computation", edited in place: its last sentence leaves "what a partly missing participant vector should mean" open; it now records the answer and why. Net word count stays within the budget, including #341's entry.

## Verified before planning

- Installed v1.0.1, v1.1.0 and `main`, and ran #337's reproduction with character and numeric IDs: every one returns a CI with 1024 of 1024 pixels `NA`. There is no version whose valid output this changes.
- A default call (`save_as_png = TRUE`) returns normally, raises 6 warnings (`no non-missing arguments to min`/`max`) and writes `ci_a.png`, which is solid black (every value 0). `ci`, `scaled` and `combined` are all `NA`.
- Only `generateCI()` takes `participants`. `batchGenerateCI()` passes `NA`; `generateCI2IFC()` has no such argument.
- The release gate's `participants` configuration (`tools/compare-harness.R:340`) uses complete IDs, so the gate should report no deviation.

## Tests

In `tests/testthat/test-ci-inputs.R`:

- `coerceTrialVectors()` stops for character, numeric and factor IDs with one `NA`, and the message gives the count and position;
- `generateCI()` with a partly missing vector stops, and writes nothing into `targetpath`;
- all-`NA` still means no grouping (existing test), and complete IDs give the same CI as before.

Each new test must fail on the current code; checked with `git stash push -- R/`.

## The step most likely to fail

**Choosing an error over dropping those trials.** Dropping them would return a plausible CI, but it decides for the user that a missing ID means "exclude", when it may be a data-entry error in a trial that belongs to a participant. The precedent is "Trial alignment errors stop computation" (#294, #300): a warning would still return a result built from a guess. The cost is that a script that used to run to a black PNG now stops, which is the point.

## Out of scope

- `NA` in `responses` or `stimuli`: `stimuli` is already validated; `responses` is a separate question, not measured here.
- #334, #338, #339: separate changes.
