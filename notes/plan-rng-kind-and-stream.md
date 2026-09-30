# Plan: record the RNG kind in the stimulus file, and give generateStimuli2IFC() back the caller's stream (#315, #189)

Both issues are about the global random stream, so they share one PR, one gate run and one Reproducibility impact entry.

## What was verified

Run here, R 4.3.3:

- Assigning a saved `.Random.seed` back to the global environment restores **both** the kind and the position. After `RNGkind("L'Ecuyer-CMRG"); set.seed(5); runif(3)`, saving `.Random.seed`, switching to Mersenne-Twister, reseeding and drawing, then assigning the saved vector back: `RNGkind()` reports L'Ecuyer-CMRG again and the next `runif(2)` is identical to the one the untouched stream gave.
- Doing the same with `RNGkind()` (save, switch, `do.call(RNGkind, old)`) does **not** give back the same draws (`FALSE`). This is the failed round trip #315 describes, and why the plan never uses it.
- `sample.kind` does not change `runif()`: `set.seed(1); runif(3)` gives identical results under Rejection and Rounding. `normal.kind` is not used anywhere on this path. Only `RNGkind()[1]` decides the stimuli and the simulated responses.
- Nothing the package computes draws after `generateStimuli2IFC()` returns without first reseeding: the golden master, the gate harness and `seedResponseStream()` all call `set.seed()` first (`tests/testthat/test-regression-baseline.R:57,102`, `tools/compare-harness.R:221,303`). So restoring the stream moves no number the gate compares.

## #315: record the kind, and replay references under it

1. **Record.** `generateStimuli2IFC()` saves a new `.Rdata` field `rng_kind <- RNGkind()`, the session's full three-element kind at generation. The stimuli are still drawn under the session's kind, so **no stimulus changes for anyone**. The file now says which kind its seed belongs to. The field is appended only (the `.Rdata` contract is append-only). No argument of any `load()`-ing function is called `rng_kind`, so it clobbers nothing; it goes in `globalVariables()` and in the README's `.Rdata` anatomy.
2. **Replay.** Every reference draw (`seedResponseStream()`: both the default stream and a `response_seed`) runs under the file's recorded uniform generator when there is one, via `set.seed(seed, kind = rng_kind[1])`.
   - **Same kind as the session** (the ordinary case): the code path is the current one, unchanged. The test that pins the caller's trailing stream after `generateReferenceDistribution2IFC()` (`test-generateReferenceDistribution2IFC.R:282-291`) still holds as written.
   - **Different kind**: draw under the file's kind, then assign the caller's saved `.Random.seed` back. That restores their kind and their position. Leaving the stream where the replay ended would switch the caller's session to another generator, which is worse than either option #315 rejected.
   - **Files without `rng_kind`** (every file written so far): unchanged behaviour. The session's kind is used, as documented today, and #314's "gave different values" message on an automatic refresh still covers the case where the kinds differ. The kind that built an old file cannot be recovered. Assuming Mersenne-Twister would be right for nearly every file but would change the reference for anyone who deliberately generated under another kind, so it is not assumed.
3. **Rejected: pinning one kind at generation.** It would make the seed alone determine the stimuli, but it changes the stimuli of anyone who generates under a non-default kind (for example L'Ecuyer-CMRG set for parallel work elsewhere in their script), and overrides a choice they made. Recording loses nothing. The stimuli are saved in the file anyway, so the seed is only needed to replay the response stream, and the recorded kind makes that replay exact.

A `user-supplied` kind is recorded like any other. Replaying it needs the same library loaded, and `set.seed()` errors otherwise, which is the right outcome.

## #189: generateStimuli2IFC() restores the caller's stream

Capture `.Random.seed` (or its absence) on entry, and on exit assign it back (or remove the one `set.seed()` created) under `on.exit()`. This is the helper `preserveRandomStream()` already uses for the automatic refresh; it moves to a new `R/random-stream.R` and both call sites share it. It restores on error and on interrupt too, so an aborted call no longer leaves the stream moved either.

## Tests

- New file: `rng_kind` is saved and equals `RNGkind()` at generation, under two different kinds.
- Replay: a file generated under L'Ecuyer-CMRG and scored under Mersenne-Twister gives the reference it gives under L'Ecuyer-CMRG, identically, for the default stream and for a `response_seed`. Also assert the **wrong** answer differs: the same file with `rng_kind` removed and scored under Mersenne-Twister gives a different reference.
- The caller's kind and position are unchanged after a cross-kind reference (`RNGkind()` and the next `runif()` compared with an untouched copy of the stream).
- The two #314 tests in `test-reference-from-saved-noise.R:785-836` currently assert that a refresh under another kind gives different values and says so. They describe files without `rng_kind`, so they run on a fixture with the field removed. A companion test asserts that a new file refreshes identically under either kind.
- #189: after `generateStimuli2IFC()`, `.Random.seed` is identical to before; with no `.Random.seed` before, none after; the same after an error mid-generation (the existing mocked-failure pattern in `test-stimulus-png-collision.R`).
- Mutation: removing the restore, and removing `kind =` from the replay, must each turn a test red.

## Reproducibility impact (NEWS.md)

- **#189**: after `generateStimuli2IFC()`, the next random draw in the caller's script differs from before, because the stream is returned where it was. Scripts that drew random numbers after generating stimuli without reseeding get different draws; anything the package computes is unchanged.
- **#315**: for files written from this version on, the InfoVal reference no longer depends on the scoring session's `RNGkind()`. Where it used to differ (scored under a different kind than generated), it now equals the reference under the generating kind. Existing files are unaffected.

`DECISIONS.md` line 34 ("Files do not record the RNG kind") is edited in place, with the recording-over-pinning reasoning in a few lines. The `?generateReferenceDistribution2IFC` Reproducibility section and `generateStimuli2IFC`'s `seed` documentation are updated to match.

## Gate

Expected: no deviation. `Rscript tools/compare-release-output.R --quick --ref="$(git rev-parse origin/main)"` should report `0 expected deviations` and no unexpected one: the harness draws everything under the default kind, where both changes are inert.

## The step most likely to fail

The cross-kind restore inside `generateReferenceDistribution2IFC()`. That function deliberately leaves the stream where the old loop did when the kinds match. Only the mismatch branch may restore it, and the two branches must not bleed into each other on error. Tested with a mocked failure inside the replay under a different kind.
