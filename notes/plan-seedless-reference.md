# Plan: stop when a stimulus file has no seed to replay (#334)

## Problem

Without a `response_seed`, the reference distribution replays the stimulus stream: `seedResponseStream()` (`R/reference-base.R`) calls `set.seed(source$seed)`. When the `.Rdata` file has no `seed`, that is `set.seed(NULL)`, which reseeds from the clock. The reference, and every InfoVal scored against it, then differs on every call, with no warning. #334 has the reproduction.

## Change

1. **Check for the seed before any simulation**, in both reference paths: `generateReferenceDistribution2IFC()` (shared parameters) and `generateBaseReference()` (independent base images), ahead of `referenceNoise()`, which at 512px is the slow step. One helper, called from both, stops when `response_seed` is `NULL` and the file has no `seed`. The message says the file has no stimulus seed to replay, and that `response_seed` draws a reproducible reference from the saved noise instead (not stored in the file).
2. **Only a missing `seed` is rejected.** Every other value `set.seed()` either already rejects with an error or seeds deterministically, so rejecting it would turn reproducible results into errors.
3. **`NEWS.md`, under "Bug fixes"**: who is affected (a stimulus file without `seed`, on 1.4.0 and 1.4.1), what they got (a different reference and InfoVal on every call), what happens now, and that a `response_seed` gives a reproducible reference. Files that have `seed` are unaffected.

No `DECISIONS.md` entry: this applies the existing "Trial alignment errors stop computation" reasoning, and the file has 2 words of headroom.

## Verified before planning

- `set.seed()` on each kind of bad value: `NULL` is the only one that runs and is not reproducible. `numeric(0)`, `NA`, `NaN`, `Inf`, `"abc"` and `list(1)` stop with "supplied seed is not a valid integer"; `"7"`, `c(1, 2)`, `1.5` and `TRUE` reproduce.
- Installed from their tags and ran #334's reproduction (each version generating its own file, then `seed` removed): v1.0.1 and v1.1.0 stop with "object 'seed' not found"; v1.4.1, the version on CRAN, and `main` return a reference that differs between two calls. `set.seed(source$seed)` arrived with #305 and is in v1.4.0 and v1.4.1, not v1.3.0.
- Every generator in the repository saves `seed`, back to the R-Forge import, and so do all three legacy fixtures. So the gate and every existing test use files with `seed`, and neither should move.
- `computeInfoVal2IFC()` passes `response_seed` through to the same paths, so the message's remedy applies to it too.

## Tests

- `generateReferenceDistribution2IFC()` on a file without `seed` stops, naming `seed` and `response_seed`, for shared and for independent base images;
- `computeInfoVal2IFC()` on such a file stops the same way;
- with a `response_seed`, a file without `seed` gives the same reference on two calls;
- the stop comes before the noise is built (checked by mocking `referenceNoise()` to fail if reached).

Each new failure test must fail on the current code; checked with `git stash push -- R/`.

## The step most likely to fail

**A cached reference that is refreshed automatically.** A file whose cached reference predates saved-noise references is rebuilt without being asked (`resolveReferenceNorms()`, "stale"). If such a file also lacks `seed`, `computeInfoVal2IFC()` used to return a new random InfoVal on each call and will now stop. That is the intended outcome, since the number could not be reproduced, but it is the one path where a call that returned a number now errors. A cached reference that does not need refreshing is still used, and no simulation runs.

## Out of scope

- #338 and #339: separate changes.
