# Plan: stop when a stimulus file has no seed to replay (#334)

## Problem

Without a `response_seed`, the reference distribution replays the stimulus stream: `seedResponseStream()` (`R/reference-base.R`) calls `set.seed(source$seed)`. When the `.Rdata` file has no `seed`, or one whose value is `NULL`, that is `set.seed(NULL)`, which reseeds from the clock. `generateStimuli2IFC(seed = NULL)` is accepted and writes such a file (named `rcic_seed__time_...`). Measured on `main`:

| path | result |
|---|---|
| `computeInfoVal2IFC()`, shared base images, no stored reference, `seed` field absent | stops, with the `dplyr` error "In argument: `ref_seed == seed`" from the empty `ref_lookup` |
| the same call on a file made with `seed = NULL` | passes the lookup (it matches nothing) and returns a different InfoVal on two fresh copies of the file |
| `computeInfoVal2IFC(force_gen_ref_dist = TRUE)`, or `generateReferenceDistribution2IFC()` | a different reference on every call |
| `computeInfoVal2IFC()`, first automatic build (independent base images, or a shared reference refreshed as stale) | a random reference, stored; later calls reuse it, so the InfoVal is stable on that file but differs on a fresh copy |

## Change

1. **Stop before any simulation when there is no `seed` and no `response_seed`.** "No `seed`" means the field is absent or `NULL`, tested with `is.null()` on the loaded value, never `exists()`. One helper, called ahead of `referenceNoise()` (the slow step at 512px) in both reference paths, `generateReferenceDistribution2IFC()` and `generateBaseReference()`, and at the top of `computeInfoVal2IFC()`'s lookup block, so the first row gets this message instead of the `dplyr` one. The message says the file has no stimulus seed to replay, and names `generateReferenceDistribution2IFC(rdata, response_seed = <n>)` as the way to store a reproducible reference, adding `baseimage = "<label>"` on the independent-base path, where that call requires it (`selectReferenceBase()`). The stored reference records its seed, and later `computeInfoVal2IFC()` calls reuse it. For a file that cannot be written, it also names `computeInfoVal2IFC(..., response_seed = <n>)`, which scores reproducibly without writing (a seeded draw is never stored from that call); storing needs a writable copy. The message makes no other storage claim, since whether a seeded reference is saved depends on the call.
2. **Only a missing `seed` is rejected.** Every other value `set.seed()` either already rejects with an error or seeds deterministically, so rejecting it would turn reproducible results into errors.
3. **Warn when a stored reference was drawn from the clock.** A file without `seed` (absent or `NULL`) whose stored reference has no `response_seed` got it from a clock-seeded draw. It is still used, so no number changes, but `computeInfoVal2IFC()` warns that it cannot be reproduced from the stimuli, and how to replace it (the call in 1).
4. **`NEWS.md`, under "Reproducibility impact"**, since forced and direct calls used to return numbers and now stop: the three rows above, for 1.4.0 and 1.4.1; what happens now; how to get a reproducible reference; and that files with `seed` are unaffected.

No `DECISIONS.md` entry: this applies the existing "Trial alignment errors stop computation" reasoning, and the file has 2 words of headroom.

## Verified before planning

- The table above, on `main` (script: two copies of one seedless file, each path called twice), for a file with the field removed and for one made with `generateStimuli2IFC(seed = NULL)`.
- `set.seed()` on each kind of bad value: `NULL` is the only one that runs and is not reproducible. `numeric(0)`, `NA`, `NaN`, `Inf`, `"abc"` and `list(1)` stop with "supplied seed is not a valid integer"; `"7"`, `c(1, 2)`, `1.5` and `TRUE` reproduce.
- Installed from their tags: v1.0.1 and v1.1.0 stop on a seedless file with "object 'seed' not found"; v1.4.1, the version on CRAN, returns a different reference on two generator calls. `set.seed(source$seed)` arrived with #305 and is in v1.4.0 and v1.4.1, not v1.3.0.
- `generateReferenceDistribution2IFC(rd, response_seed = 7)` on a seedless file stores `reference_norms` with `reference_norms_seed = 7`, and two `computeInfoVal2IFC()` calls then reuse it ("Using reference distribution found in rdata file") and return identical values.
- Every generator in the repository saves `seed`, back to the R-Forge import, and so do all three legacy fixtures, so the gate and the existing tests should not move.

## Tests

- each row of the table now stops with the new message, before `referenceNoise()` is reached (mocked to fail if called), with fixtures of both kinds: `seed` removed, and made with `seed = NULL`;
- on a file that cannot be written (`writableFile()` mocked), `computeInfoVal2IFC(response_seed = <n>)` returns the same value twice and writes nothing;
- the message's call runs as written, for shared and for independent base images, and `computeInfoVal2IFC()` then reuses the stored reference with identical values;
- a seedless file carrying a stored reference without `response_seed` returns the same InfoVal as before, with the warning; one with a recorded `response_seed` gives no warning.

Each new failure test must fail on the current code; checked with `git stash push -- R/`.

## The step most likely to fail

**The third row.** Its users got a stable number and will now get an error on a fresh copy, or a warning on the file that holds it. Both are intended: that number cannot be reproduced from the stimuli. The warning is also the one change that reaches files where nothing is recomputed, so its condition must be exactly "`seed` absent or `NULL`, and no recorded `response_seed`", never a file with a real `seed`.

## Out of scope

- #338 and #339: separate changes.
