# Plan: de-duplicate mask-import logic between generateCI() and plotZmap()

Issue #185. Originally reported as #89 (closed) — the duplication is still there and has
already produced divergent behaviour once: `plotZmap(mask = ...)` validated its mask and then
never applied it (fixed for 1.1.0; see `DECISIONS.md` → "`plotZmap(mask = ...)` was applied
rather than deprecated"). The two implementations have kept drifting since.

## Where the duplication lives today

- `R/generateCI.R`'s `applyMask(ci, mask, img_size)` (internal, tested directly via `rcicr:::`
  in `tests/testthat/test-error-paths.R` and `test-plotZmap.R`).
- `R/plotZmap.R`'s inline `if (!is.null(mask)) { ... }` block inside `plotZmap()` itself
  (lines ~125-177 on current `main`).

Per `DECISIONS.md` → "When the docs and the code disagreed about `mask`, the code was the
contract": `applyMask()` is the copy that has "been executing for a decade" and is the one
`tests/testthat/test-error-paths.R` calls the canonical guard tests against (comment there:
"plotZmap() has one of these covered; generateCI() has none -- and generateCI() is the copy
that has always been live"). So the direction is: `plotZmap()` starts calling `applyMask()`
directly — not the reverse, and not a third function in between either.

**Revision, made in review:** the first draft of this plan invented a separate `importMask()`
helper returning a boolean matrix, with `applyMask()` reduced to a thin wrapper around it. Ron
asked directly on the PR why not just call `applyMask()` from `plotZmap()` and skip the extra
layer — and there is no reason not to. `applyMask(ci, mask, img_size)` does not care what the
matrix it masks represents; it imports, validates, and sets masked cells to `NA` on whatever
numeric matrix is passed in. `plotZmap()` can call `applyMask(zmap, mask, img_size = nrow(zmap))`
exactly as `generateCI()` calls `applyMask(ci, mask, img_size)`. One existing function, one new
caller, no new file.

## Behavioural differences found while reading both, and how each is resolved

Read both implementations line by line (`R/generateCI.R:421-480`, `R/plotZmap.R:125-177` on
current `main`) rather than assuming they agree. Three differences (a fourth, the size-check
shape, disappears entirely under the revision above — see below), each checked against the test
suite for whether resolving it toward `applyMask()`'s behaviour changes anything currently
tested or reachable through documented usage:

1. **Size check no longer needs generalising.** The first draft worried that `applyMask()`
   compares against a scalar `img_size` while `plotZmap()` compares `nrow(zmap)`/`ncol(zmap)`
   separately, and proposed a `target_dims` vector to reconcile them. That concern only existed
   because of the abandoned `importMask()` indirection. Calling `applyMask(zmap, mask, img_size =
   nrow(zmap))` directly needs no such change: `img_size` is just "the dimension to check the
   mask against", and every z-map is square (built from a square CI), so `nrow(zmap)` is exactly
   the right scalar to pass. `applyMask()`'s own signature and size-check logic are untouched.

2. **Error wording mentions "stimuli" unconditionally** in `applyMask()`'s size-mismatch
   message, which is inaccurate when called from `plotZmap()`. Resolution: `applyMask()` gains
   one new optional parameter, `context = 'stimuli'`. `generateCI()`'s existing call sites pass
   nothing and get byte-for-byte the same message as today; `plotZmap()` passes `context =
   'z-map'`. `tests/testthat/test-error-paths.R:159-168` only asserts the substrings `"same
   dimensions"`, `"8"`, `"4"` — none of which depend on the context word — so this is safe
   against the pinned assertions, and the new parameter is additive so no existing call
   (positional or named) breaks.

3. **Multi-channel collapsing.** `applyMask()` hardcodes three channels
   (`list(mask_matrix[,,1], mask_matrix[,,2], mask_matrix[,,3])`), so it throws `subscript out
   of bounds` — not a clean error — on any PNG with `dim()[3] != 3`. `plotZmap()` instead loops
   `for (i in 2:dim(mask)[3])`, which works for any channel count including 2 (an 8-bit
   grayscale-plus-alpha PNG, which `png::readPNG()` decodes to a 2-channel array). Verified
   directly rather than assumed: on current `main`,
   `plotZmap(mask = <2-channel PNG, planes identical>, decoration = FALSE, ...)` **succeeds**,
   while the equivalent `applyMask()` call on the same file throws `subscript out of bounds` —
   confirmed by sourcing both files and calling each against a real 2-channel PNG written with
   `png::writePNG()`. So routing `plotZmap()` through `applyMask()` unchanged would have broken
   a call that works today — caught by review, not by the earlier draft's own reasoning.

   Resolution: fix `applyMask()`'s channel-collapse in place, since it is now the one
   implementation both functions share. Generalise to arbitrary channel count `n =
   dim(mask_matrix)[3]`: collapse when every channel is identical to channel 1
   (`all(sapply(2:n, function(i) identical(mask_matrix[,,i], mask_matrix[,,1])))`), else the
   existing "not a greyscale image" error. This is equivalent to `plotZmap()`'s
   pairwise-consecutive loop by transitivity (so its currently-working 2-channel and 3-channel
   calls are unaffected) and equivalent to `applyMask()`'s current behaviour for the tested
   3-channel case (same three channels, same comparison, same result) — a genuine bug fix for
   `generateCI()` too: a 2-channel mask PNG crashes it today, and will not after this change.
   For 4-channel (RGBA) input the new logic is *stricter* than today's hardcoded check (which
   silently ignores channel 4): the alpha channel must now match too, or the mask is rejected as
   non-greyscale. Untested either way (`test-error-paths.R:194-217` only covers 3-channel RGB) —
   reject-on-mismatch is the safer default for previously-unreachable input, matching how every
   other guard in this function already errs toward rejecting an ambiguous mask over guessing.

   Added to the test plan below: a 2-channel (greyscale + alpha) mask with identical planes must
   succeed through `applyMask()` directly (and so through both callers), and one with differing
   planes must error as non-greyscale — the exact case review caught, pinned so it cannot
   regress again.

4. **Type-check on a non-string, non-matrix mask.** `applyMask()` explicitly rejects it
   (`"neither a string nor a matrix"`, tested). `plotZmap()`'s current inline code has no such
   check — it calls `png::readPNG()` on whatever isn't `is.matrix()`, which fails with a
   low-level, unrelated error for something like `mask = TRUE` or `mask = list()`. No test pins
   `plotZmap()`'s current cryptic failure mode here. Resolution: once `plotZmap()` calls
   `applyMask()`, it inherits this clear error for free — a robustness improvement, not a change
   to any currently-succeeding call.

None of the four points changes the result of any call that currently succeeds through either
function, except point 3's 2-channel case, which changes `generateCI()` from crashing to working
(a bug fix) and leaves `plotZmap()`'s already-working behaviour unchanged. Points 3 (RGBA) and 4
change the *error path* for inputs nothing in the test suite or documented usage reaches today.
This is not a numeric-output or call-syntax change under the "guiding constraint" in
`AGENTS.md`, but it is still an observable behaviour change worth a `NEWS.md` entry, named
explicitly rather than left for someone to notice.

## The change itself

`R/generateCI.R`'s `applyMask()` keeps its exact name, its exact three existing parameters
(`ci, mask, img_size = nrow(ci)`), and its exact behaviour for every case currently tested. It
gains one additive parameter and a generalised channel-collapse loop:

```r
applyMask <- function(ci, mask, img_size = nrow(ci), context = 'stimuli') {
  # ... same import/validate logic as today, with:
  #  - the channel-collapse block generalised to arbitrary n = dim(mask_matrix)[3] (point 3)
  #  - "stimuli" in the size-mismatch message replaced by the `context` argument (point 2)
}
```

`R/plotZmap.R`'s ~50-line inline `if (!is.null(mask)) { ... }` import/validate block (current
`main` lines ~126-177) is deleted and replaced with:

```r
if (!is.null(mask)) {
  zmap <- applyMask(zmap, mask, img_size = nrow(zmap), context = 'z-map')
}
```

This also folds in the separate `zmap[abs(zmap) < threshold] <- NA` / `zmap[mask] <- NA`
two-step that currently follows the inline block — `applyMask()` already does the "set masked
cells to NA" step internally, so `plotZmap()` no longer needs its own trailing `zmap[mask] <-
NA`. Sequencing is unchanged: threshold is applied first, then the mask, exactly as today.

## Test changes

- `tests/testthat/test-plotZmap.R`'s `"mismatched mask and zmap dimensions error"` test
  currently asserts the message matches `"not the same size"` (plotZmap's own old wording).
  Update the expected pattern to `"same dimensions"` (the now-shared wording) — this is
  asserting our own error text, not a numeric contract.
- No other existing test asserts `plotZmap()`'s previous mask-import wording (checked via
  `grep` across `tests/testthat/`), so nothing else needs updating.
- `tests/testthat/test-error-paths.R`'s four `applyMask()` guard tests need no change: every
  existing call (three positional args plus `img_size =`) keeps its exact signature and exact
  message text for every case they cover, since `context` is additive with a default.
- Add one new test exercising `plotZmap()`'s newly-inherited "neither a string nor a matrix"
  guard directly (`plotZmap(mask = TRUE, ...)` or similar), since this is new observable
  behaviour for that function and nothing currently covers it.
- Add a test pinning the 2-channel (greyscale + alpha) case found during review: a mask PNG with
  identical grayscale and alpha planes succeeds through `applyMask()` directly, one with
  differing planes errors as non-greyscale — see point 3 above.

## Verification before marking ready

- `testthat::test_local()` — full suite green, including the updated and new assertions above.
- `lintr::lint_package()` — zero new lints (the `.lintr` baseline from #183 must not need
  regenerating; if it does, that is a signal something more than intended changed).
- `tools/compare-release-output.R --quick` — the release gate. Checked what its `mask` scenario
  actually calls (`tools/compare-harness.R:271-278`): only `generateCI(mask = ...)`, never
  `plotZmap(mask = ...)`. So the gate covers `applyMask()`'s import/validate logic (now shared)
  but not `plotZmap()`'s specific call site — that side's regression coverage is
  `tests/testthat/test-plotZmap.R` alone, which is why its existing mask tests (same-half,
  all-zero, matrix-vs-PNG-agree, dimension-mismatch) all have to keep passing unchanged, not
  just the new ones.
- `NEWS.md` gets an entry describing the dedup, the 2-channel bug fix, and the two resolved
  error-path differences (points 3-4 above) explicitly.
