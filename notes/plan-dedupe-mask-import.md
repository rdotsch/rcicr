# Plan: extract one mask-import helper, shared by generateCI() and plotZmap()

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
that has always been live"). So the direction is: extract `applyMask()`'s import/validate logic
into a shared helper, and make `plotZmap()` use it — not the reverse.

## Behavioural differences found while reading both, and how each is resolved

Read both implementations line by line (`R/generateCI.R:421-480`, `R/plotZmap.R:125-177` on
current `main`) rather than assuming they agree. Four differences, each checked against the
test suite for whether resolving it toward `applyMask()`'s behaviour changes anything currently
tested or reachable through documented usage:

1. **Size check.** `applyMask()` compares `dim(mask_matrix)` against a scalar `img_size` (an
   implicit square-target assumption already baked into every caller — CIs are always square).
   `plotZmap()` compares `dim(mask)` against `nrow(zmap)`/`ncol(zmap)` directly. Resolution: the
   shared helper takes `target_dims = c(nrow, ncol)` instead of a scalar, so both callers keep
   their exact current comparison (`applyMask` passes `c(img_size, img_size)`) with no behaviour
   change for any square target — which is every real one.

2. **Error wording mentions "stimuli" unconditionally** in `applyMask()`'s size-mismatch
   message. Not meaningful from `plotZmap()`. Resolution: the helper takes a `context` string
   (`'stimuli'` default, so `applyMask()`'s call site and its existing message text are
   byte-for-byte unchanged; `plotZmap()` passes `'z-map'`).
   `tests/testthat/test-error-paths.R:159-168` only asserts the substrings `"same dimensions"`,
   `"8"`, `"4"` — none of which depend on the context word — so this is safe against the pinned
   assertions.

3. **Multi-channel collapsing.** `applyMask()` hardcodes three channels
   (`list(mask_matrix[,,1], mask_matrix[,,2], mask_matrix[,,3])`), so it throws `subscript out
   of bounds` — not a clean error — on any PNG with `dim()[3] != 3`. `plotZmap()` instead loops
   `for (i in 2:dim(mask)[3])`, which works for any channel count including 2 (an 8-bit
   grayscale-plus-alpha PNG, which `png::readPNG()` decodes to a 2-channel array). Verified
   directly rather than assumed: on current `main`,
   `plotZmap(mask = <2-channel PNG, planes identical>, decoration = FALSE, ...)` **succeeds**,
   while the equivalent `applyMask()` call on the same file throws `subscript out of bounds` —
   confirmed by sourcing both files and calling each against a real 2-channel PNG written with
   `png::writePNG()`. So naively adopting `applyMask()`'s hardcoded-3 logic, as an earlier draft
   of this plan did, would have broken a call that works today on `plotZmap()` — caught by
   review, not by the earlier draft's own (insufficiently checked) reasoning.

   Resolution: neither existing implementation is right to copy as-is. The shared helper
   generalises to arbitrary channel count `n = dim(mask_matrix)[3]`: collapse when every channel
   is identical to channel 1 (`all(sapply(2:n, function(i) identical(mask_matrix[,,i],
   mask_matrix[,,1])))`), else the existing "not a greyscale image" error. This is equivalent to
   `plotZmap()`'s pairwise-consecutive loop by transitivity (so its currently-working
   2-channel and 3-channel calls are unaffected) and equivalent to `applyMask()`'s current
   behaviour for the tested 3-channel case (same three channels, same comparison, same result).
   For 4-channel (RGBA) input it is *stricter* than `applyMask()`'s current hardcoded check
   (which silently ignores channel 4): now the alpha channel must match too, or the mask is
   rejected as non-greyscale rather than silently accepted with alpha discarded. Untested either
   way (`test-error-paths.R:194-217` only covers 3-channel RGB) — reject-on-mismatch is the
   safer default for previously-unreachable input, matching how every other guard in this helper
   already errs toward rejecting an ambiguous mask rather than guessing.

   Added to the test plan below: a 2-channel (greyscale + alpha) mask with identical planes must
   still succeed through both `applyMask()` and `plotZmap()`, and one with differing planes must
   still error as non-greyscale through both — the exact case review caught, pinned so it cannot
   regress again.

4. **Type-check on a non-string, non-matrix mask.** `applyMask()` explicitly rejects it
   (`"neither a string nor a matrix"`, tested). `plotZmap()` has no such check — it calls
   `png::readPNG()` on whatever isn't `is.matrix()`, which fails with a low-level, unrelated
   error for something like `mask = TRUE` or `mask = list()`. No test pins `plotZmap()`'s
   current cryptic failure mode here. Resolution: `plotZmap()` gains `applyMask()`'s clear
   error for this case — a robustness improvement, not a change to any currently-succeeding
   call.

None of the four changes alters the result of any call that currently succeeds or that any
test currently exercises. Two (#3, #4) change the *error path* for inputs nothing in the test
suite or documented usage reaches today. This is not a numeric-output or call-syntax change
under the "guiding constraint" in `AGENTS.md` — it changes what happens for previously-broken
inputs — but it is still an observable behaviour change worth a `NEWS.md` entry under
"Internal", named explicitly rather than left for someone to notice.

## The shared helper

New file `R/importMask.R`:

```r
importMask <- function(mask, target_dims, context = 'stimuli') {
  # ... exact logic of the current applyMask() import/validate block, generalised
  # per points 1-2 above, plus the channel-count generalisation from point 3
  # (collapse when every channel equals channel 1, for any n = dim()[3] >= 2),
  # returning a logical matrix (TRUE = masked, i.e. mask == 0).
}
```

`R/generateCI.R`'s `applyMask(ci, mask, img_size = nrow(ci))` keeps its exact existing name and
signature (tests call it directly with `img_size =`) and becomes a thin wrapper:

```r
applyMask <- function(ci, mask, img_size = nrow(ci)) {
  bool_mask <- importMask(mask, target_dims = c(img_size, img_size))
  ci[bool_mask] <- NA
  return(ci)
}
```

`R/plotZmap.R`'s inline block becomes:

```r
if (!is.null(mask)) {
  mask <- importMask(mask, target_dims = c(nrow(zmap), ncol(zmap)), context = 'z-map')
}
```

leaving the existing `if (!is.null(mask)) { zmap[mask] <- NA }` below it untouched (already
expects a boolean matrix in that convention).

## Test changes

- `tests/testthat/test-plotZmap.R`'s `"mismatched mask and zmap dimensions error"` test
  currently asserts the message matches `"not the same size"` (plotZmap's own old wording).
  Update the expected pattern to `"same dimensions"` (the now-shared wording) — this is
  asserting our own error text, not a numeric contract.
- No other existing test asserts `plotZmap()`'s previous mask-import wording (checked via
  `grep` across `tests/testthat/`), so nothing else needs updating.
- `tests/testthat/test-error-paths.R`'s four `applyMask()` guard tests need no change: the
  wrapper preserves the exact call signature and exact message text for every case they cover.
- Add one new test exercising `plotZmap()`'s newly-shared "neither a string nor a matrix" guard
  directly (`plotZmap(mask = TRUE, ...)` or similar), since this is new observable behaviour for
  that function and nothing currently covers it.
- Add a test pinning the 2-channel (greyscale + alpha) case found during review: a mask PNG with
  identical grayscale and alpha planes succeeds (through both `applyMask()` and `plotZmap()`),
  one with differing planes errors as non-greyscale (through both) — see point 3 above.

## Verification before marking ready

- `testthat::test_local()` — full suite green, including the updated and new assertions above.
- `lintr::lint_package()` — zero new lints (the `.lintr` baseline from #183 must not need
  regenerating; if it does, that is a signal something more than intended changed).
- `tools/compare-release-output.R --quick` — the release gate. Checked what its `mask` scenario
  actually calls (`tools/compare-harness.R:271-278`): only `generateCI(mask = ...)`, never
  `plotZmap(mask = ...)`. So the gate covers `applyMask()`'s path but not the shared helper's
  use from `plotZmap()` — that side's regression coverage is `tests/testthat/test-plotZmap.R`
  alone, which is why its existing mask tests (same-half, all-zero, matrix-vs-PNG-agree,
  dimension-mismatch) all have to keep passing unchanged, not just the new one.
- `NEWS.md` gets an "## Internal" entry describing the dedup and naming the two resolved
  behaviour differences (points 3-4 above) explicitly.
