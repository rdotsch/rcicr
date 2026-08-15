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

1. **Size check: wrong to assume every z-map is square.** An earlier revision of this plan
   argued `plotZmap(mask = zmap, img_size = nrow(zmap))` needed no generalisation because
   "every z-map is square (built from a square CI)" — true for z-maps `generateCI()` builds
   internally, false as a constraint on `plotZmap()` itself. It is a public, exported function;
   its roxygen for `mask` says only "the same dimensions as zmap", never square, and nothing in
   its body requires `nrow(zmap) == ncol(zmap))`. A caller passing a rectangular `zmap` with a
   matching rectangular mask works on current `main` (`nrow(zmap) == dim(mask)[1] && ncol(zmap)
   == dim(mask)[2]`, checked independently) and would wrongly error under a single scalar
   `img_size`, caught by review before it was written, not after.

   Checked rather than assumed how much actually needs to change: `applyMask()`'s existing
   comparison, `dim(mask_matrix) == img_size`, already works correctly for a length-2 `img_size`
   with no edit at all — R recycles a length-1 `img_size` against `dim()`'s length-2 result
   (today's scalar behaviour, unchanged) and compares element-wise against a length-2 one
   (verified directly: `applyMask(matrix(1,4,6), matrix(0,4,6), img_size = c(4,6))` returns a
   4×6 result on unmodified `main`). Only the size-mismatch *message* hardcodes `img_size` for
   both dimensions and needs a small change — see the shared-helper code below. `plotZmap()`
   passes `img_size = c(nrow(zmap), ncol(zmap))` instead of a bare scalar; `generateCI()`'s
   existing calls keep passing a scalar, unaffected.

2. **Error wording mentions "stimuli" unconditionally** in `applyMask()`'s size-mismatch
   message, which is inaccurate when called from `plotZmap()`. Resolution: `applyMask()` gains
   one new optional parameter, `context = 'stimuli'`, used everywhere the message currently
   hardcodes the word. `generateCI()`'s existing call sites pass nothing and get the same
   substantive message as today, with one incidental correction: the current text reads "...the
   stimuli! (stimulus dimensions: ..." — plural then singular, a pre-existing inconsistency
   (checked directly against `R/generateCI.R:461-462`) — which a single `context` value used in
   both places normalizes to "stimuli" throughout. Not byte-for-byte identical, then, but not a
   claim worth engineering around either: `tests/testthat/test-error-paths.R:159-168` only
   asserts the substrings `"same dimensions"`, `"8"`, `"4"`, none of which depend on the word,
   and the new parameter is additive so no existing call (positional or named) breaks.
   `plotZmap()` passes `context = 'z-map'`.

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

   **First revision, made in review:** generalising to "all `n` channels must match channel 1"
   (as an earlier draft of this plan proposed) is *also* wrong, caught by a second review
   comment. `applyMask()`'s current channels-1:3-only check is not incidental — for an RGBA
   mask whose RGB planes carry a genuine binary pattern and whose alpha plane is a constant
   (e.g. fully opaque, the common case for a mask authored in an external tool), `applyMask()`
   **succeeds today** by ignoring alpha, and real `generateCI()` scripts can depend on that.
   Verified directly: built exactly this array (RGB planes = a varying 0/1 pattern, alpha =
   constant 1) and confirmed `applyMask()` accepts it on current `main` while `plotZmap()`
   — using its all-channels-must-match loop — rejects the *same file*. So the two functions
   already disagree on this exact input; "all channels must match" would have resolved that
   disagreement by breaking the one call that currently works, not by fixing the drift.

   Resolution: alpha is not colour information, so it should never be compared. Generalise
   using PNG's own channel semantics instead of a generic "all channels equal" rule — a
   trailing channel is alpha whenever the total is even (2 = greyscale+alpha, 4 = RGBA); every
   other channel is colour and must agree with channel 1:

   ```r
   n <- dim(mask_matrix)[3]
   n_color <- if (n %in% c(2, 4)) n - 1L else n
   if (n_color > 1 && !all(sapply(2:n_color, function(i) {
     identical(mask_matrix[,,i], mask_matrix[,,1])
   }))) {
     stop(<existing "not a greyscale image" message>)
   }
   mask_matrix <- mask_matrix[,,1]
   ```

   Checked against six cases by direct execution, not just read through: 3-channel
   identical/differing (the two existing tests — unchanged results), 2-channel
   identical/differing (Codex's first finding — both now succeed, ignoring alpha), and
   4-channel with uniform RGB but differing alpha in both orders (Codex's second finding — the
   RGBA case still succeeds through `applyMask()` exactly as today, since alpha is dropped
   before the comparison; a 4-channel case where the RGB planes themselves differ still errors).
   No case that currently succeeds in either function changes result. What changes is only:
   `applyMask()` on 2-channel input goes from an unconditional crash to success (a bug fix,
   channel content does not matter — alpha was never compared for 3- or 4-channel input either);
   and `plotZmap()` on a 2-channel mask with a differing alpha plane goes from erroring to
   succeeding (previously-erroring path only, same category as points 2 and 4).

   Added to the test plan below: a 2-channel mask (identical planes, and differing planes) must
   succeed through `applyMask()` directly, and a 4-channel mask with uniform RGB but a differing
   alpha plane must also still succeed — the two cases review caught, pinned so neither can
   regress again.

4. **Type-check on a non-string, non-matrix mask.** `applyMask()` explicitly rejects it
   (`"neither a string nor a matrix"`, tested). `plotZmap()`'s current inline code has no such
   check — it calls `png::readPNG()` on whatever isn't `is.matrix()`, which fails with a
   low-level, unrelated error for something like `mask = TRUE` or `mask = list()`. No test pins
   `plotZmap()`'s current cryptic failure mode here. Resolution: once `plotZmap()` calls
   `applyMask()`, it inherits this clear error for free — a robustness improvement, not a change
   to any currently-succeeding call.

None of the four points changes the result of any call that currently succeeds through either
function. Point 3 changes two previously-*erroring* paths to succeed (a 2-channel mask crashing
`applyMask()` unconditionally; a 2-channel mask with a differing alpha plane erroring in
`plotZmap()`) while leaving every currently-succeeding call — including the 4-channel RGBA case
found in review — at its current result. Point 4 changes the error path for an input nothing in
the test suite or documented usage reaches today. This is not a numeric-output or call-syntax
change under the "guiding constraint" in `AGENTS.md`, but it is still an observable behaviour
change worth a `NEWS.md` entry, named explicitly rather than left for someone to notice.

## The change itself

`R/generateCI.R`'s `applyMask()` keeps its exact name, its exact three existing parameters
(`ci, mask, img_size = nrow(ci)`), and its exact behaviour for every case currently tested. It
gains one additive parameter and a generalised channel-collapse loop:

```r
applyMask <- function(ci, mask, img_size = nrow(ci), context = 'stimuli') {
  # ... same import/validate logic as today, with:
  #  - the channel-collapse block generalised per point 3's alpha-aware n_color rule
  #  - "stimuli" in the size-mismatch message replaced by the `context` argument (point 2)
  #  - the size-mismatch message's two hardcoded `img_size` reads replaced by
  #    `img_size[1]` and `img_size[length(img_size)]`, so it reports correctly whether
  #    img_size is a scalar (generateCI()'s calls, unchanged output) or length-2 (point 1)
}
```

`img_size = nrow(ci)` stays the default (unchanged — every `generateCI()` call site already
passes a scalar or relies on this default, both untouched).

`R/plotZmap.R`'s ~50-line inline `if (!is.null(mask)) { ... }` import/validate block (current
`main` lines ~126-177) is deleted and replaced with:

```r
if (!is.null(mask)) {
  zmap <- applyMask(zmap, mask, img_size = c(nrow(zmap), ncol(zmap)), context = 'z-map')
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
  existing call (three positional args plus `img_size =`) keeps its exact signature, and the
  substrings they assert (`"same dimensions"`, `"8"`, `"4"`, `"other than 0 or 1"`, `"not a
  greyscale image"`, `"neither a string nor a matrix"`) are all still present — see point 2 on
  the one incidental wording normalization that does not affect any of them.
- Add a test for the rectangular-mask case found during review: `applyMask()` (and so
  `plotZmap()`) accepts a mask whose dimensions match a non-square target passed as
  `img_size = c(rows, cols)`, and still reports both dimensions correctly (not the same value
  twice) when they do not match.
- Add one new test exercising `plotZmap()`'s newly-inherited "neither a string nor a matrix"
  guard directly (`plotZmap(mask = TRUE, ...)` or similar), since this is new observable
  behaviour for that function and nothing currently covers it.
- Add tests pinning the two channel-count cases found during review, all through `applyMask()`
  directly: a 2-channel (greyscale + alpha) mask succeeds whether the planes are identical or
  differing (alpha is always ignored); a 4-channel (RGBA) mask with uniform RGB planes and a
  differing alpha plane succeeds (matching current `main`'s behaviour exactly); a 4-channel mask
  whose RGB planes themselves differ still errors as non-greyscale — see point 3 above.

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
