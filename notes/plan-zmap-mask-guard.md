# Plan — align `plotZmap()`'s mask guard with `generateCI()`'s (#246)

## What is actually true today

Measured, not assumed. `devtools::load_all()`, an 8x8 z-map, `decoration = FALSE`:

| call | today | `hasMask()` says |
|---|---|---|
| `plotZmap(mask = NA)` | **error**: `The mask argument is neither a string nor a matrix!` | `FALSE` |
| `plotZmap(mask = NaN)` | **error**: same message | `FALSE` |
| `plotZmap(mask = NULL)` | no mask applied, renders | `FALSE` |
| `plotZmap(mask = c(NA, NA))` | **error**: same message | `TRUE` |

The guard is `!is.null(mask)` at `R/plotZmap.R:136`; `generateCI()` uses `hasMask(mask)` at
`R/generateCI.R:295` and `:321`, defined at `:421`.

**Why the two diverged, which the issue does not say:** the defaults differ. `generateCI()` has
`mask = NA` (`R/generateCI.R:108`), `plotZmap()` has `mask = NULL` (`R/plotZmap.R:112`). Each
guard is correct for its own default, and `hasMask()` exists precisely because `generateCI()`'s
default is `NA` — see its comment: a bare `!is.na(mask)` returns a whole matrix when a mask *is*
supplied, which R >= 4.2 rejects.

That `generateCI(mask = NA)` silently means "no mask" needs no new measurement: `NA` is its
default, so every default-mask call in the test suite already demonstrates it.

## The change

One line, `R/plotZmap.R:136`:

```r
if (!is.null(mask)) {        ->        if (hasMask(mask)) {
```

`hasMask()` is an internal function in the same package, already reachable from `plotZmap.R`
exactly as `applyMask()` is; no `NAMESPACE` or `man/` change.

## What this changes for users, and what it cannot

`plotZmap(mask = NA)` and `plotZmap(mask = NaN)` stop erroring and render an unmasked z-map.
`NaN` is not in the issue and falls out of `hasMask()`'s `is.na()` test — worth stating in
`NEWS.md` rather than discovering later.

**No working call changes.** Both affected inputs error today, so nothing that currently
succeeds can behave differently: the change can only turn a hard error into a render. No
numeric output, no reproducibility impact, no golden-master or release-gate exposure. `mask =
NULL`, a matrix, a PNG path, and `c(NA, NA)` are all untouched.

## The step most likely to fail

Not the edit — the test. `plotZmap()` returns nothing and writes a PNG, so "did the mask get
applied" is only observable through rendered pixels, and `CONTRIBUTING.md` → Testing forbids
asserting absolute pixel values (they belong to the graphics device). The assertion must be a
*relationship between renders*: `mask = NA` must produce the same image as `mask = NULL`, and
both must differ from a real half-mask. `tests/testthat/test-plotZmap.R` already has the
`render_zmap()` helper and a uniform background for exactly this, so the test extends existing
machinery rather than inventing any.

Second risk, smaller: a test asserting only "`mask = NA` does not error" would pass vacuously if
the guard were wrong in the other direction. Hence the comparison against both a no-mask render
and a masked one.

## Rejected alternative

**Make `plotZmap()` strict *and* `generateCI()` strict** — reject `NA` in both. This is the
defensible reading of the issue in the other direction, and consistency-wise it is just as
aligned. It is not possible: `NA` is `generateCI()`'s documented default, so erroring on it
would break every call that omits `mask`. Alignment can only go the permissive way.

Worth naming because it is the genuine cost of this change: silently ignoring `NA` means a
typo'd or accidentally-`NA` mask variable now yields an unmasked z-map with no complaint. That
is already true of `generateCI()`, and matching it is the point.

## Steps

1. Edit `R/plotZmap.R:136`.
2. Extend `tests/testthat/test-plotZmap.R` with the three-way render comparison above.
3. Confirm the test fails without the fix, via `git stash push -- R/` (not a plain `git stash`,
   which would take the new test with it and pass vacuously).
4. `NEWS.md` entry under **Behaviour changes**, naming `NaN` as well as `NA`.
5. Delete this file, mark the draft ready, answer the Codex review, squash.
