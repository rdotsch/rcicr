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

Two edits, not one. The second came out of the first review round.

**1. `R/plotZmap.R:136`** — the guard the issue asks for:

```r
if (!is.null(mask)) {        ->        if (hasMask(mask)) {
```

`hasMask()` is an internal function in the same package, already reachable from `plotZmap.R`
exactly as `applyMask()` is; no `NAMESPACE` or `man/` change.

**2. `R/generateCI.R:421`** — `hasMask()` cannot currently tell the scalar sentinel from a
one-cell matrix, so adopting it in `plotZmap()` unmodified would spread that gap rather than
just align the two. Require the sentinel to be dimensionless:

```r
hasMask <- function(mask) {
  !is.null(mask) &&
    !(length(mask) == 1L && is.na(mask) && is.null(dim(mask)))
}
```

Measured across the inputs that matter — `NULL`, `NA`, `NaN`, `NA_character_` and `array(NA)`
stay `FALSE` (no mask); `matrix(NA, 1, 1)` flips `FALSE` -> `TRUE`; `matrix(1, 1, 1)` and
`matrix(NA, 2, 2)` were already `TRUE` and stay `TRUE`.

## What this changes for users, and what it cannot

**From edit 1 — nothing that works today.** `plotZmap(mask = NA)` and `plotZmap(mask = NaN)`
stop erroring and render an unmasked z-map. Both error today, so the change can only turn a
hard error into a render. `NaN` is not in the issue and falls out of `hasMask()`'s `is.na()`
test — worth stating in `NEWS.md` rather than leaving to be discovered.

**From edit 2 — one call that succeeds today becomes an error**, and this is the part with real
user-facing surface. Measured on the standard fixture:

```
generateCI(..., mask = matrix(NA, 1, 1))   ->  OK, any NA in ci: FALSE
```

It silently returns an *unmasked* CI: a mask argument the user plainly meant as a mask is
discarded without warning. That is the silent-wrong-output class this repo ranks above crashes,
and it becomes `This mask contains values other than 0 or 1!`. Erroring is the safe direction —
the alternative is aligning `plotZmap()` onto the weaker behaviour — but it is a behaviour
change to a second function beyond the issue's literal ask, so it gets its own `NEWS.md` entry
under **Behaviour changes**, above the `plotZmap` one.

For `plotZmap()` the same input errors today either way (size mismatch against a larger z-map;
`This mask contains values other than 0 or 1!` when the sizes match), so edit 2 changes nothing
there — it only stops edit 1 from making it silent.

**No numeric output moves.** Every affected call either errors today or starts erroring; none
crosses from one successful numeric result to a different one. No golden-master or release-gate
exposure. `mask = NULL`, a PNG path, a correctly-sized matrix, `matrix(1, 1, 1)` and
`c(NA, NA)` are untouched.

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

## Rejected alternatives

**Make `plotZmap()` strict *and* `generateCI()` strict** — reject the scalar `NA` in both. This
is the defensible reading of the issue in the other direction, and consistency-wise it is just
as aligned. It is not possible: `NA` is `generateCI()`'s documented default, so erroring on it
would break every call that omits `mask`. For the scalar sentinel, alignment can only go the
permissive way.

Worth naming because it is the genuine residual cost: silently ignoring a scalar `NA` means a
mask variable that is accidentally `NA` yields an unmasked z-map with no complaint. That is
already true of `generateCI()`, and matching it is the point. Edit 2 confines the silence to the
sentinel, where it is unavoidable, instead of extending it to malformed matrices, where it is
not.

**Take the one-liner alone and file the `hasMask()` gap separately.** Rejected: `plotZmap()`
would then join `generateCI()` in silently ignoring `matrix(NA, 1, 1)`, so this PR would widen a
silent-wrong-output bug in the course of closing a consistency issue. Fixing the guard and the
predicate together is one coherent change; splitting them ships the regression first.

## Steps

1. Edit `R/plotZmap.R:136` and `hasMask()` at `R/generateCI.R:421`.
2. Extend `tests/testthat/test-plotZmap.R` with the three-way render comparison above, and add
   a `matrix(NA, 1, 1)` case on both `plotZmap()` and `generateCI()` asserting the validation
   error rather than a silent pass.
3. Confirm each test fails without its fix, via `git stash push -- R/` (not a plain `git stash`,
   which would take the new tests with it and pass vacuously).
4. Two `NEWS.md` entries under **Behaviour changes**: the `generateCI()` one first, since it
   changes a call that currently succeeds; the `plotZmap()` one second, naming `NaN` as well as
   `NA`.
5. Check `test-regression-baseline.R` and the release gate stay green — expected, since no
   successful numeric result moves, but expected is not measured.
6. Delete this file, mark the draft ready, answer the Codex review, squash.
