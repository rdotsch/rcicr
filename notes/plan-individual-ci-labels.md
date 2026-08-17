# Plan — individual-CI PNGs carry the wrong participant's name (#261)

## The defect

`computeParticipantCIs()` indexes participants by sorted factor level and names the output
file by order of appearance. When those orderings disagree, every file in
`<targetpath>/individual_cis` carries the wrong participant's ID.

`R/ci-compute.R:24` — the loop index:

```r
pids <- as.numeric(factor(participants))     # sorted level order
```

`R/ci-compute.R:59` — the filename:

```r
saveToImage(..., unique(participants)[obs], antiCI)   # appearance order
```

## What I verified, and how

**The mislabelling is real end-to-end, not just index arithmetic.** Two participants with
maximally-different responses (`rep(1,6)` then `rep(-1,6)`) presented as `c("b","b",...,"a","a")`,
each written PNG compared pixelwise against that participant's CI computed on its own trials:

```
  ci_a.png:  vs a = 0.41906   vs b = 0.00196   -> holds b  <-- MISLABELLED
  ci_b.png:  vs a = 0.00196   vs b = 0.41813   -> holds a  <-- MISLABELLED
```

That 0.00196 is 8-bit PNG quantisation, 0.419 the genuine between-participant difference — a
separation of ~200x, which is what makes a content-based test viable through a PNG round-trip.

**Nothing `generateCI()` returns is affected.** The group CI is `apply(pid.cis, c(1,2), mean)`,
invariant to participant order; `pid_cis` is in factor-level order throughout, which is also the
order the `t.test` z-map consumes it in. Only the filenames are wrong.

**Present in every tagged release** — v1.0.1 through v1.2.3, checked with `git grep` at each
tag — introduced 2017-08-15 in 8a74974 at version 0.4.0. **Never shipped to CRAN**: the last
CRAN release is 0.3.4.1, which predates that commit, and before it the code wrote to a literal
`'targetpath/individual_cis<filename>'` so no file landed usefully at all.

**Why the existing test does not catch it.** `test-generateCI-paths.R:62` asserts only the set
of filenames, and uses `rep(c("a","b","c"), each = 4)` — already in sorted order, so the two
orderings coincide. A filename-only assertion cannot see this bug in any case.

## The fix

```r
sort(unique(participants))[obs]
```

Chosen over `levels(factor(participants))[obs]`, which is the same ordering: `factor()` levels
are `sort(unique(x))` coerced to character, verified `identical()` for `c(100000, 2, 10)`
including the `1e+05` formatting. `sort(unique())` preserves the caller's type, so numeric
participant IDs produce byte-identical filenames to today's — `paste0()` formats them the same
way either route. Nothing else in the loop reads the label.

Out of scope: `participants` containing some (not all) `NA` is separately broken — `pids == obs`
yields `NA` rows that index `params` with `NA`. Both the old and new expressions drop `NA` the
same way, so this change neither fixes nor worsens it.

## The test

`tests/testthat/test-fixed-bugs.R`, participants deliberately out of sorted order, asserting
**contents**, not filenames:

- read back `ci_a.png` and `ci_b.png`, compare each against that participant's CI computed
  standalone through `generateCI(scaling = "independent")`;
- tolerance `0.01` — comfortably above the 0.0039 (1/255) quantisation floor and far below the
  0.419 separation measured above;
- assert the two participants' CIs actually differ, so the test cannot pass vacuously;
- confirm it fails without the fix via `git stash push -- R/`, not by assumption.

## Reproducibility impact

Real, and it goes in `NEWS.md` under both "Reproducibility impact" and "Bug fixes", matching how
the masked-`targetci` fix is recorded. Someone who published a figure taken from `ci_p2.png`
finds different content under that name after this change — the *correct* content, but different.
Only unsorted-order data is affected; sorted-order data is byte-identical.

## The step most likely to fail

**The release gate will stay green, and that is not evidence.** `tools/compare-harness.R:246`
runs its `participants` configuration with `save_as_png = FALSE` and never sets
`save_individual_cis`, so the gate has never executed this write path. I am not adding a
configuration for it: the battery is chosen by the reference version, v1.0.1 has the bug, and a
new configuration would only manufacture an `EXPECTED` entry recording it. The suite test above
is the check that matters; the gate's silence here proves nothing either way.

Second most likely: the PNG read-back. `png::readPNG()` returns a 2D matrix for greyscale and
3D for RGB(A) — the reproduction hit this — so the test must handle both rather than indexing
`[, , 1]` blind.

## Not in this change

#87 (autoscaling for individual CIs) rides on the same code and is the reason this was found.
It needs the per-participant save moved out of the `%dopar%` loop, since the shared constant is
not known until every participant's CI exists. Separate PR, on top of this one — this fix is a
one-expression correction that should be revertable on its own.
