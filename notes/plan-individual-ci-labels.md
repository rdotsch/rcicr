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

This entry carries more than the usual bug-fix note, because a researcher may be holding
mislabelled figures right now and cannot tell from the package which case they are in. It must
say three things a bare "fixed a labelling bug" would not:

**Which versions produced wrong files.** Every GitHub release from v1.0.1 to v1.2.3 inclusive,
and every install from `main` since 2017-08-15. Never any CRAN release — the last one, 0.3.4.1,
predates the defect — so `install.packages('rcicr')` never produced a mislabelled file.

**How to tell whether your own output is affected**, without re-running: compare
`unique(participants)` against `sort(unique(participants))` on the vector you passed. Identical
means every file was already correct. That is the whole test, and it is two expressions.

**That the damage is recoverable by renaming, not re-running.** The mapping is deterministic:
the file named `unique(participants)[i]` holds the CI of `sort(unique(participants))[i]`. This is
not asserted — it is what the reproduction measured, where `unique` was `c("b","a")`, `sort(unique)`
was `c("a","b")`, and `ci_b.png` held a's CI while `ci_a.png` held b's. So the entry gives the
recovery directly:

```r
old <- unique(participants)         # the names the files were given
new <- sort(unique(participants))   # the participants they actually hold
```

Group-level results, the returned `ci`, `pid_cis` and every z-map are unaffected, and the entry
must say so plainly — otherwise a reader assumes their whole analysis is in question when only
the individual-CI filenames ever were.

Wording follows the house rules: no version number in a `##` heading (it breaks `R CMD check`'s
news parsing), and largest-impact first within each section — this outranks the entries already
under "Reproducibility impact", since it is the one that can invalidate a published figure.

## The release gate covers this

An earlier draft of this plan argued the gate should stay out of it, on the grounds that the
battery is chosen by the reference version. That was wrong, and the Codex review on #262 caught
it. The reference-version restriction is specifically about calls where the reference
**crashes** — `mask`, undecorated z-maps, small z-maps — because a crash has no old value to
compare against. v1.0.1 does not crash here: its `generateCI()` carries `save_individual_cis`,
`targetpath` and `n_cores` (checked against `git show v1.0.1:R/generateCI.R`), so it runs this
path and writes the mislabelled files. A change that alters researchers' output while the gate
reports success is exactly what the gate exists to prevent.

So: a new `individual_cis` extra in `tools/compare-harness.R`, with **no `SINCE` floor**, since
nothing crashes.

- It calls `generateCI(participants = ..., save_individual_cis = TRUE, save_as_png = FALSE,
  targetpath = ..., n_cores = 1)` with participant IDs **deliberately out of sorted order** —
  the existing `participants` extra uses `rep(paste0("p", 1:4), ...)`, whose appearance order
  already equals its sorted order, so it could never show this.
- It captures `stats::setNames(unname(tools::md5sum(pngs)), basename(pngs))` over
  `<targetpath>/individual_cis`, reusing the idiom `stimulus_pngs` already uses at
  `compare-harness.R:206`. MD5-keyed-by-basename is the right shape because the defect *is* a
  filename-to-content mapping.
- Do **not** perturb the existing `participants` extra to make it unsorted. Its grouping of
  trials to participants would change with it, moving `participants_ci` and
  `participants_scaled` and manufacturing two unrelated deviations.

Both versions write the same *set* of filenames — the labels are permuted, not different — so
the keys align across runs and the deviation is a clean set of changed values rather than
missing keys.

Then an `EXPECTED` entry per config that lists the extra, `ref = "v1.0.1"`,
`news = "Reproducibility impact"`, whose reason states that the old bytes under each name
belonged to a different participant. Adding the extra to `ALL_EXTRAS` picks it up in
`defaults-512-sinusoid` automatically; the marginal cost there is writing four PNGs, since that
config already runs the `participants` extra.

## The step most likely to fail

**The `EXPECTED` entry firing for the wrong reason.** The gate fails a *stale* expectation as
well as an undocumented deviation, so if the individual-CI MD5s were to differ for some reason
beyond the relabelling — PNG encoding, scaling, anything — the entry would still be satisfied
and would mask it. Mitigation: confirm the deviating keys are exactly the permuted pair, and
that a sorted-ID run shows no deviation at all.

Second most likely: the PNG read-back in the test. `png::readPNG()` returns a 2D matrix for
greyscale and 3D for RGB(A) — the reproduction hit this — so the test must handle both rather
than indexing `[, , 1]` blind.

## Not in this change

#87 (autoscaling for individual CIs) rides on the same code and is the reason this was found.
It needs the per-participant save moved out of the `%dopar%` loop, since the shared constant is
not known until every participant's CI exists. Separate PR, on top of this one — this fix is a
one-expression correction that should be revertable on its own.
