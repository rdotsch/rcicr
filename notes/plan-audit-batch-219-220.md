# Plan — closing out the audit batch: #219, then measuring #220

**For review before implementation.** This file is deleted in the PR that implements it, so the
plan does not outlive the work and become a second, drifting source of status: `NEWS.md`,
[`DECISIONS.md`](../DECISIONS.md) and the issues hold what survives.

## Context

The repo-wide static audit filed five issues, none of which its author ran against R. Three are
closed (#217, #218, #221). Two remain, and both are still hypotheses rather than confirmed facts.

**#219 (P2).** [`test-legacy-rdata.R`](../tests/testthat/test-legacy-rdata.R) reads real 1.0.1 and
1.1.0 `.Rdata` files and pins the missing-`nscales` fallback, but both fixtures use the default
**sinusoidal** noise. `sigma` only reaches the noise basis through `generateGabor()`
([`generateNoisePattern.R:45-49`](../R/generateNoisePattern.R#L45-L49)), so no sinusoidal fixture
can exercise the missing-`sigma` fallback at
[`generateReferenceDistribution.R:104-106`](../R/generateReferenceDistribution.R#L104-L106).
Changing that default from 25 leaves the whole suite green while silently changing the reference
distribution — and therefore the infoVal — for every legacy Gabor file.

Verified from the v1.0.1 tag rather than assumed: its `generateStimuli2IFC()` accepts
`noise_type`/`nscales`/`sigma` and passes them to `generateNoisePattern()`, and its explicit
`save()` call (`git show v1.0.1:R/generateStimuli2IFC.R`, line 171) writes `noise_type` but
neither `nscales` nor `sigma`. A 1.0.1 Gabor file is therefore the case that hits **both**
fallbacks at once, with `noise_type` correctly recorded as `gabor`.

Covering it exposes a second gap. `nscales` and `noise_type` both warn loudly when missing;
`sigma` defaults to 25 in silence — even though for a Gabor file a wrong `sigma` corrupts the
basis exactly the way a wrong `nscales` does. That asymmetry is fixed here *conditionally*, so
sinusoidal files, where `sigma` is provably inert, stay quiet.

**#220 (P3, low confidence).** `computeCumulativeCICorrelation()` does not aggregate repeated
stimuli the way `generateCI()` does ([`generateCI.R:184-191`](../R/generateCI.R#L184-L191)).
Almost certainly correct: the function walks trials in presentation order by design. The one
thing worth checking is its *self-computed* `finalCI`
([`computeCumulativeCICorrelation.R:112`](../R/computeCumulativeCICorrelation.R#L112)). This plan
**measures it and stops**, rather than proposing a change in advance of the number.

The issue's line citations (325/334/348) are stale — #224 rewrote that file, now 135 lines.

## Part 1 — #219: Gabor legacy fixture, tests, and a conditional warning

### 1a. Extend the fixture generator

[`tools/make-legacy-rdata.R`](../tools/make-legacy-rdata.R) maps tag → nscales via the `NSCALES`
vector ([line 49](../tools/make-legacy-rdata.R#L49)) and derives the output filename from
`DESCRIPTION`'s `Version` alone ([line 142](../tools/make-legacy-rdata.R#L142)).

- Replace `NSCALES` with a **`(tag, noise_type, nscales)` table**, one row per fixture: the
  existing `v1.0.1`/sinusoid/5 and `v1.1.0`/sinusoid/1 unchanged, plus `v1.0.1`/gabor/5.
- **`nscales = 5` and the era's default `sigma = 25` are load-bearing for the new row**, for the
  reason the existing comment already gives for `nscales`: 1.0.1 saves neither, so generating at
  the fallback values is what makes the fallbacks *correct*, and is the situation a returning
  researcher is actually in. Extend that comment rather than adding a second one.
- Derive the destination as `legacy-rdata-<version>[-<noise_type>].Rdata`, suffixing only
  non-sinusoidal rows so the two committed filenames are untouched.
- Keep the `die()` guard for a tag with no table entry.
- Pass `noise_type` and `sigma` explicitly into the generated `gen-<tag>.R` script.
- Several paths are keyed on `tag` alone (`worktree`, `lib`, `outdir`, `install-<tag>.log`); two
  rows now share `v1.0.1`, so key them on a per-row id instead.

Then `Rscript tools/make-legacy-rdata.R --versions=v1.0.1`, and commit the fixture.

**Primary execution risk:** nothing has yet confirmed v1.0.1's Gabor path *runs* at 32px. At
`nscales = 5` the smallest patch is 2px (`size <- patchSize[scale, img_size]`), so
`generateGabor(2, 1.5, ...)` must produce something usable. The existing sinusoidal fixture proves
the geometry works; only the patch generator differs. If it errors or yields a degenerate basis,
stop and report — do not quietly retune `nscales`, which would defeat the fixture's purpose.

Expect ~200 KB, matching the existing 1.0.1 sinusoid fixture (209 KB), under pre-commit's 1024 KB
`check-added-large-files` ceiling. Fixtures ship in the tarball, so confirm its size stays sane.

### 1b. Warn on a missing `sigma`, but only for Gabor files

In `generateReferenceDistribution2IFC()`
([`generateReferenceDistribution.R:104-106`](../R/generateReferenceDistribution.R#L104-L106)):
keep the silent default of 25, and warn **only when `noise_type` is `gabor`**. Word it like the
existing `nscales` warning — the assumed value, that rcicr did not save it before 1.1.0, and what
to do about it.

Ordering matters: the `noise_type` fallback sits *below* the `sigma` block
([lines 113-121](../R/generateReferenceDistribution.R#L113-L121)), so resolve `noise_type` first
or move the `sigma` check after it. A 1.0.1 Gabor file records `noise_type`, so it reaches the new
branch either way — but a file lacking both must not read an undefined `noise_type`.

Update the comment at [lines 107-112](../R/generateReferenceDistribution.R#L107-L112), which
describes `sigma` as taking "the silent default" and would otherwise be wrong in place.

### 1c. Tests

Extend [`test-legacy-rdata.R`](../tests/testthat/test-legacy-rdata.R), mirroring the `nscales`
tests and reusing `local_fixture_copy()`, `reference_norms_for()` and `mutate_rdata()`
([`helper-fixtures.R:53`](../tests/testthat/helper-fixtures.R#L53)) unchanged. `legacy_versions`
and `legacy_fixture()` are keyed on version strings and need a variant-aware form.

- The Gabor 1.0.1 fixture lacks `nscales` and `sigma` and carries `noise_type == "gabor"`.
- The current `generateCI()` reads it and returns a real CI.
- **The missing-`sigma` fallback picks 25 specifically**: the fallback null is `identical()` to an
  explicit `sigma = 25` and differs from an explicit `sigma = 10`. The *differs* half is what
  stops the test passing vacuously — finite norms would not, per CONTRIBUTING's "a test's title is
  a claim".
- The new warning fires for the Gabor file and does **not** fire for the sinusoidal 1.0.1 fixture.
  Match a fragment of the real `conditionMessage()`.

Prove the tests fail without the change with `git stash push -- R/` — never a plain `git stash`,
which reverts the new tests too and passes vacuously. The fixture and tests are both new, so also
check each new assertion individually rather than resting on the stash alone.

### 1d. `NEWS.md`

One entry under the development heading for the new warning. No numeric output changes, so **no
"Reproducibility impact" section**: the golden master and the release gate both stay green, since
the gate's battery generates files with versions that save `sigma` and never reaches this path.
Ordered per the largest-impact-first rule, it ranks below anything touching numbers.

## Part 2 — #220: measure only, then report

No files change in this part. Build a design presenting one stimulus an **unequal** number of
times, and compare:

1. `computeCumulativeCICorrelation()`'s self-computed `finalCI`
   ([line 112](../R/computeCumulativeCICorrelation.R#L112)), against
2. the CI `generateCI()` returns for the same data, which aggregates first
   ([`generateCI.R:184-191`](../R/generateCI.R#L184-L191)).

`generateCINoise()` takes `colMeans(stimuli * responses)`
([`generateCINoise.R:20-28`](../R/generateCINoise.R#L20-L28)), so with equal repeat counts the two
are algebraically identical and with unequal counts they diverge. Confirm that empirically and
quantify it — correlation and max absolute pixel difference between the two CIs, and the effect on
the curve's endpoint, which is what a user would actually notice. Check the equal-count case
agrees exactly, so the result distinguishes "only under unequal counts" from "always".

**Then stop and report the numbers**, before any documentation, [`DECISIONS.md`](../DECISIONS.md)
entry or test is written. What follows is the maintainer's call, resting on the measurement: a
doc note plus a `DECISIONS.md` entry recording a deliberate asymmetry if the difference is what
the issue expects, or a correction if the numbers show the self-computed `finalCI` is simply
wrong.

A correction is **not** ruled out in advance. The guiding constraint forbids changing numeric
output *silently*, not changing it — a documented fix carries a `NEWS.md` "Reproducibility
impact" entry and has to survive the release gate. What this step rules out is deciding either
way before the numbers exist.

## Verification

- `devtools::install()` first: tests calling `generateStimuli2IFC()` need the package installed,
  not just `load_all()`-loaded, because its parallel workers `library(rcicr)` themselves.
- `devtools::test()`, full suite green — `test-regression-baseline.R` and `test-legacy-rdata.R`
  matter most here.
- `Rscript tools/compare-release-output.R --quick`, to confirm the release gate is unaffected.
- A tarball build, then `R CMD check --as-cran` on it, watching for an installed-size NOTE from
  the added fixture. This is an ad-hoc build rather than a release one, so run it as:

  ```sh
  pkg=$(pwd)                                        # at the repo root
  tmp=$(mktemp -d) && (cd "$tmp" && R CMD build "$pkg")
  ```

  `R CMD build` takes the package directory as an argument and writes the tarball into the
  *working* directory, so this builds the real repo root — not a copy and not a worktree, so
  `.git` is a directory and `^\.git$` in `.Rbuildignore` behaves exactly as it does for a release
  — while the tarball lands outside the root. That matters because `R CMD build .` at the root
  overwrites any externally-checked tarball sitting there, which has happened. The release
  procedure is unchanged and still `R CMD build .` at the repo root
  ([`CONTRIBUTING.md`](../CONTRIBUTING.md) → "Releasing").
