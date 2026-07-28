# Contributing to rcicr

Contributions, thoughts and criticisms are welcome. This file records the few conventions
that are specific to this package — the rest is ordinary R package practice.

## The one constraint that shapes everything else

**Researchers re-run old analysis scripts years later, and publish the numbers that come
out.** A change that silently alters what `generateCI()` or `computeInfoVal2IFC()` returns
can invalidate a published result without anyone noticing. So:

- Do not change existing call syntax, argument meanings, or numeric output silently.
  Deprecate rather than delete.
- Treat the `.Rdata` file's contents as **append-only**: add fields, never rename or
  repurpose them. (Its layout is documented in the README.)
- If a change *does* alter numeric output, that is not automatically wrong — but it must be
  deliberate, and it must be written up in `NEWS.md` under a "Reproducibility impact"
  heading saying who is affected and what to do about it.

`tests/testthat/test-regression-baseline.R` is a golden master pinning the default
pipeline's output. **If your change turns it red, that change alters researchers' results.**
Document it; do not update the baseline to match.

`tools/compare-release-output.R` asks the same question from the outside: it installs a
released version of the package from its own commit, runs both versions over a battery of
configurations, and compares. Every release has to pass it — see **Releasing** below.

## Getting set up

```r
devtools::install()   # not load_all() -- see below
devtools::test()
devtools::check()
```

`generateStimuli2IFC()` spawns parallel workers that each call `library(rcicr)`, so anything
touching it needs the package **actually installed**. Under `devtools::load_all()` alone the
workers fail with `there is no package called 'rcicr'`.

## Reporting a bug

The most useful report includes the `.Rdata` file's `img_size`, `nscales`, `noise_type` and
the rcicr version, plus a reproducible snippet. A synthetic base image works fine as a stand-
in — `png::writePNG(matrix(runif(32*32), 32, 32), "base.png")` — and is preferable to
attaching a real face photo.

## Pull requests

- **Branch from `main`**, one self-contained change per PR: the fix, its test, and its
  `NEWS.md` entry together. There is no `develop` branch.
- **Tests are expected.** Every pure function has unit tests; the I/O-heavy ones have lighter
  targeted tests. A bug fix should come with a test that fails without it — check this,
  rather than assuming, with `git stash push -- R/` (a plain `git stash` reverts your new
  tests too, so the suite passes vacuously and proves nothing).
- **A test's title is a claim; make the assertions actually support it.** Asserting that a
  result has the right shape is not the same as asserting it is correct. Where it is cheap,
  also assert that the *wrong* answer differs, so the test cannot pass vacuously.
- **`NEWS.md` entries go largest-impact first** within each section: changes to numeric
  output or return values, then behaviour changes, then bug fixes that only ever produced
  errors, then message-only fixes. Someone who stops reading after three bullets should have
  read the three that could change their results.
- **Any new top-level file must be added to `.Rbuildignore`** unless it genuinely belongs in
  the built package, or `R CMD check` raises a "non-standard file/directory found at top
  level" NOTE. Note the file holds `^`-anchored *regexes*, not globs.
- Run `git diff --stat main...HEAD` before opening the PR. `R CMD check` leaves a full copy
  of the package behind, and it has been committed by accident before.
- PRs are merged to `main` with **squash merges**, so put anything a future reader needs —
  measurements, rejected alternatives, reproducibility impact — in the PR description or
  `NEWS.md`, not only in individual branch commits.

CI runs `R CMD check` on the current R release and devel, reports coverage to Codecov, and
runs a small set of whitespace/YAML pre-commit hooks. `styler` and `lintr` are deliberately
**not** run: they would reformat nearly every file in one sweep and destroy `git blame`.

## Releasing

Releases are squash-merged onto `main` and marked with a tag; there is no `develop` branch.
The version in `DESCRIPTION` carries a `.9000` suffix between releases, and the release
commit is the one that drops it.

**Step 1 is the reproducibility gate, and a release does not go out while it is red.**

```sh
Rscript tools/compare-release-output.R                 # vs v1.0.1 -- the published baseline
Rscript tools/compare-release-output.R --ref=v1.1.0    # vs the last release
```

Both runs are required, and they answer different questions. The first asks whether this
tree still produces the numbers that are *in the literature* — it stays pinned at v1.0.1,
the last version on CRAN, and never advances, because a reference that moves forward with
each release lets a tree drift away from those numbers one tolerated epsilon at a time. The
second asks whether anything broke since the last release, and covers more ground: calls
that used to crash (masks, small z-maps, undecorated z-maps) have a comparable value in
v1.1.0 and none at all in v1.0.1.

Exit status is `0` for a clean run, `1` for a difference that is not accounted for, and `2`
if the comparison could not be run at all. On a difference there are exactly two honest
outcomes:

- **it was not intended** — fix it; or
- **it was intended** — write it up in `NEWS.md` under "Reproducibility impact", saying who
  is affected and what to do, and add an entry to `EXPECTED` in the script naming the
  reference, the affected output, the reason, and that `NEWS.md` heading. The script fails
  if an `EXPECTED` entry stops firing (so the list cannot rot) and if it fires while
  `NEWS.md` says nothing (so a deviation cannot be quiet).

Silently widening a tolerance, deleting a configuration from the battery, or skipping the
run is none of those things. The whole point is that a change to researchers' numbers can
only get out deliberately and loudly.

The rest of the checklist:

1. Reproducibility gate green (above), and `devtools::check()` clean.
2. Rename `NEWS.md`'s `# rcicr (development version)` heading to `# rcicr X.Y.Z (date)`.
3. Drop the `.9000` from `DESCRIPTION`. This is what makes CI run the *full* battery on the
   release PR rather than the everyday `--quick` one, so it has to happen before the PR is
   opened, not after it is approved.
4. Squash-merge, then `git tag -a vX.Y.Z <commit>` and cut a GitHub release.
5. **Build the CRAN tarball from the tag, never from `main` HEAD** — that is also what keeps
   the development suffix out of the submitted version.
6. Bump `DESCRIPTION` to `X.Y.Z.9000` and start a fresh `# rcicr (development version)`
   heading in `NEWS.md`.

`.github/workflows/reproducibility.yaml` runs the gate on every PR to `main` (`--quick`,
which skips the 512px configuration) and in full on release PRs, version tags and manual
dispatch. Note it only *blocks* a merge if it is listed as a required status check in the
repository's branch protection — the workflow cannot enforce that on its own.

## Larger changes

`BACKLOG.md` is the prioritized list of known work, with each item's evidence and any
decisions already taken — including things that look like bugs but are intentional, and
should not be "fixed". Read it before starting anything substantial, and if you are picking
up an item, say so on the issue first so the work is not duplicated.

## Code of conduct

Be decent to each other. Problems can be raised privately with the maintainer at the address
in `DESCRIPTION`.
