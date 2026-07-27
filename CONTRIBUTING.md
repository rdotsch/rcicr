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

## Larger changes

`BACKLOG.md` is the prioritized list of known work, with each item's evidence and any
decisions already taken — including things that look like bugs but are intentional, and
should not be "fixed". Read it before starting anything substantial, and if you are picking
up an item, say so on the issue first so the work is not duplicated.

## Code of conduct

Be decent to each other. Problems can be raised privately with the maintainer at the address
in `DESCRIPTION`.
