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

### 1. Prepare the release PR

Branch from `main` (`release-X.Y.Z` by convention) and make these four edits together:

```sh
Rscript tools/compare-release-output.R              # both gate runs green
Rscript tools/compare-release-output.R --ref=vA.B.C # the previous release
R CMD build . && R CMD check --as-cran rcicr_*.tar.gz
```

- **`NEWS.md`**: rename `# rcicr (development version)` to `# rcicr X.Y.Z (YYYY-MM-DD)`.
  Until you do, none of this release's entries are in the news database at all — R indexes
  only sections under a heading it can read a version from. Keep version numbers *out* of
  `##` section headings, or the whole file stops parsing and `R CMD check` NOTEs.
- **`DESCRIPTION`**: drop the `.9000`. This is also what switches CI from the everyday
  `--quick` gate to the full battery, so it has to happen **before** the PR is opened, not
  after review.
- **`ChangeLog`**: add a dated pointer entry deferring to `NEWS.md` — not a duplicate. See
  `DECISIONS.md` for why.
- **`cran-comments.md`**: re-read every claim in it rather than trusting it. It has gone
  stale *within a single day* before, by describing a URL reference that another PR had
  just removed, one step short of a CRAN reviewer reading it.

Open the PR. The full battery runs on it, takes ~20 minutes, and is a required check.

### 2. Merge and tag

```sh
gh pr merge <n> --squash --delete-branch
git checkout main && git pull && git fetch --prune
git tag -a vX.Y.Z -m "rcicr X.Y.Z" && git push origin vX.Y.Z
gh release create vX.Y.Z --title "rcicr X.Y.Z" --notes-file <(...)   # NEWS section
```

The tag push re-runs the full gate against the released tree, which is the copy that
matters. Squash-merge, always: one commit per PR on `main`, and the squash message is where
measurements and rejected alternatives have to live, because the branch commits disappear.

### 3. Submit to CRAN

```sh
git switch --detach vX.Y.Z      # never build from main HEAD
R CMD build .
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
git switch -                    # back to where you were
```

Building from the tag is what keeps the development suffix out of the submitted version:
`Version contains large components` is only a CRAN blocker if the *tarball* carries it. Two
NOTEs are expected and explained in `cran-comments.md` — CRAN incoming feasibility (new
submission of an archived package) and future file timestamps (no network clock in a
sandbox).

Then run the two external checks and record their results in the `<!-- TODO -->` markers in
`cran-comments.md`. Neither is run casually: win-builder mails the maintainer, and both put
the tarball in front of a third party.

- **win-builder** — `devtools::check_win_devel()` and `devtools::check_win_release()`, run
  from a checkout of the tag. Each FTPs the tarball to `win-builder.r-project.org` and mails
  a check log to the maintainer address within the hour. The server **will not overwrite an
  existing upload**: a second attempt at the same filename fails with a bare
  `Failed FTP upload: 550`, which means "already queued", not "rejected". Confirm with
  `curl --list-only ftp://win-builder.r-project.org/R-devel/` before re-uploading, or you
  will conclude the upload failed when it succeeded.
- **R-hub** — `.github/workflows/rhub.yaml` is the stock R-hub v2 workflow, so the check runs
  on this repository's own Actions rather than uploading anywhere. Trigger it from the
  Actions tab ("R-hub" → Run workflow), which needs nothing but repository access.
  `rhub::rhub_check()` does the same over the API but requires a GitHub PAT with the `repo`
  scope, stored via `gitcreds::gitcreds_set()` — `rhub::rhub_doctor()` reports whether yours
  has it. Because the workflow is `workflow_dispatch`-only, it must be merged to the
  **default branch** before it can be triggered at all, and R-hub wants the file on the
  default branch *and* on whichever branch is being checked.

  **`RHUB_TOKEN` is deliberately unset, and nothing is missing.** The stock workflow passes
  `secrets.RHUB_TOKEN` to four actions, which makes it look like a prerequisite. It is an
  optional slot for your own PAT to reach private repositories — not a credential R-hub
  issues — and as of `r-hub/actions@v1` none of those four actions references the input in
  any step. This package is public; an unset secret expands to an empty string and the jobs
  run normally.

**Ron submits to CRAN personally.** CRAN emails the maintainer address to confirm, and for a
package archived over an undeliverable address, that confirmation *is* the point of the
resubmission.

### 4. Reopen development

Bump `DESCRIPTION` to `X.Y.Z.9000` and start a fresh `# rcicr (development version)` heading
in `NEWS.md`. `main` is protected, so this goes through a small PR like anything else.

`.github/workflows/reproducibility.yaml` runs the gate on every PR to `main` (`--quick`,
which skips the 512px configuration) and in full on release PRs, version tags and manual
dispatch. The `compare` job is a **required status check** on `main`, alongside
`R CMD check` on release and devel, so a red gate blocks the merge rather than merely
reporting. That is a repository ruleset, not something the workflow can enforce on itself.

## Larger changes

`BACKLOG.md` is the prioritized list of known work, with each item's evidence and any
decisions already taken — including things that look like bugs but are intentional, and
should not be "fixed". Read it before starting anything substantial, and if you are picking
up an item, say so on the issue first so the work is not duplicated.

## Code of conduct

Be decent to each other. Problems can be raised privately with the maintainer at the address
in `DESCRIPTION`.
