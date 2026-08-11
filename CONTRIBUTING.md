# Contributing to rcicr

Contributions, thoughts and criticisms are welcome. This file records
the few conventions that are specific to this package — the rest is
ordinary R package practice.

## The one constraint that shapes everything else

**Researchers re-run old analysis scripts years later, and publish the
numbers that come out.** A change that silently alters what
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
or
[`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
returns can invalidate a published result without anyone noticing. So:

- Do not change existing call syntax, argument meanings, or numeric
  output silently. Deprecate rather than delete.
- Treat the `.Rdata` file’s contents as **append-only**: add fields,
  never rename or repurpose them. (Its layout is documented in the
  README.)
- If a change *does* alter numeric output, that is not automatically
  wrong — but it must be deliberate, and it must be written up in
  `NEWS.md` under a “Reproducibility impact” heading saying who is
  affected and what to do about it.

`tests/testthat/test-regression-baseline.R` is a golden master pinning
the default pipeline’s output. **If your change turns it red, that
change alters researchers’ results.** Document it; do not update the
baseline to match.

`tools/compare-release-output.R` asks the same question from the
outside: it installs a released version of the package from its own
commit, runs both versions over a battery of configurations, and
compares. Every release has to pass it — see **Releasing** below.

## Getting set up

``` r

devtools::install()      # not load_all() -- see below
devtools::test()
devtools::check()
roxygen2::roxygenise()   # after editing any roxygen comment; commit man/ and NAMESPACE
```

`man/` and `NAMESPACE` are generated and tracked, and
`ubuntu-latest (release)` regenerates them and fails on a difference.
Use the roxygen2 version in `DESCRIPTION`’s `RoxygenNote`; the CI step
is pinned to it and says so when they disagree.

[`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
spawns parallel workers that each call
[`library(rcicr)`](https://rdotsch.github.io/rcicr/), so anything
touching it needs the package **actually installed**. Under
`devtools::load_all()` alone the workers fail with
`there is no package called 'rcicr'`.

**An `--as-cran` run needs a complete check environment, or it reports
failures the package does not have.** Install `texlive` and `tidy`, and
set both incoming variables:

``` sh
_R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE \
  R CMD check --as-cran rcicr_X.Y.Z.tar.gz
```

Without the toolchain a run here reported 1 ERROR + 1 WARNING + 4 NOTEs
that were entirely the sandbox — no `pdflatex`, no `tidy`, and a
leftover `rcicr-manual.tex` from the failed PDF build. Installing it was
the only change needed to reach 2 NOTEs. **Do not reach for
`--no-manual` to make the manual checks go away**: it skips them rather
than passing them, and a note here once recorded that as the problem
being *resolved*. A clean run shows
`checking PDF version of manual ... OK` and
`checking HTML version of manual ... OK` — if you do not see those
lines, you have not checked the manual.

## Reporting a bug

The most useful report includes the `.Rdata` file’s `img_size`,
`nscales`, `noise_type` and the rcicr version, plus a reproducible
snippet. A synthetic base image works fine as a stand- in —
`png::writePNG(matrix(runif(32*32), 32, 32), "base.png")` — and is
preferable to attaching a real face photo.

## Pull requests

- **Branch from `main`**, one self-contained change per PR: the fix, its
  test, and its `NEWS.md` entry together. There is no `develop` branch.
- **Tests are expected.** Every pure function has unit tests; the
  I/O-heavy ones have lighter targeted tests. A bug fix should come with
  a test that fails without it — check this, rather than assuming, with
  `git stash push -- R/` (a plain `git stash` reverts your new tests
  too, so the suite passes vacuously and proves nothing).
- **A test’s title is a claim; make the assertions actually support
  it.** Asserting that a result has the right shape is not the same as
  asserting it is correct. Where it is cheap, also assert that the
  *wrong* answer differs, so the test cannot pass vacuously.
- **`NEWS.md` entries go largest-impact first** within each section:
  changes to numeric output or return values, then behaviour changes,
  then bug fixes that only ever produced errors, then message-only
  fixes. Someone who stops reading after three bullets should have read
  the three that could change their results.
- **Any new top-level file must be added to `.Rbuildignore`** unless it
  genuinely belongs in the built package, or `R CMD check` raises a
  “non-standard file/directory found at top level” NOTE. Note the file
  holds `^`-anchored *regexes*, not globs.
- Run `git diff --stat main...HEAD` before opening the PR. `R CMD check`
  leaves a full copy of the package behind, and it has been committed by
  accident before.
- PRs are merged to `main` with **squash merges**, so put anything a
  future reader needs — measurements, rejected alternatives,
  reproducibility impact — in the PR description or `NEWS.md`, not only
  in individual branch commits.
- **Write for the state the change ends in, not the route you took to
  it.** Commit messages, PR comments and code comments describe what the
  change is and why it is right. How many attempts it took, and what
  each one got wrong, is invisible in the squashed history and in the
  merged code — nobody arrives at the file that way. Keep a rejected
  alternative only when someone would otherwise retry it, and give it a
  line, not a chronology.

### Testing

- **`skip_if_not_installed()` is for `Suggests` only.** `withr` may
  genuinely be absent, so skipping on it is honest. A package in
  `Imports` — `png`, `jpeg`, `matlab` — cannot be: the package will not
  load without it. The skip can never fire for a real reason, and if it
  somehow did it would hide the test rather than fail it. A skip reads
  as “this passed” at a glance, which is the opposite of what you want
  to know.
- **Going beyond the one-off check above — mutating the code repeatedly
  — keep `cp` backups in a scratchpad, and never restore with
  `git checkout <file>`.** The `git stash push -- R/` in the checklist
  is right for proving a single fix; it is the *restore* step that
  bites, because `git checkout <file>` discards unstaged work and has
  destroyed an in-progress implementation here. Guard each mutation with
  a `grep -q MUTANT` check too — one that silently failed to apply looks
  exactly like a surviving mutant.
- **When a test reads pixels back from a graphics device, assert only
  relationships between renders.** Every absolute property of those
  pixels belongs to the device: the channel count (cairo writes RGB,
  macOS quartz writes RGBA) and the values alike (quartz renders a 0.5
  background at roughly 0.573 where cairo gives 0.502). Render onto a
  *uniform* background so “drew nothing” becomes “the image is one flat
  value”, count distinct values over the colour channels only, and
  compare two renders rather than pinning a number.
  [`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.html#pixel-assertions-have-measured-the-graphics-device-twice)
  has the two failures that taught this.
- **Check figures by looking at them.** Three problems in the
  walkthrough vignette were invisible to every assertion and obvious on
  sight; they are listed in
  [`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.html#three-vignette-figures-were-wrong-in-ways-only-viewing-them-showed).
- **Base images in tests and vignettes are always synthetic**, never a
  real photograph. Avoiding licensing and consent questions entirely is
  worth more than a realistic-looking figure.

### The Codex review

Codex reviews pull requests here and has caught real errors (#170,
\#172, \#173, and three on the PR that wrote this section). Nothing in
the merge path makes you notice it: it submits as `COMMENTED`, so it
never blocks and `reviewDecision` stays empty; `gh pr checks` is green
regardless; and `gh pr view --comments` shows only the review wrapper,
never the findings.

**It must not become something that can block.** If it has been switched
off, or simply is not answering, merge on the other checks.

Two conditions clear a squash.

1.  **Push everything, then trigger, keeping the timestamp.** A push
    does not re-trigger the review — only this comment, opening the PR,
    or marking a draft ready does. Posting it returns the anchor you
    need:

    ``` sh
    trig=$(gh api repos/rdotsch/rcicr/issues/<n>/comments -f body='@codex review' --jq '.created_at')
    ```

2.  **Wait for a 👍 newer than `$trig`.** The reaction on the PR body
    tracks the *latest run* rather than the PR’s history: 👀 while it
    runs, 👍 when it finishes with nothing to say, and no reaction at
    all when it has findings.

    ``` sh
    gh api repos/rdotsch/rcicr/issues/<n>/reactions \
      --jq '.[] | select(.user.id == 199175422) | "\(.user.login) \(.content) \(.created_at)"'
    ```

    Only a `+1` dated after `$trig` clears it. 👀, nothing, and an older
    `+1` all mean not cleared. The filter is Codex’s numeric account id
    (`chatgpt-codex-connector[bot]`, printed back so you can see whose
    reaction cleared the gate) because reactions are open to anyone with
    read access: a substring match on the login would let any account
    containing “codex” clear a review that is still running, and an id
    survives a rename besides. When not cleared, read the findings and
    answer each one:

    ``` sh
    gh api --paginate repos/rdotsch/rcicr/pulls/<n>/comments --jq '.[] | "\(.path): \(.body)"'
    ```

3.  **Resolve every thread you have answered.** The 👍 speaks only for
    the latest run, so a finding left unanswered from an earlier round
    does not stop a later clean one. This half is **enforced rather than
    checked**: the `main` ruleset sets
    `required_review_thread_resolution`, so an unresolved thread blocks
    the squash server-side and no client-side count can clear it by
    mistake. Resolve one by passing its `id` to the
    `resolveReviewThread` mutation.

    To read what is still open — for answering, not for deciding:

    ``` sh
    gh api graphql -f query='query { repository(owner: "rdotsch", name: "rcicr") {
      pullRequest(number: <n>) { reviewThreads(first: 100) { nodes { id isResolved path } } } } }' \
      --jq '.data.repository.pullRequest.reviewThreads.nodes[] | select(.isResolved == false)'
    ```

    That reads the first 100 threads and does not paginate, which is why
    it must not be the gate: on a PR with more it would report nothing
    while an unresolved thread sat on page two. GitHub is the gate; this
    is a convenience.

`--paginate` on the findings listing is not optional: that endpoint
pages at 30, and a review-heavy PR truncates silently without it. The
filter there is per-element, so it prints every finding whichever way
the pages arrive.

One green state, everything else not green: anything other than a 👍
newer than your trigger sends you to read the findings, whose worst case
is a wasted look. Earlier versions of this section computed completion
from review objects and `commit_id` and were wrong six separate ways —
reporting a clean run as unfinished, an unrelated bot comment as clean,
and hiding a reply because replies inherit the parent thread’s commit.
The lesson each time was the same: a client-side predicate that answers
“is this safe to merge” fails open. The reaction is the one signal Codex
sets deliberately, and thread resolution is enforced by GitHub rather
than counted here.

CI runs `R CMD check` on the current R release and devel, reports
coverage to Codecov, and runs a small set of whitespace/YAML pre-commit
hooks. `styler` and `lintr` are deliberately **not** run: they would
reformat nearly every file in one sweep and destroy `git blame`.

## Code conventions

These describe what the code already does. **They apply to new and
modified code**; there is no expectation that anyone reformats untouched
files to match, and no sweep is planned — see issue \#194 for the one
that would be, and why it needs a major version.

**Two of these are frozen by the constraint above, not chosen:**

- **Exported functions are camelCase** —
  [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md),
  [`batchGenerateCI2IFC()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md),
  [`plotZmap()`](https://rdotsch.github.io/rcicr/reference/plotZmap.md).
  This is *not* what a modern R style guide would say, and it cannot be
  changed: these names are in researchers’ stored scripts and in
  published methods sections.
- **Arguments and list fields are snake_case** — `img_size`, `n_trials`,
  `noise_type`, `save_as_png`. Also frozen; people call them by name.
  Note this means arguments already match tidyverse style and only
  function names diverge.

The rest are ordinary consistency, and internal helpers are free to
change because nothing outside the package can call them:

|  | Convention | Notes |
|----|----|----|
| Assignment | `<-`, never `=` | Already universal — there are zero `=` assignments in `R/`. |
| Internal helpers | camelCase, matching the exported style | Currently 7 camelCase (`applyMask`, `saveToImage`, `startBackend`, …) against one snake_case (`default_ncores`). New helpers follow the majority. |
| Strings | single quotes | Roughly 3:2 in favour today; not worth churn to unify, but write new code single-quoted. |
| Indentation | 2 spaces, no tabs | Already consistent; there are no tabs in `R/`. |
| Booleans | `TRUE`/`FALSE`, never `T`/`F` | `T` and `F` are rebindable variables. In package code they resolve through the namespace and never reach a user’s globals, so this is style rather than a hazard — but it is free, and the reasoning has to be re-derived every time someone notices. All occurrences were replaced in 1.2.2. |
| Sequences | [`seq_len()`](https://rdrr.io/r/base/seq.html)/[`seq_along()`](https://rdrr.io/r/base/seq.html), not `1:n` | `1:0` counts *backwards*, so `1:length(x)` on an empty vector iterates twice. |
| Returns | explicit [`return()`](https://rdrr.io/r/base/function.html) at the end of exported functions | The existing style throughout. |
| Files | one file per exported function, named after it | `R/generateCI.R` holds [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md). `zzz.R` holds the [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html) declarations. |
| Namespacing | prefer `@importFrom pkg fn` over `@import pkg` | `@import matlab` masks [`base::sum()`](https://rdrr.io/r/base/sum.html) with MATLAB semantics across six files — a live trap, issue \#182. |

Line length is not enforced; 35 lines in `R/` already exceed 100
characters. Wrap new code at something reasonable rather than reflowing
what is there.

**Comments: as few as will do, and only where the code cannot speak for
itself.** Write one when the reason would otherwise have to be
re-derived — a rejected alternative, a constraint, a measured number, a
trap. Never narrate what the next line does, restate the error message
below it, or re-tell an issue’s history. The long comments already in
`R/` are the first kind and earn their length; the test is whether a
reader would otherwise get it wrong, not whether the code looks bare.
Where the explanation is about the package rather than about the line,
it belongs in
[`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.md).

**If the package is ever run through `styler`, it goes in as a commit of
its own** — never as a side effect of other work, so `git blame`
survives it. Why the pre-commit hooks stay minimal and language-agnostic
is above.

Three rules that are about this package rather than about R:

- **Do not write a package-qualified call as code — in backticks or
  `\code{}` — for a package the docs only *mention*.** The pkgdown site
  resolves such a link by loading that package, and a package that is
  installed but cannot load takes the whole site build down with it.
  This is not hypothetical: after `raster` was dropped in \#186 the
  runners kept it in their restored library while the sysreqs scan
  stopped installing its GDAL/PROJ system libraries, so the build died
  with `libproj.so.25: cannot open shared object file` — once from the
  roxygen, then again from the `NEWS.md` entry describing the same
  change. Name it in prose instead: “the raster package’s plot method”.
  A package that is simply absent is fine, which is why `rhub::` and
  `gitcreds::` appear above without trouble. `README.md`, `NEWS.md`,
  `DECISIONS.md`, this file and the vignettes are all rendered.
- **Add new names loaded from an `.Rdata` file to `R/zzz.R`’s
  [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)**,
  or `R CMD check` NOTEs about undefined globals.
- **Check any new function argument against the names saved in the
  `.Rdata`**, and any new saved field against existing argument names.
  [`load()`](https://rdrr.io/r/base/load.html) assigns into the calling
  frame, so a collision silently overwrites the argument — this has
  caused three separate bugs, most recently a saved `sigma` capturing
  [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)’s
  z-map blur `sigma`.
  [`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.html#load-assigns-into-the-calling-frame--check-every-new-argument-against-saved-names)
  has the details.

## Releasing

Releases are squash-merged onto `main` and marked with a tag; there is
no `develop` branch. The version in `DESCRIPTION` carries a `.9000`
suffix between releases, and the release commit is the one that drops
it.

**Step 1 is the reproducibility gate, and a release does not go out
while it is red.**

``` sh
Rscript tools/compare-release-output.R                 # vs v1.0.1 -- the published baseline
Rscript tools/compare-release-output.R --ref=v1.1.0    # vs the last release
```

Both runs are required, and they answer different questions. The first
asks whether this tree still produces the numbers that are *in the
literature* — it stays pinned at v1.0.1, the last version on CRAN, and
never advances, because a reference that moves forward with each release
lets a tree drift away from those numbers one tolerated epsilon at a
time. The second asks whether anything broke since the last release, and
covers more ground: calls that used to crash (masks, small z-maps,
undecorated z-maps) have a comparable value in v1.1.0 and none at all in
v1.0.1.

Exit status is `0` for a clean run, `1` for a difference that is not
accounted for, and `2` if the comparison could not be run at all. On a
difference there are exactly two honest outcomes:

- **it was not intended** — fix it; or
- **it was intended** — write it up in `NEWS.md` under “Reproducibility
  impact”, saying who is affected and what to do, and add an entry to
  `EXPECTED` in the script naming the reference, the affected output,
  the reason, and that `NEWS.md` heading. The script fails if an
  `EXPECTED` entry stops firing (so the list cannot rot) and if it fires
  while `NEWS.md` says nothing (so a deviation cannot be quiet).

Silently widening a tolerance, deleting a configuration from the
battery, or skipping the run is none of those things. The whole point is
that a change to researchers’ numbers can only get out deliberately and
loudly.

### 1. Prepare the release PR

Branch from `main` (`release-X.Y.Z` by convention) and make these five
edits together:

``` sh
Rscript tools/compare-release-output.R              # both gate runs green
Rscript tools/compare-release-output.R --ref=vA.B.C # the previous release
R CMD build . && R CMD check --as-cran rcicr_*.tar.gz
```

- **`NEWS.md`**: rename `# rcicr (development version)` to
  `# rcicr X.Y.Z (YYYY-MM-DD)`. Until you do, none of this release’s
  entries are in the news database at all — R indexes only sections
  under a heading it can read a version from. Keep version numbers *out*
  of `##` section headings, or the whole file stops parsing and
  `R CMD check` NOTEs.
- **`DESCRIPTION`**: drop the `.9000`. This is also what switches CI
  from the everyday `--quick` gate to the full battery, so it has to
  happen **before** the PR is opened, not after review.
- **`ChangeLog`**: add a dated pointer entry deferring to `NEWS.md` —
  not a duplicate. See
  [`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.html#changelog-gets-a-pointer-entry-not-a-duplicate)
  for why.
- **`CITATION.cff`**: regenerate it, after the `DESCRIPTION` edit above,
  with
  `cffr::cff_write("DESCRIPTION", dependencies = FALSE, gh_keywords = FALSE)`
  — it carries the version. `ubuntu-latest (release)` regenerates and
  compares, so it fails until you do.
- **`cran-comments.md`**: re-read every claim in it rather than trusting
  it. It has gone stale *within a single day* before, by describing a
  URL reference that another PR had just removed, one step short of a
  CRAN reviewer reading it. Leave the test-environment results for step
  2; everything else is written now.

Open the PR. The full battery runs on it, takes ~20 minutes, and is a
required check.

### 2. Run the external checks — on the release branch, before merging

win-builder and R-hub run against the release branch, and their results
go into `cran-comments.md` **in the same PR**, so the release commit
carries its own evidence.

``` sh
R CMD build .                                # at the repo root, never in a worktree
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-devel/
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-release/
```

Then trigger R-hub from the Actions tab against the release branch. Two
NOTEs are expected locally and explained in `cran-comments.md` — CRAN
incoming feasibility (new submission of an archived package) and future
file timestamps (no network clock in a sandbox).

This ordering works because **`cran-comments.md` is `.Rbuildignore`d and
never reaches the tarball**, so writing the results into it cannot
invalidate the checks those results describe — the tarball is unchanged.
And `DESCRIPTION` on the release branch already carries the clean
version, so the checks see the exact version string that will be
submitted.

Neither check is run casually: win-builder mails the maintainer, and
both put the tarball in front of a third party.

- **win-builder** — `devtools::check_win_devel()` and
  `devtools::check_win_release()` do the same thing as the `curl` lines
  above. Each FTPs the tarball to `win-builder.r-project.org` and mails
  a check log to the maintainer address within the hour; `226` means the
  transfer completed. The server **will not overwrite an existing
  upload**: a second attempt at the same filename fails with a bare
  `Failed FTP upload: 550`, which means “already queued”, not
  “rejected”. Confirm with
  `curl --list-only ftp://win-builder.r-project.org/R-devel/` before
  re-uploading, or you will conclude the upload failed when it
  succeeded. The mail carries a result URL;
  `curl -s https://win-builder.r-project.org/<key>/00check.log` fetches
  the full log, which beats transcribing the one-line summary. The pages
  expire after ~72 hours, so `cran-comments.md` is the durable copy.

- **R-hub** — `.github/workflows/rhub.yaml` is the stock R-hub v2
  workflow, so the check runs on this repository’s own Actions rather
  than uploading anywhere. Trigger it from the Actions tab (“R-hub” →
  Run workflow), which needs nothing but repository access.
  `rhub::rhub_check()` does the same over the API but requires a GitHub
  PAT with the `repo` scope, stored via `gitcreds::gitcreds_set()` —
  `rhub::rhub_doctor()` reports whether yours has it. Because the
  workflow is `workflow_dispatch`-only, it must be merged to the
  **default branch** before it can be triggered at all, and R-hub wants
  the file on the default branch *and* on whichever branch is being
  checked.

  To read the results, **download the artifacts — do not use
  `gh run view --log`**, which truncates long jobs: at 1.2.1 the
  38-minute macOS job returned 1548 lines ending four minutes in, with
  no `Status:` line anywhere. `gh run download <run-id> -D <dir>` and
  read `*/rcicr.Rcheck/00check.log`. Expect `Status: OK` where
  win-builder reports 1 NOTE — R-hub does not run the CRAN incoming
  feasibility check, so that is agreement, not a discrepancy.

  **`RHUB_TOKEN` is deliberately unset, and nothing is missing.** The
  stock workflow passes `secrets.RHUB_TOKEN` to four actions, which
  makes it look like a prerequisite. It is an optional slot for your own
  PAT to reach private repositories — not a credential R-hub issues —
  and as of `r-hub/actions@v1` none of those four actions references the
  input in any step. This package is public; an unset secret expands to
  an empty string and the jobs run normally.

### 3. Merge and tag

``` sh
gh pr merge <n> --squash --delete-branch
git checkout main && git pull && git fetch --prune
git tag -a vX.Y.Z -m "rcicr X.Y.Z" && git push origin vX.Y.Z
gh release create vX.Y.Z --title "rcicr X.Y.Z" --notes-file <(...)   # NEWS section
```

The tag push re-runs the full gate against the released tree, which is
the copy that matters. Squash-merge, always: one commit per PR on
`main`, and the squash message is where measurements and rejected
alternatives have to live, because the branch commits disappear.

**Confirm the squash did not change the tree** the external checks ran
on:

``` sh
git diff <pr-head-sha> vX.Y.Z --stat    # must be empty
```

Empty output is what lets step 2’s results stand for the tagged tree.

### 4. Submit to CRAN

``` sh
git switch --detach vX.Y.Z      # never build from main HEAD, never in a worktree
R CMD build .
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
git switch -                    # back to where you were
```

Building from the tag is what keeps the development suffix out of the
submitted version: `Version contains large components` is only a CRAN
blocker if the *tarball* carries it.

Submit at <https://cran.r-project.org/submit.html>, pasting the body of
`cran-comments.md` into the “Optional comment” field — that file is
`.Rbuildignore`d and reaches CRAN only this way, never inside the
tarball. `devtools::submit_cran()` does the same thing.

**Paste the tag’s copy, not `main`’s** —
`git show vX.Y.Z:cran-comments.md`. Submission can trail tagging by
weeks, and `main` moves in between; its copy would describe check
results from a tree that is not the tarball, and nothing in the
submission would contradict it.

**Ron submits to CRAN personally.** CRAN emails the maintainer address
to confirm, and for a package archived over an undeliverable address,
that confirmation *is* the point of the resubmission.

### 5. Reopen development

Bump `DESCRIPTION` to `X.Y.Z.9000` and start a fresh
`# rcicr (development version)` heading in `NEWS.md`. `main` is
protected, so this goes through a small PR like anything else.

`.github/workflows/reproducibility.yaml` runs the gate on every PR to
`main` (`--quick`, which skips the 512px configuration) and in full on
release PRs, version tags and manual dispatch. The `compare` job is a
**required status check** on `main`, alongside `R CMD check` on release
and devel, so a red gate blocks the merge rather than merely reporting.
That is a repository ruleset, not something the workflow can enforce on
itself.

## Larger changes

Known work lives in [GitHub
Issues](https://github.com/rdotsch/rcicr/issues), prioritized with the
`P0`–`P3` labels: P0 correctness and availability, P1 dependencies and
toolchain, P2 usability and maintainability, P3 user-requested features.
Each carries the evidence behind it. Read the tracker before starting
anything substantial, and comment on the issue before you begin so the
work is not duplicated.

`DECISIONS.md` is the companion: it records decisions already taken,
**including things that look like bugs and are intentional and must not
be “fixed”**. Check it before changing something that looks wrong.

### Plan first, in the same pull request

**When a change touches behaviour, numbers or a contract, the first
commit on the branch is a plan, and it gets reviewed before any of the
change is written.** That means anything altering `R/` behaviour,
numeric output, the `.Rdata` contract, test fixtures, or the release and
CI machinery. It does *not* mean prose, `man/`, `NEWS.md` wording or
comment-only edits — roughly the inert set
`.github/workflows/reproducibility.yaml` already allowlists.

1.  Branch from `main` and commit the plan to `notes/plan-<topic>.md`.
    `notes/` is already `.Rbuildignore`d via `^notes$` and on that same
    inert allowlist, so a plan-only diff needs no `.Rbuildignore` entry
    and the gate reports green without doing work.
2.  **Open the PR as a draft** and request the review, exactly as in
    “The Codex review” above. The same two conditions clear it.
3.  Implement on that same branch, and **delete the plan file there** as
    part of the work.
4.  **Mark the draft ready**, which re-triggers the review on the full
    diff — a push alone does not.
5.  Squash as usual.

Because the plan file is added and deleted within the branch, the squash
gives `main` one commit — the change itself, carrying no plan file —
while the plan and both review rounds stay on the PR thread. The plan
does not outlive the work and cannot become a second, drifting source of
status; `NEWS.md`, `DECISIONS.md` and the tracker hold what survives.

A plan is worth reviewing only if it can be wrong. State what you
verified and how, name the step most likely to fail, and where the
change rests on a claim about behaviour, measure it rather than
asserting it — an issue’s proposed fix is a hypothesis until it has been
run.

## Code of conduct

Be decent to each other. Problems can be raised privately with the
maintainer at the address in `DESCRIPTION`.
