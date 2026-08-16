# Contributing to rcicr

Contributions, thoughts and criticisms are welcome. This file records
the few conventions that are specific to this package — the rest is
ordinary R package practice. Releases are in `RELEASING.md`, the
repository’s automation in `MAINTENANCE.md`, and why the package behaves
as it does in `DECISIONS.md`.

**Keep this file under 2800 words.** A convention nobody reaches is not
a convention; over budget, something comes out before something goes in.

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
compares. Every release has to pass it — see `RELEASING.md`.

## Getting set up

On a fresh Ubuntu machine with only R installed (no compiler, no package
library) – `tools/dev-setup.sh` is Ubuntu-specific, it shells out to
`apt-get` – it installs the system toolchain and every Imports/Suggests
package from source, then installs the package itself. Idempotent –
rerun it after pulling a dependency change.

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
that were entirely the sandbox; installing it was the only change needed
to reach 2. **Do not reach for `--no-manual` to make the manual checks
go away**: it skips them rather than passing them, and a note here once
recorded that as the problem being *resolved*. A clean run shows
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

Codex reviews pull requests here and has caught real errors. Nothing in
the merge path makes you notice it: it submits as `COMMENTED`, so
`gh pr checks` stays green and `gh pr view --comments` shows only the
wrapper, never the findings.

**It must not become something that can block.** If it is switched off,
erroring, or simply not answering, merge on the other checks.

Push everything first — a push does not re-trigger the review; only this
comment, opening the PR, or marking a draft ready does. Keep the
timestamp it returns:

``` sh
trig=$(gh api repos/rdotsch/rcicr/issues/<n>/comments -f body='@codex review' --jq '.created_at')
```

Two conditions clear a squash.

1.  **A 👍 newer than `$trig`.** The reaction tracks the latest run: 👀
    running, 👍 clean, none when it has findings. Filter on the numeric
    account id — reactions are open to anyone with read access, and an
    id survives a rename.

    ``` sh
    gh api repos/rdotsch/rcicr/issues/<n>/reactions \
      --jq '.[] | select(.user.id == 199175422) | "\(.content) \(.created_at)"'
    ```

    Anything else — 👀, nothing, an older `+1` — means read the findings
    and answer each. `--paginate` is not optional; that endpoint pages
    at 30 and truncates silently.

    ``` sh
    gh api --paginate repos/rdotsch/rcicr/pulls/<n>/comments --jq '.[] | "\(.path): \(.body)"'
    ```

2.  **Every answered thread resolved**, via the `resolveReviewThread`
    mutation. This half is **enforced rather than checked**: the `main`
    ruleset sets `required_review_thread_resolution`, so an unresolved
    thread blocks the squash server-side.

Do not compute “is this safe to merge” client-side. Earlier versions of
this section derived it from review objects and `commit_id` and were
wrong six separate ways, each failing *open*. The reaction is the one
signal Codex sets deliberately; thread resolution is GitHub’s to
enforce.

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
| Roxygen | exported functions only | `man/` holds exactly the 17 exports plus the package doc; internal helpers use plain `#` comments, including the ones sharing a file with an exported function. Roxygen on a non-exported function either publishes a help page for something no user can call, or needs `@noRd` — a comment in stricter syntax. There are none. |
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
