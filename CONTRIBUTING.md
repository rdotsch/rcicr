# Contributing to rcicr

Contributions, thoughts and criticisms are welcome. This file holds the conventions specific to this package; everything else is ordinary R package practice. Releases are covered in `RELEASING.md`, the repository's automation in `MAINTENANCE.md`, and why the package behaves as it does in `DECISIONS.md`.

**Keep this file under 2800 words.** Over budget, something comes out before something goes in.

## The one constraint that shapes everything else

**Researchers re-run old analysis scripts years later, and publish the numbers that come out.** A change that silently alters what `generateCI()` or `computeInfoVal2IFC()` returns can invalidate a published result without anyone noticing. So:

- Do not change existing call syntax, argument meanings or numeric output silently. Deprecate rather than delete.
- Treat the `.Rdata` file's contents as **append-only**: add fields, never rename them or change their meaning. The README documents its layout.
- A change that *does* alter numeric output is not automatically wrong, but it must be deliberate and described in `NEWS.md` under a "Reproducibility impact" heading: who is affected, and what they should do.

`tests/testthat/test-regression-baseline.R` is a golden master that pins the default pipeline's output. **If your change turns it red, your change alters researchers' results.** Document it; do not update the baseline to match.

`tools/compare-release-output.R` asks the same question from outside. It installs a released version of the package from its own commit, runs both versions over a battery of configurations, and compares the results. Every release has to pass it; see `RELEASING.md`.

## Getting set up

`tools/dev-setup.sh` sets up a fresh Ubuntu machine that has only R. It installs the system toolchain with `apt-get`, then every Imports and Suggests package from source, then the package itself. It is safe to re-run after a dependency change.

```r
devtools::install()      # not load_all(); see below
devtools::test()
devtools::check()
roxygen2::roxygenise()   # after editing any roxygen comment; commit man/ and NAMESPACE
```

`man/` and `NAMESPACE` are generated and tracked. `ubuntu-latest (release)` regenerates them and fails on any difference. Use the roxygen2 version named in `DESCRIPTION`'s `RoxygenNote`; the CI step is pinned to it and says so when they disagree.

`generateStimuli2IFC()` starts parallel workers that each call `library(rcicr)`, so anything that reaches it needs the package **actually installed**. Under `devtools::load_all()` alone the workers fail with `there is no package called 'rcicr'`.

**An `--as-cran` run needs a complete check environment, or it reports failures the package does not have.** Install `texlive` and `tidy`, and set both incoming variables:

```sh
_R_CHECK_CRAN_INCOMING_=TRUE _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE \
  R CMD check --as-cran rcicr_X.Y.Z.tar.gz
```

Without that toolchain, a run reported 1 ERROR, 1 WARNING and 4 NOTEs that all came from the environment; installing it brought the result down to 2 NOTEs. **Do not use `--no-manual` to make the manual checks go away**: it skips them rather than passing them. A clean run shows `checking PDF version of manual ... OK` and its HTML equivalent; without those two lines, the manual has not been checked.

## Reporting a bug

The most useful report gives the `.Rdata` file's `img_size`, `nscales` and `noise_type`, the rcicr version, and a reproducible snippet. A synthetic base image works fine as a stand-in and is better than attaching a real face photo: `png::writePNG(matrix(runif(32*32), 32, 32), "base.png")`.

## Pull requests

- **Branch from `main`**, one self-contained change per PR: the fix, its test and its `NEWS.md` entry together. There is no `develop` branch.
- **Tests are expected.** Every pure function has unit tests; the I/O-heavy ones have lighter, targeted tests. A bug fix comes with a test that fails without it. Check that with `git stash push -- R/`; a plain `git stash` also stashes your new tests, so the suite passes without testing anything.
- **A test's title is a claim; make the assertions support it.** Asserting that a result has the right shape is not asserting that it is correct. Where it is cheap, also assert that the *wrong* answer differs, so the test cannot pass vacuously.
- **Order `NEWS.md` entries largest impact first** within each section: changes to numeric output or return values, then behaviour changes, then fixes for bugs that only ever produced errors, then message-only fixes. Someone who stops after three bullets should have read the three that could change their results.
- **Add every new top-level file to `.Rbuildignore`** unless it belongs in the built package, or `R CMD check` NOTEs "non-standard file/directory found at top level". The file holds `^`-anchored *regular expressions*, not globs.
- **Run `git diff --stat main...HEAD` before opening the PR.** `R CMD check` leaves a full copy of the package behind, and one has been committed by accident.
- PRs are **squash-merged** to `main`, so put what a future reader needs (measurements, rejected alternatives, reproducibility impact) in the PR description or `NEWS.md`, not only in branch commits.
- **Write for the state the change ends in, not the route you took.** Commit messages, PR descriptions and code comments say what the change is and why it is right. How many attempts it took, and what each got wrong, disappears in the squash and cannot be seen in the merged code. Keep a rejected alternative only when someone would otherwise try it again, and give it one line.

### Testing

- **`skip_if_not_installed()` is for `Suggests` packages only.** `withr` may genuinely be absent, so skipping on it is honest. A package in `Imports` (`png`, `jpeg`, `matlab`) cannot be absent, because the package will not load without it. Such a skip never fires for a real reason, and if it somehow did, it would hide the test instead of failing it. At a glance a skip reads as "passed", which is the opposite of what you need to know.
- **When mutating the code repeatedly, keep backups with `cp` in a scratch directory, and never restore with `git checkout <file>`.** That discards unstaged work, and it has destroyed an in-progress implementation here. `git stash push -- R/` is right for proving a single fix; it is the restore that bites. Guard each mutation with a `grep -q MUTANT` check too: a mutation that silently failed to apply looks like a surviving mutant.
- **When a test reads pixels back from a graphics device, assert only relationships between renders.** Every absolute property of those pixels belongs to the device: the channel count (cairo writes RGB, macOS quartz RGBA) and the values (quartz renders a 0.5 background at about 0.573 where cairo gives 0.502). Render onto a *uniform* background, so "drew nothing" becomes "the image is one flat value"; count distinct values over the colour channels only; and compare two renders instead of pinning a number. [`DECISIONS.md`](DECISIONS.md#pixel-assertions-have-measured-the-graphics-device-twice) has the two failures behind this rule.
- **Check figures by looking at them.** Three problems in the walkthrough vignette passed every assertion and were obvious on sight; [`DECISIONS.md`](DECISIONS.md#three-vignette-figures-were-wrong-in-ways-only-viewing-them-showed) lists them.
- **Base images in tests and vignettes are always synthetic**, never a real photograph. Avoiding licensing and consent questions entirely is better than a realistic-looking figure.

### The Codex review

Codex reviews pull requests here and has caught real errors. Nothing in the merge path makes you notice it: it submits as `COMMENTED`, so `gh pr checks` stays green, and `gh pr view --comments` shows only the wrapper, never the findings.

**It must never become something that blocks.** If it is switched off, erroring or not answering, merge on the other checks.

Push everything first; a push never re-triggers the review. **Marking a draft ready cannot be relied on either**: once it triggered a review within four minutes, and once not at all in thirty-five. Post the request yourself and keep its timestamp:

```sh
trig=$(gh api repos/rdotsch/rcicr/issues/<n>/comments -f body='@codex review' --jq '.created_at')
```

Two conditions clear a squash:

1. **A 👍 newer than `$trig`.** Codex's reaction tracks its latest run: 👀 while running, 👍 when clean, none when it has findings. Filter on the numeric account id: anyone with read access can add reactions, and an id survives a rename.
   ```sh
   gh api repos/rdotsch/rcicr/issues/<n>/reactions \
     --jq '.[] | select(.user.id == 199175422) | "\(.content) \(.created_at)"'
   ```
   Anything else (👀, nothing, an older `+1`) means: read the findings and answer each one. Use `--paginate`: the endpoint pages at 30 and truncates silently.
   ```sh
   gh api --paginate repos/rdotsch/rcicr/pulls/<n>/comments --jq '.[] | "\(.path): \(.body)"'
   ```
2. **Every answered thread resolved**, with the `resolveReviewThread` mutation. GitHub **enforces** this half: the `main` ruleset sets `required_review_thread_resolution`, so an unresolved thread blocks the squash.

Do not work out "is this safe to merge" on the client side. Earlier attempts derived it from review objects and `commit_id` and were wrong in six ways, each of which let a merge through. The reaction is the one signal Codex sets on purpose, and thread resolution is GitHub's to enforce.

## Code conventions

These describe what the code already does. **They apply to new and modified code.** Nobody is expected to reformat untouched files, and no sweep is planned; issue #194 describes the one that would be, and why it needs a major version.

**Two conventions are frozen by the constraint above, not chosen:**

- **Exported functions are camelCase**: `generateCI()`, `batchGenerateCI2IFC()`, `plotZmap()`. A modern R style guide would say otherwise, but these names are in researchers' stored scripts and in published methods sections.
- **Arguments and list fields are snake_case**: `img_size`, `n_trials`, `noise_type`, `save_as_png`. Also frozen, since people call them by name. Arguments therefore already match tidyverse style; only function names differ.

The rest is ordinary consistency. Internal helpers are free to change, because nothing outside the package can call them.

| | Convention | Notes |
|---|---|---|
| Assignment | `<-`, never `=` | Already universal in `R/`. |
| Internal helpers | camelCase, like the exported functions | All but `default_ncores()` follow this; new helpers follow the majority. |
| Strings | single quotes | Mixed today; not worth the churn to unify, but write new code with single quotes. |
| Indentation | 2 spaces, no tabs | Already consistent. |
| Booleans | `TRUE`/`FALSE`, never `T`/`F` | `T` and `F` can be reassigned. Inside the package they resolve through the namespace, so this is style rather than a hazard, but it costs nothing and settles the question. |
| Sequences | `seq_len()`/`seq_along()`, not `1:n` | `1:0` counts *backwards*, so `1:length(x)` on an empty vector runs twice. |
| Returns | explicit `return()` at the end of exported functions | The existing style throughout. |
| Files | one file per exported function, named after it | `R/generateCI.R` holds `generateCI()`; `zzz.R` holds the `globalVariables()` declarations. |
| Roxygen | exported functions only | `man/` holds exactly the exports plus the package page. Internal helpers use plain `#` comments, even when they share a file with an export: roxygen on an unexported function either publishes a page no user can reach or needs `@noRd`. |
| Namespacing | `pkg::fn()` or `@importFrom pkg fn`, not `@import pkg` | `@import matlab` once masked `base::sum()` with MATLAB semantics across six files (#182); the package now calls `matlab::` explicitly. |

Line length is not enforced, and some lines in `R/` exceed 100 characters. Wrap new code at something reasonable instead of reflowing what is there.

**Do not line-wrap Markdown prose: write one line per paragraph and let the reader's editor wrap it. Avoid em dashes**; use a full stop, comma, colon or parentheses instead.

**Show why something matters; never just say that it does.** Cut "It is worth …", "it matters because", "in ways that matter" and similar, and state the consequence. Not "worth running first" but "a `TRUE` ends the matter".

**Comments: as few as will do, and only where the code cannot speak for itself.** Write one when the reason would otherwise have to be worked out again: a rejected alternative, a constraint, a measured number, a trap. Never narrate what the next line does, repeat the error message below it, or retell an issue's history. The long comments already in `R/` are of the first kind and earn their length; the test is whether a reader would otherwise get it wrong, not whether the code looks bare. An explanation about the package rather than the line belongs in [`DECISIONS.md`](DECISIONS.md).

**If the package is ever run through `styler`, that is a commit of its own**, never a side effect of other work, so `git blame` stays useful. `MAINTENANCE.md` explains why `styler` is not a pre-commit hook.

Three rules about this package rather than about R:

- **Do not write a package-qualified call as code (in backticks or `\code{}`) for a package the docs only *mention*.** The pkgdown site resolves such a link by loading that package, and a package that is installed but cannot load takes the whole site build down. After `raster` was dropped in #186, the CI runners kept it in their cached library but stopped installing its GDAL/PROJ system libraries, and the build died with `libproj.so.25: cannot open shared object file`: once from the roxygen, and again from the `NEWS.md` entry describing the change. Name it in prose instead: "the raster package's plot method". A package that is simply absent is fine. `README.md`, `NEWS.md`, `DECISIONS.md`, this file and the vignettes are all rendered.
- **Add new names loaded from an `.Rdata` file to `globalVariables()` in `R/zzz.R`**, or `R CMD check` NOTEs about undefined globals.
- **Check every new function argument against the names saved in the `.Rdata` file**, and every new saved field against existing argument names. `load()` assigns into the calling frame, so a name collision silently overwrites the argument. This has caused three separate bugs, most recently a saved `sigma` replacing `generateCI()`'s z-map blur `sigma`. [`DECISIONS.md`](DECISIONS.md#load-assigns-into-the-calling-frame--check-every-new-argument-against-saved-names) has the details.

## Larger changes

Known work lives in [GitHub Issues](https://github.com/rdotsch/rcicr/issues), prioritised with the `P0`–`P3` labels: P0 correctness and availability, P1 dependencies and toolchain, P2 usability and maintainability, P3 user-requested features. Each issue carries the evidence behind it. Read the tracker before starting anything substantial, and comment on the issue before you begin so the work is not duplicated.

`DECISIONS.md` records decisions already taken, **including things that look like bugs, are intentional, and must not be "fixed"**. Check it before changing something that looks wrong.

### Plan first, in the same pull request

**When a change touches behaviour, numbers or a contract, the branch's first commit is a plan, reviewed before any of the change is written.** That covers changes to `R/` behaviour, numeric output, the `.Rdata` contract, test fixtures, and the release and CI machinery. It does *not* cover prose, `man/`, `NEWS.md` wording or comment-only edits: roughly the inert set that `.github/workflows/reproducibility.yaml` already allowlists.

1. Branch from `main` and commit the plan as `notes/plan-<topic>.md`. `notes/` is already `.Rbuildignore`d (`^notes$`) and on the inert allowlist, so a plan-only diff needs no `.Rbuildignore` entry and the gate passes without running.
2. **Open the PR as a draft** and request the review as described in "The Codex review" above. The same two conditions clear it.
3. Implement on the same branch, and **delete the plan file there** as part of the work.
4. **Mark the draft ready**, then post the `@codex review` comment.
5. Squash as usual.

Because the plan file is added and deleted on the branch, the squash leaves `main` one commit, the change itself, with no plan file. The plan and both review rounds stay on the PR thread. The plan cannot outlive the work or become a second source of status that drifts; `NEWS.md`, `DECISIONS.md` and the tracker hold what lasts.

A plan is worth reviewing only if it can be wrong. State what you verified and how, and name the step most likely to fail. Where the change rests on a claim about behaviour, measure it: an issue's proposed fix is a hypothesis until it has been run.

## Code of conduct

Be decent to each other. Raise problems privately with the maintainer at the address in `DESCRIPTION`.
