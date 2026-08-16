# Maintaining the repository

How this repository's automation is wired, and why. `CONTRIBUTING.md` is how to contribute,
`RELEASING.md` how to cut a release, `DECISIONS.md` why the package behaves as it does — this
file is the machinery around them.

**Keep this file under 1800 words.**

---

## The workflows

| workflow | what it does |
|---|---|
| `R-CMD-check.yaml` | `R CMD check` on `ubuntu-latest` (release + devel), `macos-latest`, `windows-latest`. macOS and Windows are the platforms CRAN gates on, and each has caught a failure green on Linux for months. |
| `reproducibility.yaml` | the release gate: `--quick` on every PR to `main`, full on release PRs, `v*` tags and dispatch. It tells them apart from `DESCRIPTION` — no `.9000` *is* a release — so the version bump turns the full gate on by itself. |
| `test-coverage.yaml` | Codecov; needs a `CODECOV_TOKEN`. `codecov.yml` sets lenient thresholds because coverage is deliberately partial. |
| `pkgdown.yaml` | builds the site on every PR, deploys on push to `main`. |
| `rhub.yaml` | stock R-hub v2, `workflow_dispatch` only. |
| `lint.yaml` | `lintr::lint_package()`, gated by `.lintr`'s `exclusions:` baseline so only a *new* lint fails it. A required check (below); `lintr` is pinned in the workflow, and `tools/regenerate-lintr-baseline.R` regenerates the baseline. |

**The required status checks on `main`** are `compare`, `ubuntu-latest (release)`,
`ubuntu-latest (devel)`, `macos-latest (release)`, `windows-latest (release)` and `lint` — query
`gh api repos/rdotsch/rcicr/rules/branches/main` for the live list rather than trust a count
here, since it drifts by one every time a check is added or removed. They are enforced by a
**ruleset**, not classic branch protection, so `gh api repos/rdotsch/rcicr/branches/main` reports
no required contexts and looks unconfigured. `lint` was added last, once the
maintainer added its check *name* to the ruleset by hand via GitHub → Settings → Rules →
Rulesets: adding a name to the ruleset is a repo-settings write agents here cannot make, so a
newly added workflow's check reports on PRs but does not block one until the maintainer does
this — a UI edit rather than a hand-built API payload, since a malformed `PUT` risks dropping
unrelated ruleset fields.

Two consequences worth knowing before editing a workflow. **Required checks are matched by
name**, so rows can be added to a matrix freely but never renamed — a renamed check reads as
pending forever and blocks every PR. And **the gate must never become a `paths:` filter**: a
skipped required check never reports at all, which GitHub also reads as pending forever, making
every docs-only PR unmergeable. Instead the job always runs and exits early on an allowlist of
inert paths, which is in the workflow; anything unrecognised runs the gate.

`.pre-commit-config.yaml` drives pre-commit.ci with minimal language-agnostic hooks only.
`styler` is deliberately left out: it would reformat nearly every file in one sweep and destroy
`git blame`. `lintr` runs as `lint.yaml` instead of a pre-commit hook — `object_usage_linter`
needs the package *installed*, not just parsed, which the language-agnostic hooks here never
need and a pre-commit.ci hook would have to carry on every commit.

**Exclusions are never pinned to line numbers** — they do not survive an edit above them, and
once shifted seven at once, failing CI far below. Fix a lint where you can: a clean
file leaves `exclusions:` entirely. Six dense legacy files keep
`<linter> = Inf`; a later cosmetic lint there goes unnoticed, the right trade for code
deliberately not reformatted. Where a hit means something is *wrong* — `object_usage`,
`commented_code`, `object_name`, `object_length`, `seq` — the line carries
`# nolint: <linter>.`, which `tools/regenerate-lintr-baseline.R` refuses to baseline.

---

## CI, checks and packaging

### `fail_ci_if_error: false` on the coverage workflow
Chosen over getting a Codecov token or deleting the workflow. **This does not make coverage
reporting work** — no token means no badge, no per-PR comments, and `codecov.yml`'s thresholds
stay inert. What it buys is narrower and more useful: a red `main` now means the package is
broken. The `token:` input is left wired so adding the secret later restores reporting.

### `^\.git$` is in `.Rbuildignore` because a worktree's `.git` is a *file*
`R CMD build` drops `.git` on its own, so the entry looks redundant. It is not: the built-in
exclusion matches a **directory** named `.git`, and in a `git worktree` checkout `.git` is a
49-byte text file holding a `gitdir:` pointer. It ships. The tarball then draws
`checking for hidden files and directories ... NOTE`, naming `.git`, from a tree that looks
identical to a clean one.

This is not hypothetical — it reached win-builder. Building the tag in a worktree (the
obvious way to honour "never build from `main` HEAD") returned 2 NOTEs where the same commit
built at the repo root returned 1. Reproduced both ways before believing it: worktree tarball
contains `rcicr/.git`, repo-root tarball contains no match.

Separately: `.Rbuildignore` was itself listed in `.gitignore` from 2016 until 2026, an
RStudio-template leftover, so it never reached CI or other contributors. Do not re-add it.

### R-hub runs on `workflow_dispatch` only, never on push
The stock file `rhub::rhub_setup()` writes, kept unmodified so it can be refreshed upstream. It
answers a release-time question — does this build under CRAN's own compiler flags, on platforms
nobody develops on — so running it per-PR would spend a platform matrix on a question nobody
asked, competing with the reproducibility gate that *is* required.

**R-hub does not stand in for per-platform coverage.** The everyday matrix once varied only the
R version against a pinned `ubuntu-latest` — one platform twice — and the hole stayed invisible
until an R-hub dispatch failed a `plotZmap()` test on macOS that had been green on Linux for
months. CRAN builds on macOS, so it would have landed as a submission ERROR. **"Already covered
by X" is a claim about what X runs, worth reading X to check.** Both platforms CRAN gates on are
in the matrix now.

Cost of dispatch-only: the workflow is invisible until it reaches the **default branch** —
there is no "Run workflow" button for a file on a feature branch — so it merges ahead of the
release it checks.

### The pkgdown site deploys via a `gh-pages` branch, not the Actions-native route
Pages was already configured here — `gh api repos/rdotsch/rcicr/pages` reported `status:
built`, `build_type: legacy`, source branch `gh-pages` — and the only thing missing was that
no such branch had ever been pushed, which is why the URL 404'd. Publishing that branch from
a workflow therefore needed **no repository setting to change**. The Actions-native deploy
would have been the more modern choice and was rejected on one fact: it requires flipping
`build_type`, a repo-settings write that agents here cannot make (`Resource not accessible by
integration`), turning a self-contained PR into one that stalls on the maintainer.

`docs/` is the built site and is **git-ignored as well as `.Rbuildignore`d** — the two are not
interchangeable, and a locally built site has already been committed by accident once.

### The stale-`man/` gate is a step in an existing job, and a pre-commit hook was rejected
`R CMD check` already fails on documentation that disagrees with a function's *signature* —
`tools::codoc()`, confirmed by renaming `deg2rad`'s formal without regenerating `man/`. What it
cannot see is prose drift: an edited `@description` or `@examples` never regenerated. The same
edit leaves `codoc()` silent and still reaches the pkgdown site, which builds from `man/`. So
the gate re-runs roxygen and diffs `man/` and `NAMESPACE`.

Three things about where it lives:

- **A step in `ubuntu-latest (release)`, not a new job.** Required checks are matched by name,
  and adding a name to the ruleset is a repo-settings write agents here cannot make — a new job
  would report without ever blocking. It runs last so a failure cannot cost us the check results.
- **Not a pre-commit hook.** The check job already has R and every dependency installed where a
  hook would install them per run, and an optional local hook is not a gate.
- **roxygen2 is pinned to `RoxygenNote`, not `any::`.** `.Rd` output varies between generator
  versions, so `any::` lets a CRAN release of roxygen2 redden a required check on formatting
  alone — an external event breaking the gate, in an environment that cannot re-run a job. The
  step asserts the installed version matches `DESCRIPTION`, so bumping one means bumping both.

### `CITATION.cff` is generated by `cffr` and compared, not hand-written and field-checked
A hand-written file with a field-by-field comparator was built first and **rejected**: review
found five defects in four rounds, every one the same shape — a field the comparison did not
reach. Unordered URL membership let `repository-code` and `url` swap; comparing authors by name
left the *address* unchecked, on a package CRAN archived for an undeliverable address. A
hand-written comparator is a list of the fields someone remembered, and it fails silently.

`cffr::cff_create()` derives the whole file, so the comparison is structural and nothing depends
on remembering a field. It also validates against the CFF schema — GitHub silently declines to
render the "Cite this repository" button on a file it cannot parse.

Two generation settings are load-bearing. `dependencies = TRUE` emits 380 lines of dependency
authors with years read from whichever versions happen to be installed, so a different machine
produces a different file: off. And `preferred-citation`'s year comes from `inst/CITATION`, which
reads the clock when `DESCRIPTION` has no publication date — that field changes on 1 January
with nothing else, so the comparison excludes it or a required check reddens on a calendar
boundary. `gh_keywords = FALSE` keeps generation off the network.

---


---

## Repository documents

### A work tracker carries no project state
`BACKLOG.md`'s "state as of" block drifted four times, each within a day of a release — what a
duplicated fact does, not what a careless author does. The fourth fix was to write it more
carefully, and it was wrong within the hour. Applies now to the issue tracker: a *hold
condition* belongs on the issue it holds, the project's current position nowhere in it. Same
reasoning as deleting the hand-maintained `Version:`/`Date:` table from `man/rcicr-package.Rd`.

The tracker is a **working surface, not a curated public one**, so internal maintenance work
sits in it alongside user-visible bugs. Splitting chores elsewhere was rejected: two backlogs is
the same duplication one layer up. What does *not* become an issue is an "already correct — do
not re-fix" decision about the package; those go in `DECISIONS.md`.

### `CLAUDE.md` is a stub that imports `AGENTS.md`, not a symlink
Claude Code reads `CLAUDE.md` and not `AGENTS.md`, so the rename in #166 silently stopped the
conventions loading in every agent session for months. **The symlink form was rejected**: it
needs Administrator privileges or Developer Mode on Windows, and this package has Windows
contributors and a Windows CI runner. The `@AGENTS.md` import is a plain file every platform
treats identically.

---
