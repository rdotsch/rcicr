# Maintaining the repository

How this repository's automation is wired, and why; `AGENTS.md` says which document owns what.

**Keep this file under 1800 words.**

---

## The workflows

| workflow | what it does |
|---|---|
| `R-CMD-check.yaml` | `R CMD check` on `ubuntu-latest` (release and devel), `macos-latest` and `windows-latest`. macOS and Windows are the platforms CRAN gates on, and each has caught a failure that had been green on Linux for months. Also runs weekly and on dispatch, so a break caused from outside the repository shows up even when nobody is pushing; the comment on its `schedule:` lists the cron's own failure modes. |
| `reproducibility.yaml` | The release gate: `--quick` on every PR to `main`; the full battery on release PRs, `v*` tags and dispatch. It tells them apart by `DESCRIPTION`: a version without `.9000` *is* a release, so the version bump turns on the full gate by itself. |
| `test-coverage.yaml` | Codecov; needs a `CODECOV_TOKEN`. `codecov.yml` sets lenient thresholds because coverage is deliberately partial. |
| `pkgdown.yaml` | Builds the site on every PR and deploys it on push to `main`. |
| `rhub.yaml` | Stock R-hub v2, `workflow_dispatch` only. |
| `lint.yaml` | `lintr::lint_package()`, which must come back clean; `.lintr` holds the config and an empty `exclusions:`. A required check (below), with `lintr` pinned; `tools/regenerate-lintr-baseline.R` rebuilds `exclusions:` if that is ever needed. |

**The required status checks on `main`** are `compare`, `ubuntu-latest (release)`, `ubuntu-latest (devel)`, `macos-latest (release)`, `windows-latest (release)` and `lint`. For the live list, query `gh api repos/rdotsch/rcicr/rules/branches/main`. A **ruleset** enforces them, not classic branch protection, so `gh api repos/rdotsch/rcicr/branches/main` reports no required contexts and looks unconfigured. Agents here cannot add a name to that ruleset, so a new workflow's check reports on PRs but blocks none until the maintainer adds it under GitHub → Settings → Rules → Rulesets. Use the web interface: a malformed `PUT` risks dropping unrelated ruleset fields. A separate **Tags** ruleset blocks creating, moving and deleting tags; the maintainer (Repository admin role) bypasses it to tag releases. A ruleset also refuses force-pushes to `claude/*` branches, so catch a stale branch up by merging `main` into it, not by rebasing.

Two rules for editing a workflow:

- **Never rename a job.** Required checks are matched by name. Rows can be added to a matrix freely, but a renamed check stays "pending" forever and blocks every PR.
- **Never turn the gate into a `paths:` filter.** A skipped required check never reports, which GitHub also reads as pending forever, so every docs-only PR would become unmergeable. Instead the job always runs and exits early when every changed path is on an allowlist of inert paths, kept in the workflow. Anything not on the list runs the gate.

`.pre-commit-config.yaml` drives pre-commit.ci with minimal, language-agnostic hooks only. `lintr` runs as `lint.yaml` rather than as a hook: `object_usage_linter` needs the package *installed*, not just parsed, which a pre-commit.ci hook would have to do on every commit. **`styler` is a tool, not a hook.** As a hook it needs an R toolchain, and `tidyverse_style()` rewrites single quotes to double, undoing the choice `.lintr` records by disabling `quotes_linter`. Scope it to `I(c("spaces", "indention"))`, in a commit of its own; #255 tracks the rest.

**Keep `exclusions:` empty; fix the lint instead.** Never pin an exclusion to a line number: it does not survive an edit above it, and seven once shifted at the same time, failing CI far below the edit. For a false positive, put `# nolint: <linter>.` on the line. `tools/regenerate-lintr-baseline.R` refuses to baseline the linters that mean something is genuinely wrong (`object_usage`, `commented_code`, `object_name`, `seq`).

---

## CI, checks and packaging

### `fail_ci_if_error: false` on the coverage workflow
Chosen over getting a Codecov token or deleting the workflow. **Coverage reporting still does not work**: without a token there is no badge, no per-PR comment, and `codecov.yml`'s thresholds do nothing. What the setting buys is narrower: a red `main` now means the package is broken. The `token:` input stays wired, so adding the secret later restores reporting.

### `^\.git$` is in `.Rbuildignore` because a worktree's `.git` is a *file*
`R CMD build` drops `.git` by itself, so the entry looks redundant. It is not. The built-in exclusion matches a **directory** named `.git`, but in a `git worktree` checkout `.git` is a 49-byte text file holding a `gitdir:` pointer, and it ships. The tarball then gets `checking for hidden files and directories ... NOTE` naming `.git`, from a tree that looks identical to a clean one.

This reached win-builder. Building the tag in a worktree (the obvious way to follow "never build from `main` HEAD") gave 2 NOTEs where the same commit built at the repository root gave 1. The worktree tarball contains `rcicr/.git`; the root tarball does not.

`.Rbuildignore` itself was listed in `.gitignore` from 2016 until 2026, a leftover from the RStudio template, so it never reached CI or other contributors. Do not add it back.

### R-hub runs on `workflow_dispatch` only, never on push
This is the stock file `rhub::rhub_setup()` writes, left unmodified so it can be refreshed from upstream. It answers a release-time question (does the package build under CRAN's own compiler flags, on platforms nobody develops on), so running it on every PR would spend a platform matrix on a question nobody asked, alongside the gate that *is* required.

**R-hub does not replace per-platform checks.** The everyday matrix once varied only the R version on `ubuntu-latest`: one platform, twice. The gap stayed invisible until an R-hub run failed a `plotZmap()` test on macOS that had been green on Linux for months. CRAN builds on macOS, so it would have become a submission ERROR. **"Already covered by X" is a claim about what X runs; read X to check it.** Both platforms CRAN gates on are now in the matrix.

The cost of dispatch-only: GitHub offers no "Run workflow" button for a file on a feature branch, so the workflow must be merged to the **default branch** before the release it checks.

### The pkgdown site deploys via a `gh-pages` branch, not the Actions-native route
GitHub Pages was already configured to serve a `gh-pages` branch (`gh api repos/rdotsch/rcicr/pages` reports `build_type: legacy`); the branch just had never been pushed, which is why the URL returned 404. Publishing that branch from a workflow therefore needed **no change to repository settings**. The Actions-native deploy is more modern, but it requires changing `build_type`, a settings write agents here cannot make (`Resource not accessible by integration`), which would stall the PR on the maintainer.

`docs/` holds the built site and is **both git-ignored and `.Rbuildignore`d**. The two are not interchangeable, and a locally built site has been committed by accident once.

### The stale-`man/` gate is a step in an existing job, and a pre-commit hook was rejected
`R CMD check` already fails when documentation disagrees with a function's *signature* (`tools::codoc()`). It cannot see prose drift: an edited `@description` or `@examples` that was never regenerated. That drift still reaches the pkgdown site, which builds from `man/`. So the gate re-runs roxygen and diffs `man/` and `NAMESPACE`.

- **A step in `ubuntu-latest (release)`, not a new job.** A new job's check would report without ever blocking, for the ruleset reason above. It runs last, so a failure cannot hide the check results.
- **Not a pre-commit hook.** The check job already has R and every dependency installed, where a hook would install them on every run, and an optional local hook is not a gate.
- **roxygen2 is pinned to `RoxygenNote`, not `any::`.** `.Rd` output varies between roxygen2 versions, so `any::` would let a new roxygen2 release turn a required check red on formatting alone. The step asserts that the installed version matches `DESCRIPTION`, so bumping one means bumping both.

### `CITATION.cff` is generated by `cffr` and compared, not hand-written and field-checked
A hand-written file with a field-by-field comparison was **rejected**: every defect found in it was a field the comparison did not cover. Unordered URL matching let `repository-code` and `url` swap; comparing authors by name left the email *address* unchecked, on a package CRAN archived for an undeliverable address. Such a comparison only covers the fields someone remembered, and it fails silently.

`cffr::cff_create()` derives the whole file, so the comparison covers everything. It also validates against the CFF schema; GitHub silently hides the "Cite this repository" button for a file it cannot parse.

Three generation settings matter:

- `dependencies = FALSE`. With `TRUE` the file gains 380 lines of dependency authors, with years taken from whichever versions are installed, so each machine produces a different file.
- `gh_keywords = FALSE`, which keeps generation off the network.
- The comparison excludes `preferred-citation`'s year. It comes from `inst/CITATION`, which reads the clock when `DESCRIPTION` has no publication date, so it changes every 1 January and would otherwise turn a required check red on a calendar boundary.

---

## Repository documents

### A work tracker carries no project state
`BACKLOG.md`'s "state as of" block drifted four times, each within a day of a release. That is what a duplicated fact does, not what a careless author does: the fourth fix was to write it more carefully, and it was wrong within the hour. The rule now applies to the issue tracker. A *hold condition* belongs on the issue it holds; the project's current position belongs nowhere in it. The same reasoning removed the hand-maintained `Version:`/`Date:` table from `man/rcicr-package.Rd`.

The tracker is a **working surface, not a curated public one**, so internal maintenance work sits next to user-visible bugs. Moving chores elsewhere was rejected: two backlogs is the same duplication one level up. What does *not* become an issue is an "already correct, do not re-fix" decision about the package; that goes in `DECISIONS.md`.

### `CLAUDE.md` is a stub that imports `AGENTS.md`, not a symlink
Claude Code reads `CLAUDE.md`, not `AGENTS.md`, so renaming the file in #166 silently stopped the conventions loading in every agent session for months. **A symlink was rejected**: on Windows it needs Administrator rights or Developer Mode, and this package has Windows contributors and a Windows CI runner. The `@AGENTS.md` import is a plain file that every platform treats the same.
