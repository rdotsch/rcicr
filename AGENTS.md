# AGENTS.md

This file provides guidance to AI coding agents when working with code in this repository.
It is the single source of truth for them; put conventions here.

**Do not delete `CLAUDE.md`** — it is a stub that `@`-imports this file, and Claude Code
loads only `CLAUDE.md`. Keep this file under 200 lines; adherence drops above that.

## What this is

`rcicr` is an R package (CRAN-style) implementing the **reverse correlation image classification** technique from psychophysics: generating noise-based stimuli for 2-image-forced-choice (2IFC) perceptual tasks and computing "classification images" (CIs) from participant response data to visualize internal mental representations (e.g., of faces).

## Common commands

A standard R package — roxygen2 docs, a testthat suite under `tests/testthat/`, GitHub Actions CI — so the usual `roxygen2::roxygenise()` / `devtools::load_all()` / `test()` / `check()` / `install()` workflow applies unchanged, run from the package root. Two things about it are *not* standard:

- `generateStimuli2IFC()` spawns parallel workers via `parallel::makeCluster()` that each `library(rcicr)` themselves (`.packages='rcicr'` in its `foreach` call) — any test or script that calls it needs the package actually **installed** (`devtools::install()`), not just `load_all()`-loaded, or workers fail with "there is no package called 'rcicr'".
- Version and dependency metadata lives in `DESCRIPTION`; user-facing changes go in `NEWS.md`. `ChangeLog` predates `NEWS.md` and takes only a dated pointer entry per release — see `CONTRIBUTING.md`.

## Testing and CI

- `tests/testthat/` has unit tests for every pure/deterministic function (`deg2rad`, `generateSinusoid`, `generateGabor`, `generateNoisePattern`, `generateNoiseImage`, `generateCINoise`), light/targeted tests for the I/O-heavy ones, and an end-to-end smoke test (`test-smoke-pipeline.R`).
- `tests/testthat/helper-fixtures.R` provides shared fixtures: `make_square_png()` (synthetic base face — never a real photo), `make_fixture_rdata()` (runs a tiny `generateStimuli2IFC()` and returns the `.Rdata` path), `seed_reference_norms()` (pre-seeds a `reference_norms` vector so `computeInfoVal2IFC()` skips the expensive/interactive reference-distribution path).
- `tests/testthat/test-fixed-bugs.R` holds regression tests for the P0 bugs in `BACKLOG.md`. They assert *intended* behaviour — never the buggy output they replaced. If one fails, that is a regression; do **not** "fix" it by asserting the broken result, which locks the bug back in.
- `tests/testthat/test-regression-baseline.R` is a golden master pinning the numeric output of the default pipeline — noise basis, classification image, scaling, infoVal — to the values produced *before* the P0 fixes. **If a change turns it red, that change alters researchers' results** and must be documented in `NEWS.md` under "Reproducibility impact" before merging. It is not a test to casually update.
- `tools/compare-release-output.R` is the **release gate**: it installs a released version from its own commit (default `v1.0.1`, tagged retroactively at `b6ab269`) into a temporary library, runs it and the working tree through `tools/compare-harness.R`, and compares every output. It answers what the golden master cannot — that test pins values *this repo computed for itself*, whereas the gate runs the actual old code. A difference is allowed only with an `EXPECTED` entry in the script **and** a matching `NEWS.md` "Reproducibility impact" entry; both are checked, and a stale `EXPECTED` entry fails too. Use `--quick` (~2 min, skips 512px) while iterating. Full checklist in `CONTRIBUTING.md` → "Releasing".
  - The battery is chosen by the **reference** version, not the current one (`RCICR_COMPARE_REF_VERSION`): calls that used to crash — `mask`, z-maps below 512px, undecorated z-maps — have no old value to compare against. Before adding a configuration, check the reference version can execute it; a crash on the reference side aborts the whole gate.
  - It needs ~1.5 GB of RAM at 512px and the reference version's own dependencies — `--install-deps` puts them in a throwaway library rather than the user's.
- Both `devtools::test()` and `testthat::test_local()` set `NOT_CRAN=true` themselves, so neither can be used to verify that a `skip_on_cran()` actually fires.
- `.github/workflows/reproducibility.yaml` runs the gate: `--quick` on every PR to `main`, full on release PRs, `v*` tags and manual dispatch. It tells them apart from `DESCRIPTION` — a version without the `.9000` suffix *is* a release — so the version bump turns the full gate on by itself. On a PR touching only inert files (prose, `man/`, `vignettes/`, `tests/`) the job runs, reports green and does no work; the allowlist is in the workflow, and anything unrecognised runs the gate. **Never convert this to a `paths:` filter** — a skipped required check never reports at all, which GitHub reads as pending forever, making every docs-only PR unmergeable.
- **Five required status checks on `main`**: `compare`, `ubuntu-latest (release)`, `ubuntu-latest (devel)`, `macos-latest (release)`, `windows-latest (release)`. They are enforced by a **ruleset**, not classic branch protection, so `gh api repos/rdotsch/rcicr/branches/main` reports no required contexts and looks unconfigured. Query `gh api repos/rdotsch/rcicr/rules/branches/main` instead.
- `.github/workflows/R-CMD-check.yaml` runs `R CMD check` over four jobs: `ubuntu-latest` on R `release` and `devel`, plus `macos-latest` and `windows-latest` on `release`. macOS and Windows are the two platforms CRAN gates on, and each has caught a failure that was green on Linux for months. **The job names are required checks matched by name** — add rows to the matrix freely, but never rename one, or the renamed check reads as pending forever and blocks every PR.
- Because the checks are required, an infrastructure flake blocks a merge. **This environment cannot re-run one** — `gh run rerun` returns `Resource not accessible by integration`, as does dispatching R-hub — so re-runs and workflow dispatches need the maintainer and the Actions tab.
- **The workflows only trigger on PRs targeting `main`**, so a stacked PR based on another branch gets pre-commit and nothing else — no `R CMD check`. Retargeting an existing PR does *not* re-fire them; close and reopen it.
- `.github/workflows/test-coverage.yaml` uploads to Codecov (needs a `CODECOV_TOKEN` secret); `codecov.yml` sets lenient thresholds because coverage is deliberately partial. `.pre-commit-config.yaml` drives pre-commit.ci with minimal language-agnostic hooks only — `styler`/`lintr`/`roxygenize` were left out because `styler` would reformat nearly every file and destroy `git blame`.
- **Any new top-level file must be added to `.Rbuildignore`** unless it genuinely belongs in the built package, or `R CMD check` NOTEs "non-standard file/directory found at top level" (and a separate NOTE for dotfiles). Add it in the same commit as the file. It is a set of *regexes anchored with `^`*, not globs (e.g. `^BACKLOG\.md$`); read the file for the current list.
  - `^\.git$` is **not** redundant with `R CMD build`'s own exclusion — see [`DECISIONS.md`](DECISIONS.md#git-is-in-rbuildignore-because-a-worktrees-git-is-a-file). Build the release tarball at the repo root, not in a `git worktree`.
  - `.Rbuildignore` and `.gitignore` are **not** interchangeable. `R CMD build` works from the *working directory*, not from git, so a file that is git-ignored but present on disk still ships in the tarball unless `.Rbuildignore` also excludes it. Untracking a file is never a substitute for `.Rbuildignore`ing it.

## DECISIONS.md

`DECISIONS.md` records **why** the package is the way it is: the measurement that ruled an
option out, the alternative that looked obvious and was wrong, the thing that looks like a
bug and is deliberate. It is tracked in git and organised **by theme, not by date**. It
replaced a chronological session log; do not recreate one.

- **Update it as decisions are made**, not in a sweep at the end — the reasoning and the
  numbers are only cheap to write down while they are still to hand.
- **Add an entry when a decision was not obvious.** Rejecting a plausible alternative,
  measuring something surprising, or deliberately *not* fixing something all qualify. Routine
  changes do not.
- **Edit entries in place** when they stop being true. It is not an append-only log, and it
  carries no dates, state or next-steps sections — `BACKLOG.md` holds what is left to do and
  `NEWS.md` holds what changed for users.
- **Do not restate in `DECISIONS.md` what this file, `CONTRIBUTING.md` or `BACKLOG.md`
  already says.** Reference it instead. A duplicated rule drifts, and the copy a reader
  happens to hit first is then wrong.

## Git and merge strategy

- **Merge pull requests to `main` with squash merges** (`gh pr merge <n> --squash`). One
  commit per PR keeps history readable and makes `git revert` of a whole change
  straightforward, which matters here because a PR is usually one self-contained fix plus its
  test and its `NEWS.md` entry. This is repo convention, not enforced by GitHub settings — it
  is on whoever merges to pick the right one.
- The squash commit message is the place to preserve *why* a change was made — the per-commit
  detail on the branch disappears, so measurements, rejected alternatives and reproducibility
  impact belong in the squash message or in `NEWS.md`, not only in the branch commits.
- Delete merged branches. `--delete-branch` on `gh pr merge` handles it; `git fetch --prune`
  clears stale remote refs locally.

## Releases and versioning

Trunk-based with tags — the standard R-package layout. There is **no `develop` branch** and
there should not be one: CRAN has no concept of it, this is a single-maintainer package, and
it would add a permanent second merge direction for no gain. Feature branches → PR → squash
onto `main`, and releases are marked by tags.

- **The reproducibility gate is a release blocker.** `CONTRIBUTING.md` → "Releasing" holds the
  full checklist. The v1.0.1 reference is **pinned and does not advance** with each release —
  see [`DECISIONS.md`](DECISIONS.md#the-v101-reference-is-pinned-the-previous-release-is-a-second-run-not-a-replacement).
- **`main` carries a `.9000` development version between releases.** Right after a release,
  `DESCRIPTION` goes to `<released>.9000`; the release commit drops it to the clean number.
  `NEWS.md` accumulates entries under a `# rcicr (development version)` heading, renamed to
  `# rcicr X.Y.Z (date)` at release time.
- **Tag every release** — `git tag -a vX.Y.Z <release commit>` plus a GitHub release.
- **Log every CRAN reply verbatim** in `notes/cran-review-<version>.md`, named for the
  version whose tarball it reviews — `cran-review-1.2.1.md` reviews the 1.2.1 submission,
  answered by 1.2.2 and 1.2.3. Add a new file per reply rather than editing an old one, and
  **answer from the file, not from a summary of it**: a summary dropped one of the two
  filenames in a point, and three drafts of the reply then told the reviewer we could not
  find what she had named.
- **Build the CRAN tarball from the tag, never from `main` HEAD.** This is also what keeps the
  `.9000` suffix safe: `Version contains large components` is only a CRAN blocker if the
  *submitted tarball* carries it, and a tarball built from a tag never does. Do not drop the
  development-version convention to avoid that NOTE.
- **Never put a version number in a `NEWS.md` section heading.** `R CMD check` parses the file
  to build the news database, and a `##` heading containing something version-shaped makes it
  treat `##` as the *version* level — after which every other section title fails to yield a
  version and the whole file NOTEs. Name the version in the body text instead.
- **Write claims that survive on someone else's machine.** A bare wall-clock number measured
  here ("the example set runs in about nine seconds") can be contradicted by the reader's own
  log, and has been. Give a ratio, a comparison to a fixed bar (CRAN's five-second
  per-example limit), or just "faster". Absolute times are fine where the ratio is the point
  — "about 6x faster, 1.66s to 0.28s per call" — and where a user will feel the difference;
  they are noise when they only describe our hardware. Same rule for `cran-comments.md`,
  where the reviewer *has* their own log.
- **Order `NEWS.md` entries largest-impact first** within each section — changes to numeric
  output or return values, then behaviour changes, then bug fixes that only ever produced
  errors, then message-only fixes. Someone who stops reading after three bullets should have
  read the three that could change their results.

## Backlog

`BACKLOG.md` is the prioritized modernization backlog — read it before starting substantial work, and update it when you close something out.

**It carries what is left to do, never where the project currently stands.** No "state as of" block, no release narrative, no next-steps list — that belongs to `NEWS.md`, `cran-comments.md`, `CONTRIBUTING.md` → Releasing, and the open PRs. A hold condition belongs on the item it holds; the project's position does not.

It also records the **guiding constraint**: researchers re-run old analysis scripts years later, so never change existing call syntax, argument meanings, or the numeric output of a function silently. Deprecate rather than delete; treat the `.Rdata` contract as append-only. And it records what is **already correct and should not be "re-fixed"** — notably the infoVal formula, which already matches the published Schmitz et al. erratum.

## Architecture

The package has two pipeline halves that share state only through an `.RData` file written at
stimulus-generation time. The per-function walkthrough, the data-flow diagram and the anatomy
of the `.RData` file live in `README.md` — sections "How it works" and "Anatomy of the
`.Rdata` file". Read them there rather than maintaining two copies that drift.

### Notes on conventions in this codebase

- Functions use `save_as_png=TRUE` / `save_rdata=TRUE`-style side-effecting defaults — most analysis functions write PNGs to disk in addition to returning data structures. **The destination is always a required argument** (`stimulus_path`, `targetpath`, `zmaptargetpath`): none has a default, because a default path writes to the user's filespace uninvited and CRAN policy forbids it. Never reintroduce one — not even `tempdir()`; [`DECISIONS.md`](DECISIONS.md#write-paths-are-required-arguments-not-defaults-of-tempdir) records why.
- Scaling of CI pixel intensities (`none`, `constant`, `matched`, `independent`) is a key user-facing decision, documented at length in `generateCI.R`'s roxygen header — read it before changing scaling logic.
- `computeInfoVal2IFC()`'s `ref_lookup` tibble looks like a cache and is not one: **it has been empty since 2018**, its rows having been measured under the pre-erratum infoVal formula. Every lookup misses and the reference distribution is always regenerated. Do not describe it as a working cache; the matching machinery is kept only so the table can be repopulated cheaply.
- `pre_0.3.0` / `generator_version` fields keep backward compatibility with `.Rdata` files from older versions (index-counter starts at 0 vs 1) — do not remove without understanding this.
- `R/zzz.R` declares `utils::globalVariables()` for names that only exist at runtime after `load()`-ing an `.Rdata` file (e.g. `base_faces`, `stimuli_params`, `p`, `seed`) or inside `foreach` loops (`obs`) — add new ones here to avoid `R CMD check` NOTEs.
- **`load()` assigns into the calling function's frame**, so a field in the `.Rdata` can silently replace a function argument of the same name. Every `load()` site keeps copies of its arguments across the call; preserve that when adding either an argument or an `.Rdata` field.
- Parallelism is via base `parallel` + `doParallel`/`foreach`, not newer alternatives (e.g. `future`) — match this pattern for new parallel code, and remember worker processes need `.packages='rcicr'` set on `foreach` calls.
