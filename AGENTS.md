# AGENTS.md

The single source of conventions for AI coding agents in this repository.

**Do not delete `CLAUDE.md`.** It is a stub that `@`-imports this file, and Claude Code loads only `CLAUDE.md`. **Keep this file under 2800 words.**

## What this is

`rcicr` is an R package for **reverse correlation image classification**. It generates noise-based stimuli for two-image forced-choice (2IFC) tasks, and computes classification images (CIs) from the responses to visualize mental representations.

## Check, don't assume

**Verify a claim against the thing itself before writing it down or acting on it.**

- **Run the command, read the file, query the API.** Never write "should be" or "as expected"; if you did not run a check, say so instead of predicting its result.
- **An empty result may mean you asked the wrong question**: wrong endpoint, branch, path or scope. Rule that out before reporting that nothing was found.
- **Put the evidence next to the claim**: the command, the file and line, the actual output.
- **Numbers in prose come from the run that made the tables**, comparisons included. Hand-typed figures drift.
- **Claim only what the estimate supports.** A non-significant difference is not an absence, a mean is not a bound, and a ratio of magnitudes is not a decomposition. *Conservative* and *worst case* need evidence beside them. After fixing such a claim, grep for it elsewhere.

## Common commands

This is a standard R package (roxygen2 docs, testthat suite under `tests/testthat/`, GitHub Actions CI), so the usual `roxygen2::roxygenise()` / `devtools::load_all()` / `test()` / `check()` / `install()` workflow applies, run from the package root. Two things are *not* standard:

- `generateStimuli2IFC()` starts parallel workers via `parallel::makeCluster()`, and each worker runs `library(rcicr)` (`.packages='rcicr'` in its `foreach` call). Any test or script that calls it needs the package **installed** (`devtools::install()`); with only `load_all()`, workers fail with "there is no package called 'rcicr'".
- Version and dependency metadata lives in `DESCRIPTION`; user-facing changes go in `NEWS.md`. `ChangeLog` is the frozen archive from before `NEWS.md`, ending at 1.0.1. Never add to it.

## Testing and CI

**`MAINTENANCE.md` lists the workflows and the required check names, and gives the two rules for editing them: never rename a job, never turn the gate into a `paths:` filter.** Read it before touching `.github/workflows/`.

- `tests/testthat/` has unit tests for every pure function (`deg2rad`, `generateSinusoid`, `generateGabor`, `generateNoisePattern`, `generateNoiseImage`, `generateCINoise`), lighter targeted tests for the I/O-heavy ones, and an end-to-end smoke test (`test-smoke-pipeline.R`).
- `tests/testthat/helper-fixtures.R` provides shared fixtures: `make_square_png()` (a synthetic base face, never a real photo), `make_fixture_rdata()` (runs a tiny `generateStimuli2IFC()` and returns the `.Rdata` path) and `seed_reference_norms()` (pre-seeds `reference_norms` so `computeInfoVal2IFC()` skips the slow reference simulation).
- `tests/testthat/test-fixed-bugs.R` holds regression tests for the P0 bugs fixed in the modernization pass. They assert the *intended* behaviour, never the buggy output they replaced. If one fails, that is a regression: do **not** "fix" it by asserting the broken result, which locks the bug back in.
- `tests/testthat/test-regression-baseline.R` is a golden master pinning the default pipeline's numeric output (noise basis, classification image, scaling, infoVal) to the values from *before* the P0 fixes. **A change that turns it red alters researchers' results** and must be documented in `NEWS.md` under "Reproducibility impact" before merging. Do not casually update it.
- `tools/compare-release-output.R` is the **release gate**. It installs a released version (by default `v1.0.1`, tagged retroactively at `b6ab269`) into a temporary library, runs it and the working tree through `tools/compare-harness.R`, and compares every output. The golden master pins values *this repository computed for itself*; the gate runs the actual old code. Every difference needs both an `EXPECTED` entry in the script and a matching `NEWS.md` "Reproducibility impact" entry, and a stale `EXPECTED` entry fails too. Use `--quick` (about 2 minutes, skips 512px) while iterating. The full checklist is in `RELEASING.md`.
  - The battery is chosen by the **reference** version, not the current one (`RCICR_COMPARE_REF_VERSION`). Calls that used to crash (`mask`, z-maps below 512px, undecorated z-maps) have no old value to compare against. Before adding a configuration, check that the reference version can run it: a crash on the reference side aborts the whole gate.
  - It needs about 1.5 GB of RAM at 512px, plus the reference version's own dependencies; `--install-deps` puts those in a throwaway library rather than the user's.
- `devtools::test()` and `testthat::test_local()` both set `NOT_CRAN=true` themselves, so neither can verify that a `skip_on_cran()` actually fires.
- The checks are required, so an infrastructure flake blocks a merge. **Re-run it from here**: `gh run rerun <id> --failed` works, and so does a workflow dispatch. To test dispatch permission, use a bogus ref: 422 means allowed (bad ref), 403 means no permission.
- **The workflows only trigger on PRs targeting `main`.** A stacked PR based on another branch gets pre-commit and nothing else, so no `R CMD check`. Retargeting an existing PR does *not* re-trigger them; close and reopen it.
- **Add every new top-level file to `.Rbuildignore`** in the same commit, unless it belongs in the built package (`CONTRIBUTING.md` → "Pull requests"). Entries are `^`-anchored regular expressions, not globs (e.g. `^DECISIONS\.md$`).
  - `^\.git$` is **not** redundant with `R CMD build`'s own exclusion; see [`MAINTENANCE.md`](MAINTENANCE.md#git-is-in-rbuildignore-because-a-worktrees-git-is-a-file). Build the release tarball at the repository root, not in a `git worktree`.
  - `.Rbuildignore` and `.gitignore` are **not** interchangeable. `R CMD build` works from the *working directory*, not from git, so a file that is git-ignored but on disk still ships unless `.Rbuildignore` excludes it. Untracking a file never replaces `.Rbuildignore`ing it.

## Which file a thing goes in

Each document has one job and **a word budget**; over budget, something comes out before something goes in. Put each fact in exactly one of them and refer to it from the others. A duplicated rule drifts, and the copy a reader finds first is then wrong.

| file | its job | budget |
|---|---|---|
| `AGENTS.md` | conventions for agents (this file) | 2800 |
| `CONTRIBUTING.md` | how to contribute: setup, tests, PRs, code conventions | 2800 |
| `RELEASING.md` | the release checklist and why it is in that order | 1600 |
| `MAINTENANCE.md` | how the repository's CI, gates and generated files are wired | 1800 |
| `SECURITY.md` | vulnerability reporting and dependency posture | 600 |
| `DECISIONS.md` | why the **package** behaves as it does | 5200 |
| `NEWS.md` | what changed for users | none; trimmed at each release |

Check word counts with `LC_ALL=C.UTF-8 wc -w` (Unicode whitespace counts as a separator).

`DECISIONS.md` is where things are most often misfiled. Its subject is what `generateCI()` returns and why a number cannot change, **not** how CI is wired or how a release is cut. The test: would this still matter if the package were maintained somewhere else entirely?

- **Update it when a decision is made**, not in a sweep at the end, while the reasoning and the numbers are still to hand.
- **Add an entry when a decision was not obvious**: rejecting a plausible alternative, measuring something surprising, or deliberately *not* fixing something. Routine changes do not qualify.
- **Edit entries in place** when they stop being true. The file is organised by theme and has no dates, status or next-steps sections. It replaced a chronological session log; do not recreate one.
- **An edge case in a third-party tool is not a decision.** Work it out again if it recurs.

## Git and merge strategy

- **Re-read `CONTRIBUTING.md` → "Pull requests" every time you open a PR; do not work from memory.** Its steps fail silently rather than with an error, so a half-remembered version looks like it worked. The rule most often lost is one of omission: a PR body or commit message states the end result, never the route taken. A hook in `.claude/settings.json` repeats it at the moment of writing.
- **Plan first when a change touches behaviour, numbers or a contract**: `R/` behaviour, numeric output, the `.Rdata` contract, fixtures, or the release and CI machinery. The plan is the branch's first commit, reviewed as a **draft** PR before any of the change is written, then deleted on the same branch, so the squash leaves `main` one commit and no plan file. The full procedure is in `CONTRIBUTING.md` → "Plan first, in the same pull request". Prose, `man/`, `NEWS.md` wording and comment-only edits are exempt.
- **For a PR's implementation work alongside others in progress, use `subagent_type: "fork"`**, not a generic subagent. A fork inherits the full context; other types start cold and need the PR and prior decisions handed over.
- **Merge pull requests to `main` with squash merges** (`gh pr merge <n> --squash`). One commit per PR keeps history readable and makes `git revert` of a whole change straightforward; a PR here is usually one fix plus its test and its `NEWS.md` entry. The `main` ruleset enforces it (`allowed_merge_methods` is `["squash"]`).
- The branch's individual commits disappear in the squash, so the squash commit message is where *why* is kept. Measurements, rejected alternatives and reproducibility impact go in the squash message or in `NEWS.md`, not only in branch commits.
- **Read the Codex review before squashing, and answer it.** It is easy to merge past: it never blocks (not a required check, submitted as `COMMENTED`), and neither `gh pr checks` nor `gh pr view --comments` shows its findings. Two things clear a squash: a 👍 on the PR body dated after your own `@codex review` comment, and no unresolved review threads. `CONTRIBUTING.md` → "The Codex review" has one command for each.
- Delete merged branches: `--delete-branch` on `gh pr merge` does it, and `git fetch --prune` clears stale remote refs locally.

## Releases and versioning

Trunk-based with tags, the standard R package layout. There is **no `develop` branch**, and there should not be one: CRAN has no concept of it, this is a single-maintainer package, and it would add a permanent second merge direction for no gain. Feature branches go through a PR and are squashed onto `main`; tags mark releases.

- **The reproducibility gate blocks a release.** `RELEASING.md` has the full checklist. The v1.0.1 reference is **pinned and does not advance** with each release; see [`DECISIONS.md`](DECISIONS.md#the-v101-reference-is-pinned-the-previous-release-is-a-second-run-not-a-replacement).
- **`main` carries a `.9000` development version between releases.** Right after a release, `DESCRIPTION` goes to `<released>.9000`; the release commit drops it to the clean number. `NEWS.md` collects entries under `# rcicr (development version)`, renamed to `# rcicr X.Y.Z (date)` at release time.
- **Tag every release** with `git tag -a vX.Y.Z <release commit>` plus a GitHub release. A tag marks a *release*, not CRAN acceptance; `RELEASING.md` → "Reopen development" says why, and where acceptance is recorded instead.
- **Log every CRAN reply verbatim** in `notes/cran-review-<version>.md`, named for the version whose tarball it reviews: `cran-review-1.2.1.md` reviews the 1.2.1 submission, answered by 1.2.2 and 1.2.3. Add a new file per reply instead of editing an old one, and **answer from the file, not from a summary of it**. A summary once dropped one of two filenames in a point, and three drafts of the reply then told the reviewer we could not find what she had named.
- **Build the CRAN tarball from the tag, never from `main` HEAD.** That is also what makes the `.9000` suffix safe: `Version contains large components` only blocks CRAN if the *submitted tarball* carries it, and a tarball built from a tag never does. Do not drop the development-version convention to avoid that NOTE.
- **Never put a version number in a `NEWS.md` section heading.** `R CMD check` parses the file to build the news database. A `##` heading containing something version-shaped makes it treat `##` as the *version* level, after which every other section title fails to yield a version and the whole file NOTEs. Name the version in the body text instead.
- **Write claims that hold on someone else's machine.** A bare wall-clock time ("runs in about nine seconds") can be contradicted by the reader's own log, and has been. Give a ratio, a comparison with a fixed bar (CRAN's five-second limit per example), or just "faster". Absolute times are fine only where the ratio is the point ("about 6x faster, 1.66s to 0.28s"). The same goes for `cran-comments.md`, where the reviewer has their own log.
- **Order `NEWS.md` entries largest impact first** within each section; the full order is in `CONTRIBUTING.md` → "Pull requests".

## The guiding constraint

**Researchers re-run old analysis scripts years later.** Never silently change existing call syntax, argument meanings or a function's numeric output. Deprecate rather than delete, and treat the `.Rdata` contract as append-only. When output does change, it goes in `NEWS.md` under "Reproducibility impact", and the release gate has to agree (see "Testing and CI" above).

Some things **are already correct and must not be "re-fixed"**, notably the infoVal formula, which matches the published Schmitz et al. erratum. [`DECISIONS.md`](DECISIONS.md) records them with the reasoning; read it before changing something that looks wrong.

## Work is tracked in GitHub Issues

Use `gh issue list`; priorities are the `P0`–`P3` labels: P0 correctness and availability, P1 dependencies and toolchain, P2 usability and maintainability, P3 user-requested features. Read the tracker before starting substantial work, and close the issue when you finish.

The tracker replaced `BACKLOG.md`, an in-repository file whose status lived in two hand-maintained places that disagreed six times. **Do not recreate it.** An issue's state *is* its status, so that failure cannot happen. Anything settled, including a rejected alternative or a deliberate non-fix, belongs in `DECISIONS.md`, not in an issue.

## Architecture

The package has two halves that share state only through an `.Rdata` file written when the stimuli are generated. The per-function walkthrough, the data-flow diagram and the anatomy of the `.Rdata` file are in `README.md`, sections "How it works" and "Anatomy of the `.Rdata` file". Read them there; two copies would drift.

### Conventions in this codebase

- **Comment sparingly**: only where the reason would otherwise have to be worked out again, never to narrate the next line. This is the convention most often broken here; the full rule is in `CONTRIBUTING.md` → "Code conventions".
- Most analysis functions write PNGs in addition to returning data, through `save_as_png = TRUE` / `save_rdata = TRUE`-style defaults. **The destination is always a required argument** (`stimulus_path`, `targetpath`, `zmaptargetpath`) with no default: a default path writes into the user's files uninvited, which CRAN policy forbids. Never reintroduce one, not even `tempdir()`; [`DECISIONS.md`](DECISIONS.md#write-paths-are-required-arguments-not-defaults-of-tempdir) records why.
- How CI pixel intensities are scaled (`none`, `constant`, `matched`, `independent`) is a key user-facing choice, documented in the roxygen header of `generateCI.R`. Read it before changing the scaling logic.
- `computeInfoVal2IFC()`'s `ref_lookup` tibble looks like a cache and is not one: **it has been empty since 2018**, because its rows were measured under the pre-erratum infoVal formula. Every lookup misses, and the reference distribution is always simulated. Do not describe it as a working cache; the matching code is kept only so the table can be repopulated cheaply.
- The `pre_0.3.0` and `generator_version` fields keep `.Rdata` files from older versions working (their index counter starts at 0 instead of 1). Do not remove them without understanding this.
- `R/zzz.R` declares `utils::globalVariables()` for names that only exist at run time, after `load()`-ing an `.Rdata` file (e.g. `base_faces`, `stimuli_params`, `p`, `seed`) or inside `foreach` loops (`obs`). Add new ones there to avoid `R CMD check` NOTEs.
- **`load()` assigns into the calling function's frame**, so a field in the `.Rdata` file can silently replace a function argument of the same name. Every `load()` site keeps copies of its arguments across the call; preserve that when adding an argument or an `.Rdata` field.
- Parallelism uses base `parallel` plus `doSNOW`/`foreach`, not newer alternatives such as `future`. Match that pattern in new parallel code, and remember that workers need `.packages='rcicr'` on `foreach` calls. Tick progress bars from the parent with `.options.snow = progressOption(pb, cl)`; a tick inside the loop body only updates a worker's private copy.
