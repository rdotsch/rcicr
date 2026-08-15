# AGENTS.md

This file provides guidance to AI coding agents when working with code
in this repository. It is the single source of truth for them; put
conventions here.

**Do not delete `CLAUDE.md`** — it is a stub that `@`-imports this file,
and Claude Code loads only `CLAUDE.md`. **Keep this file under 2800
words** — counted in words, not lines, because a line budget is defeated
by longer lines and this file has been the worst offender.

## What this is

`rcicr` is an R package (CRAN-style) implementing the **reverse
correlation image classification** technique from psychophysics:
generating noise-based stimuli for 2-image-forced-choice (2IFC)
perceptual tasks and computing “classification images” (CIs) from
participant response data to visualize internal mental representations
(e.g., of faces).

## Check, don’t assume

**Verify a claim against the thing itself before writing it down or
acting on it.** Several rules below are scar tissue from not doing that
— the branch-protection API, answering CRAN from a summary, quoting this
machine’s clock — and each keeps its detail where someone doing that
task will hit it. The general form:

- **Run the command, read the file, query the API.** No “should be” or
  “as expected”; if a check was not run, say so rather than predicting
  its result.
- **An empty result may mean you asked the wrong question** — wrong
  endpoint, branch, path or scope. Rule that out before reporting
  nothing found.
- **Put the evidence next to the claim**: the command, the file and
  line, the actual output.

## Common commands

A standard R package — roxygen2 docs, a testthat suite under
`tests/testthat/`, GitHub Actions CI — so the usual
[`roxygen2::roxygenise()`](https://roxygen2.r-lib.org/reference/roxygenize.html)
/ `devtools::load_all()` / `test()` / `check()` / `install()` workflow
applies unchanged, run from the package root. Two things about it are
*not* standard:

- [`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
  spawns parallel workers via
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
  that each [`library(rcicr)`](https://rdotsch.github.io/rcicr/)
  themselves (`.packages='rcicr'` in its `foreach` call) — any test or
  script that calls it needs the package actually **installed**
  (`devtools::install()`), not just `load_all()`-loaded, or workers fail
  with “there is no package called ‘rcicr’”.
- Version and dependency metadata lives in `DESCRIPTION`; user-facing
  changes go in `NEWS.md`. `ChangeLog` predates `NEWS.md` and takes only
  a dated pointer entry per release — see `CONTRIBUTING.md`.

## Testing and CI

**The workflow inventory, the required check names, and the two rules
for editing them (never rename a job, never convert the gate to a
`paths:` filter) are in `MAINTENANCE.md`.** Read it before touching
`.github/workflows/`.

- `tests/testthat/` has unit tests for every pure/deterministic function
  (`deg2rad`, `generateSinusoid`, `generateGabor`,
  `generateNoisePattern`, `generateNoiseImage`, `generateCINoise`),
  light/targeted tests for the I/O-heavy ones, and an end-to-end smoke
  test (`test-smoke-pipeline.R`).
- `tests/testthat/helper-fixtures.R` provides shared fixtures:
  `make_square_png()` (synthetic base face — never a real photo),
  `make_fixture_rdata()` (runs a tiny
  [`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
  and returns the `.Rdata` path), `seed_reference_norms()` (pre-seeds a
  `reference_norms` vector so
  [`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
  skips the expensive/interactive reference-distribution path).
- `tests/testthat/test-fixed-bugs.R` holds regression tests for the P0
  bugs found in the modernization pass. They assert *intended* behaviour
  — never the buggy output they replaced. If one fails, that is a
  regression; do **not** “fix” it by asserting the broken result, which
  locks the bug back in.
- `tests/testthat/test-regression-baseline.R` is a golden master pinning
  the numeric output of the default pipeline — noise basis,
  classification image, scaling, infoVal — to the values produced
  *before* the P0 fixes. **If a change turns it red, that change alters
  researchers’ results** and must be documented in `NEWS.md` under
  “Reproducibility impact” before merging. It is not a test to casually
  update.
- `tools/compare-release-output.R` is the **release gate**: it installs
  a released version (default `v1.0.1`, tagged retroactively at
  `b6ab269`) into a temporary library, runs it and the working tree
  through `tools/compare-harness.R`, and compares every output. It
  answers what the golden master cannot — that test pins values *this
  repo computed for itself*, whereas the gate runs the actual old code.
  A difference needs both an `EXPECTED` entry in the script and a
  matching `NEWS.md` “Reproducibility impact” entry — a stale `EXPECTED`
  entry fails too. Use `--quick` (~2 min, skips 512px) while iterating.
  Full checklist in `RELEASING.md`.
  - The battery is chosen by the **reference** version, not the current
    one (`RCICR_COMPARE_REF_VERSION`): calls that used to crash —
    `mask`, z-maps below 512px, undecorated z-maps — have no old value
    to compare against. Before adding a configuration, check the
    reference version can execute it; a crash on the reference side
    aborts the whole gate.
  - It needs ~1.5 GB of RAM at 512px and the reference version’s own
    dependencies — `--install-deps` puts them in a throwaway library
    rather than the user’s.
- Both `devtools::test()` and
  [`testthat::test_local()`](https://testthat.r-lib.org/reference/test_package.html)
  set `NOT_CRAN=true` themselves, so neither can verify a
  `skip_on_cran()` actually fires.
- Because the checks are required, an infrastructure flake blocks a
  merge. **This environment cannot re-run one** — `gh run rerun` returns
  `Resource not accessible by integration`, as does dispatching R-hub —
  so re-runs and workflow dispatches need the maintainer and the Actions
  tab.
- **The workflows only trigger on PRs targeting `main`**, so a stacked
  PR based on another branch gets pre-commit and nothing else — no
  `R CMD check`. Retargeting an existing PR does *not* re-fire them;
  close and reopen it.
- **Any new top-level file must be added to `.Rbuildignore`** unless it
  genuinely belongs in the built package, or `R CMD check` NOTEs
  “non-standard file/directory found at top level” (and a separate NOTE
  for dotfiles). Add it in the same commit as the file. It is a set of
  *regexes anchored with `^`*, not globs (e.g. `^DECISIONS\.md$`); read
  the file for the current list.
  - `^\.git$` is **not** redundant with `R CMD build`’s own exclusion —
    see
    [`MAINTENANCE.md`](https://rdotsch.github.io/rcicr/MAINTENANCE.html#git-is-in-rbuildignore-because-a-worktrees-git-is-a-file).
    Build the release tarball at the repo root, not in a `git worktree`.
  - `.Rbuildignore` and `.gitignore` are **not** interchangeable.
    `R CMD build` works from the *working directory*, not from git, so a
    file that is git-ignored but present on disk still ships in the
    tarball unless `.Rbuildignore` also excludes it. Untracking a file
    is never a substitute for `.Rbuildignore`ing it.

## Which file a thing goes in

Each doc has one job, and **each states a word budget** — over it,
something comes out before something goes in. Put a fact in exactly one
of them and reference it from the others; a duplicated rule drifts, and
the copy a reader hits first is then wrong.

| file | its job | budget |
|----|----|----|
| `AGENTS.md` | conventions for agents (this file) | 2800 |
| `CONTRIBUTING.md` | how to contribute: setup, tests, PRs, code conventions | 2800 |
| `RELEASING.md` | the release checklist and why it is ordered that way | 1600 |
| `MAINTENANCE.md` | how the repo’s CI, gates and generated files are wired | 1800 |
| `SECURITY.md` | vulnerability reporting and dependency posture | 600 |
| `DECISIONS.md` | why the **package** behaves as it does | 5200 |
| `NEWS.md` | what changed for users | none — trimmed at each release |

Budgets are in **words, not lines**: a line budget is defeated by
writing longer lines, which is exactly how this file grew to hold more
words than `CONTRIBUTING.md` in two-thirds the lines. `wc -w` is the
check.

`DECISIONS.md` is the one most often misfiled into. Its subject is what
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
returns and why a number cannot change — **not** how CI is wired or how
a release is cut. The test: would this still matter if the package were
maintained somewhere else entirely?

- **Update it as decisions are made**, not in a sweep at the end — the
  reasoning and the numbers are only cheap to write down while they are
  still to hand.
- **Add an entry when a decision was not obvious.** Rejecting a
  plausible alternative, measuring something surprising, or deliberately
  *not* fixing something all qualify. Routine changes do not.
- **Edit entries in place** when they stop being true. It is organised
  by theme, carries no dates, state or next-steps sections, and replaced
  a chronological session log — do not recreate one.
- **An edge case in a third-party tool is not a decision.** Work it out
  again if it recurs.

## Git and merge strategy

- **Re-read `CONTRIBUTING.md` → “Pull requests” when opening a PR, every
  time — not from memory.** Its steps fail by producing silence rather
  than an error, so a half-remembered version looks like it worked. Two
  of them have now been re-learned the hard way, after being written
  down.
- **Plan first when a change touches behaviour, numbers or a contract**
  — `R/` behaviour, numeric output, the `.Rdata` contract, fixtures, or
  the release and CI machinery. The plan is the branch’s first commit,
  reviewed as a **draft** PR before any of the change is written, then
  deleted on that same branch so the squash leaves `main` one commit and
  no plan file. Full procedure in `CONTRIBUTING.md` → “Plan first, in
  the same pull request”; read it there. Prose, `man/`, `NEWS.md`
  wording and comment-only edits are exempt.
- **Use `subagent_type: "fork"` (not a generic subagent) for a PR’s
  implementation work alongside others in flight.** A fork inherits full
  context for free; other types start cold and need the PR and prior
  decisions handed over.
- **Merge pull requests to `main` with squash merges**
  (`gh pr merge <n> --squash`). One commit per PR keeps history readable
  and makes `git revert` of a whole change straightforward, which
  matters here because a PR is usually one self-contained fix plus its
  test and its `NEWS.md` entry. The `main` ruleset enforces it —
  `allowed_merge_methods` is `["squash"]`, so the other buttons are not
  offered.
- The squash commit message is the place to preserve *why* a change was
  made — the per-commit detail on the branch disappears, so
  measurements, rejected alternatives and reproducibility impact belong
  in the squash message or in `NEWS.md`, not only in the branch commits.
- **Read the Codex review before squashing, and answer it.** It is easy
  to merge past: it never blocks — not a required check, submits as
  `COMMENTED` — and neither `gh pr checks` nor `gh pr view --comments`
  shows its findings. Two things clear a squash: a 👍 on the PR body
  dated after your own `@codex review` comment, and no unresolved review
  threads. One command each, in `CONTRIBUTING.md` → “The Codex review”.
- Delete merged branches. `--delete-branch` on `gh pr merge` handles it;
  `git fetch --prune` clears stale remote refs locally.

## Releases and versioning

Trunk-based with tags — the standard R-package layout. There is **no
`develop` branch** and there should not be one: CRAN has no concept of
it, this is a single-maintainer package, and it would add a permanent
second merge direction for no gain. Feature branches → PR → squash onto
`main`, and releases are marked by tags.

- **The reproducibility gate is a release blocker.** `RELEASING.md`
  holds the full checklist. The v1.0.1 reference is **pinned and does
  not advance** with each release — see
  [`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.html#the-v101-reference-is-pinned-the-previous-release-is-a-second-run-not-a-replacement).
- **`main` carries a `.9000` development version between releases.**
  Right after a release, `DESCRIPTION` goes to `<released>.9000`; the
  release commit drops it to the clean number. `NEWS.md` accumulates
  entries under a `# rcicr (development version)` heading, renamed to
  `# rcicr X.Y.Z (date)` at release time.
- **Tag every release** — `git tag -a vX.Y.Z <release commit>` plus a
  GitHub release.
- **Log every CRAN reply verbatim** in `notes/cran-review-<version>.md`,
  named for the version whose tarball it reviews —
  `cran-review-1.2.1.md` reviews the 1.2.1 submission, answered by 1.2.2
  and 1.2.3. Add a new file per reply rather than editing an old one,
  and **answer from the file, not from a summary of it**: a summary
  dropped one of the two filenames in a point, and three drafts of the
  reply then told the reviewer we could not find what she had named.
- **Build the CRAN tarball from the tag, never from `main` HEAD.** This
  is also what keeps the `.9000` suffix safe:
  `Version contains large components` is only a CRAN blocker if the
  *submitted tarball* carries it, and a tarball built from a tag never
  does. Do not drop the development-version convention to avoid that
  NOTE.
- **Never put a version number in a `NEWS.md` section heading.**
  `R CMD check` parses the file to build the news database, and a `##`
  heading containing something version-shaped makes it treat `##` as the
  *version* level — after which every other section title fails to yield
  a version and the whole file NOTEs. Name the version in the body text
  instead.
- **Write claims that survive on someone else’s machine.** A bare
  wall-clock number (“runs in about nine seconds”) can be contradicted
  by the reader’s own log, and has been. Give a ratio, a comparison to a
  fixed bar (CRAN’s five-second per-example limit), or just “faster” —
  absolute times are fine only where the ratio itself is the point
  (“about 6x faster, 1.66s to 0.28s”). Same rule for `cran-comments.md`,
  where the reviewer has their own log.
- **Order `NEWS.md` entries largest-impact first** within each section —
  changes to numeric output or return values, then behaviour changes,
  then bug fixes that only ever produced errors, then message-only
  fixes. Someone who stops reading after three bullets should have read
  the three that could change their results.

## The guiding constraint

**Researchers re-run old analysis scripts years later.** Never change
existing call syntax, argument meanings, or the numeric output of a
function silently. Deprecate rather than delete; treat the `.Rdata`
contract as append-only. When output does change it goes in `NEWS.md`
under “Reproducibility impact”, and the release gate has to agree — see
“Testing and CI” above.

Some things **are already correct and must not be “re-fixed”** — notably
the infoVal formula, which matches the published Schmitz et al. erratum.
[`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.md) records
those with their reasoning; read it before changing something that looks
wrong.

## Work is tracked in GitHub Issues

`gh issue list`, prioritized with the `P0`–`P3` labels: P0 correctness
and availability, P1 dependencies and toolchain, P2 usability and
maintainability, P3 user-requested features. Read the tracker before
starting substantial work, and close the issue when you finish.

This replaced `BACKLOG.md`, an in-repo file whose status lived in two
hand-maintained places that disagreed six times. **Do not recreate it**
— an issue’s state *is* its status, so that failure mode disappears
rather than being managed. Anything settled, including a rejected
alternative or a deliberate non-fix, belongs in `DECISIONS.md` rather
than in an issue.

## Architecture

The package has two pipeline halves that share state only through an
`.RData` file written at stimulus-generation time. The per-function
walkthrough, the data-flow diagram and the anatomy of the `.RData` file
live in `README.md` — sections “How it works” and “Anatomy of the
`.Rdata` file”. Read them there rather than maintaining two copies that
drift.

### Notes on conventions in this codebase

- **Comment sparingly.** Only where the reason would otherwise have to
  be re-derived, never to narrate the next line. This is the convention
  most often broken here — full rule in `CONTRIBUTING.md` → “Code
  conventions”.
- Functions use `save_as_png=TRUE` / `save_rdata=TRUE`-style
  side-effecting defaults — most analysis functions write PNGs to disk
  in addition to returning data structures. **The destination is always
  a required argument** (`stimulus_path`, `targetpath`,
  `zmaptargetpath`): none has a default, because a default path writes
  to the user’s filespace uninvited and CRAN policy forbids it. Never
  reintroduce one — not even
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html);
  [`DECISIONS.md`](https://rdotsch.github.io/rcicr/DECISIONS.html#write-paths-are-required-arguments-not-defaults-of-tempdir)
  records why.
- Scaling of CI pixel intensities (`none`, `constant`, `matched`,
  `independent`) is a key user-facing decision, documented at length in
  `generateCI.R`’s roxygen header — read it before changing scaling
  logic.
- [`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)’s
  `ref_lookup` tibble looks like a cache and is not one: **it has been
  empty since 2018**, its rows having been measured under the
  pre-erratum infoVal formula. Every lookup misses and the reference
  distribution is always regenerated. Do not describe it as a working
  cache; the matching machinery is kept only so the table can be
  repopulated cheaply.
- `pre_0.3.0` / `generator_version` fields keep backward compatibility
  with `.Rdata` files from older versions (index-counter starts at 0
  vs 1) — do not remove without understanding this.
- `R/zzz.R` declares
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html)
  for names that only exist at runtime after
  [`load()`](https://rdrr.io/r/base/load.html)-ing an `.Rdata` file
  (e.g. `base_faces`, `stimuli_params`, `p`, `seed`) or inside `foreach`
  loops (`obs`) — add new ones here to avoid `R CMD check` NOTEs.
- **[`load()`](https://rdrr.io/r/base/load.html) assigns into the
  calling function’s frame**, so a field in the `.Rdata` can silently
  replace a function argument of the same name. Every
  [`load()`](https://rdrr.io/r/base/load.html) site keeps copies of its
  arguments across the call; preserve that when adding either an
  argument or an `.Rdata` field.
- Parallelism is via base `parallel` + `doSNOW`/`foreach`, not newer
  alternatives (e.g. `future`) — match this pattern for new parallel
  code, and remember worker processes need `.packages='rcicr'` set on
  `foreach` calls. Tick progress bars from the parent with
  `.options.snow = progressOption(pb, cl)`; an in-body tick only updates
  a worker’s private copy.
