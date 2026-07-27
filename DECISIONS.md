# Decisions

Why `rcicr` is the way it is. Each entry records a decision, the reasoning behind it, and —
where it matters — what would have to be true to revisit it.

This exists because the other files cannot hold this. `NEWS.md` says what changed for users,
`BACKLOG.md` says what is left, `CLAUDE.md` says what the conventions are. None of them has
room for *why*, and the why is what gets lost: the measurement that ruled an option out, the
alternative that looked obvious and was wrong, the thing that looks like a bug and is not.

**Add to this file when a decision was not obvious** — when you rejected a plausible
alternative, measured something surprising, or deliberately did not fix something. Not every
change needs an entry. Entries are grouped by theme, not by date, and are meant to be edited
in place when they stop being true.

> This file was a chronological session log (`.session-log.md`) until 2026-07-27 and was
> reorganised by theme. The original narrative, with dates and intermediate states, is in
> git history up to commit `887aea4`.

---

## The constraint that shapes everything

**Researchers re-run old analysis scripts years later and publish what comes out.** Every
other decision here is downstream of that.

Concretely: do not change call syntax, argument meanings, or numeric output silently.
Deprecate rather than delete. Treat the `.Rdata` file as **append-only** — add fields, never
rename or repurpose. When a change does alter numeric output, that is not automatically
wrong, but it must be deliberate and written up in `NEWS.md` under "Reproducibility impact"
naming who is affected and what to do.

`tests/testthat/test-regression-baseline.R` is the enforcement: a golden master pinning the
default pipeline's output to its pre-fix values. **If it goes red, the change alters
researchers' results.** Document it; never edit the numbers to match.

---

## Numerics and the random number stream

### `purrr::rbernoulli()` was replaced with `runif()`, not `rbinom()`
`BACKLOG.md` originally recommended `stats::rbinom(n, 1, p)`. **That advice was wrong.**
`rbernoulli(n, p)` is internally `runif(n) > (1 - p)`; `rbinom` draws from the stream
differently — verified across 150 seed/probability combinations. Swapping it in would have
silently changed every reference distribution, and therefore every InfoVal, computed from a
given seed. The `runif` form is bit-identical to the old behaviour.

**The general rule this stands for: check the random *stream*, not just the distribution.**
Two functions with the same marginal distribution are not interchangeable inside a seeded
pipeline.

### `rowMeans(x, dims = 2)` was adopted despite not being bit-identical
The patch-averaging step in `generateNoiseImage()` moved from `apply(..., 1:2, mean)` to
`rowMeans(..., dims = 2)` — about 6x faster end-to-end. The two sum in a different order, so
they differ by roughly 1 ULP (~1e-19 absolute on pixel values of order 0.01).

Adopted because this was checked against an **independent oracle** — the average written as
an explicit triple loop, using neither function — across noise types, spatial scales and
seeds: both forms sit ~5.6e-17 from the oracle. Neither is "more correct". At the golden
master's own configuration the results came out bit-identical.

This is a different class of thing from the `rbinom` case above, and `NEWS.md` says so
explicitly: that one would have changed the random *stream*, a large systematic divergence;
this is floating-point summation order.

Note the bug that made this subtle: **`rowMeans()` on a 3-D array defaults to `dims = 1`**,
collapsing dims 2 *and* 3, and the short result is then silently recycled by `array()`. The
version originally submitted did exactly that — measured max deviation 0.21 against data
with SD ~0.01, not the ~1e-17 it claimed.

### `set.seed()` in `generateStimuli2IFC()` is load-bearing far beyond stimulus generation
`generateReferenceDistribution2IFC()` rebuilds the stimuli through `generateStimuli2IFC()`,
so that `set.seed(seed)` call lands *before* the simulation loop's `runif()` draws. That is
what makes InfoVal reproducible from a stimulus file alone, independent of ambient RNG state
and of `ncores`.

This was originally **emergent, not designed** — nobody chose it, and moving that one line
would have silently changed every InfoVal ever computed without touching
`computeInfoVal2IFC()`. It is now documented as a guarantee in
`?generateReferenceDistribution2IFC`, pinned by a test, and the `set.seed()` call carries a
comment saying what depends on it.

### `response_seed`, not `seed`
`generateReferenceDistribution2IFC()` gained a way to vary its null distribution. The
argument is called `response_seed` because `seed` is an *object in the `.Rdata` file*, and
the function re-saves its own frame — an argument of that name would overwrite the stimulus
seed and write it back, corrupting the record of how the stimuli were generated.

It seeds the **responses**, applied after the stimulus rebuild rather than forwarded into
it. Forwarding would rebuild a different stimulus set, so the null would describe stimuli
the participants never saw. `NULL` issues no `set.seed()` call at all, so the default path
is byte-identical to before rather than merely equivalent-looking.

`computeInfoVal2IFC(response_seed=)` additionally **forces regeneration and never caches**.
Without the forced regeneration it would be silently ignored on every file after the first,
since the generator writes `reference_norms` into the file. And caching a one-off Monte
Carlo check would silently redefine what every later InfoVal from that stimulus set means.

### `load()` assigns into the calling frame — check every new argument against saved names
An object in an `.Rdata` file silently overwrites a function argument of the same name.
`generateReferenceDistribution2IFC()` re-saved its whole frame, so the files it wrote
contained `rdata` and `ncores`; a second call then ignored the caller's `ncores` and wrote
back to the path recorded during the first call. Fixed at the source by excluding the
function's own arguments from the save, *and* defensively on read for files already written
by older versions.

**Any new argument added to a function that does `load(rdata)` must be checked against the
saved names.** This is also why `response_seed` is named as it is.

### The InfoVal formula is already correct — do not "fix" it
It matches the corrected version in the erratum to Schmitz et al. (2019): the Euclidean
norm, with the *k* constant supplied by R's `mad()` (`constant = 1.4826`). Fixed years ago,
now covered by a regression test.

---

## Things that look like bugs and are not

### `autoscale()` leaves `$combined` untouched — intentional
Flagged as a bug and changed; **Ron overruled, and the change was reverted.** `$combined`
must stay as the caller supplied it so a combination made before autoscaling survives the
call and existing scripts keep plotting the same image. `$scaled` carries the autoscaled
result.

The real user-facing problem was never the code — nothing said which field to look at. The
fix was documentation (`?autoscale`, plus the vignette) and two tests, one of which asserts
`$combined` comes back identical and warns against "fixing" it.

The trap worth knowing: after `batchGenerateCI()` (which scales `'none'` first) `$combined`
is an overlay of *unscaled* noise and looks almost blank. Build the overlay with
`(ci$scaled + ci$base) / 2`, which is what `save_as_pngs = TRUE` writes.

**How this stayed cheap to undo:** the change was flagged as arguable and filed under its own
"Behaviour change" heading rather than slipped in with unambiguous bug fixes. When a fix is
debatable, say so instead of relabelling it a bug.

### `computeInfoVal2IFC`'s test oracle mirrors the implementation — deliberate
The test recomputes `(norm - median)/mad`, the same expression as the implementation, so it
pins the *implementation* rather than the published definition. Left alone: the formula was
hand-checked against the Schmitz erratum and the golden master pins the resulting number, so
the risk is covered from two other directions. Recorded so nobody mistakes it for the
independent check it resembles. The genuinely independent oracle in the suite is in
`test-generateNoiseImage.R`.

### `_R_CHECK_LIMIT_CORES_` handling caps cores under check only
`default_ncores()` returns 2 when that variable is set, and `detectCores() - 1` otherwise.
CRAN policy caps packages at two cores during checks; the variable is only ever set by the
checker, so user behaviour at the console is unchanged. The known trade-off — `detectCores()`
over-subscribes under cgroup limits — is accepted deliberately rather than unnoticed.

### `generator_version` in old `.Rdata` files is not trustworthy
It was a hardcoded `'0.4.0'` string from 2016 until 1.1.0.9000, so every file written across
that entire range claims to be 0.4.0. `p$generator_version` held the real value all along.
Anything reading the field must treat `'0.4.0'` as "unknown, somewhere in that range", and
must accept both a character string and a `package_version` object. Compare with
`numeric_version()` semantics, never as text — `'0.10.0' < '0.4.0'` is `TRUE` as strings.

The `pre_0.3.0` compatibility path does **not** key off this field; it detects the old
`sinusoids`/`sinIdx` layout structurally. Just as well, given the field was lying.

---

## Testing

### Every new test must be shown to fail without its fix
Use **`git stash push -- R/`**, not a plain `git stash`. A plain stash reverts the new
*tests* as well, so the run passes vacuously and proves nothing.

For source changes, prefer `cp` backups in a scratchpad over git for mutation testing:
`git checkout <file>` discards unstaged work, which has destroyed an in-progress
implementation here. Guard each mutation with a `grep -q MUTANT` check confirming the patch
actually applied — a mutation that silently failed to apply looks exactly like a surviving
mutant.

### A test's title is a claim; check the assertions support it
Two `batchGenerateCI*` tests were titled "computes one CI per group" while asserting only
length, names and `dim`. **Grouping was never checked** — a bug feeding all trials to every
group passed green. When auditing, read the title and the `expect_*` calls as separate
things.

**Make the discriminating case explicit.** Every grouping and threshold test now also asserts
that the *wrong* answer differs (a positional split ≠ a by-column split), so a passing test
cannot be vacuous. Two vacuous assertions have been written here and caught by mutation
testing; both looked fine on the page. One compared `p`/`stimuli_params` in an `.Rdata` file
before and after a call that runs with `save_rdata = FALSE` and never writes them back — it
could not have failed.

### To test a rendered image, render onto a uniform background
Then "drew nothing" is `expect_length(unique(as.vector(png)), 1)` — a comparison rather than
an eyeball. `plotZmap`'s tests were three-quarters `file.exists()` before this; a function
writing a uniformly blank PNG passed them all.

### The recovery test uses a permutation null, not a parametric one
`test-recovery.R` simulates an observer with a known template and asserts `generateCI()`
recovers it. The statistic is `cor(vec(CI), vec(template))`; the null comes from permuting
response labels across trials, which preserves the noise images and the 1/-1 balance and
destroys only the response-to-stimulus pairing.

**Why not a t-test on Pearson r:** it would use df = n_pixels − 2 = 1022 and be wildly
anticonservative. The basis has 60 patches, so a 32×32 CI has effective dimensionality ~60,
not 1024, and neighbouring pixels are autocorrelated by construction — measured null
correlations reach |r| = 0.454, which a parametric test on 1022 df calls overwhelmingly
significant. Every null CI is built from the same basis, so the null inherits the identical
autocorrelation. That is what makes it correct, not merely convenient.

Sensitivity established by mutation: correct pairing 0.71, sign flip −0.74, responses
reversed −0.18, shifted one trial 0.24. **Documented limitation:** it cannot catch a
transformation applied consistently inside `generateNoiseImage()`, because the template is
built through the same function and the error cancels. The oracle test covers that.

### `NOT_CRAN` makes `skip_on_cran()` unverifiable the easy way
Both `devtools::test()` and `testthat::test_local()` set `NOT_CRAN=true` themselves, so
neither can be used to check that a skip fires.

### The package must be installed, not just `load_all()`-ed
`generateStimuli2IFC()` spawns parallel workers that each call `library(rcicr)`. Under
`load_all()` alone they fail with `there is no package called 'rcicr'`.

---

## Performance and parallelism

### `ncores == 1` runs in-process instead of building a one-worker cluster
`startBackend()` in `zzz.R` registers `doSEQ` when `ncores < 2`, so the same `%dopar%` loops
run in the current process and **no loop body changed**. The test suite went from 140s to 4s;
under `R CMD check`, tests went from `[8s/126s]` to `[8s/37s]` — eight seconds of CPU against
126 elapsed was worker startup, not computation.

**Safe because neither parallel loop draws random numbers** — parameters are drawn under
`set.seed()` before the loop — so there is no per-worker RNG stream to diverge. Verified, not
asserted: `ncores = 1` and `ncores = 2` give bit-identical stimulus parameters, noise basis
and CIs. Pinned by `test-parallel-equivalence.R`, which also documents what would make this
unsafe, so whoever adds a random draw inside those loops finds out.

### Parallelism stays on `parallel` + `doParallel`/`foreach`
Not `future` or `snowfall`. Issue #66 asked specifically for snowfall as a single-core
fallback; the outcome it wanted is what `registerDoSEQ()` provides, without a dependency.
Match the existing pattern for new parallel code, and remember workers need
`.packages = 'rcicr'` on `foreach` calls.

### Memory: issue #12 was solved by deletion, not by chunking
The issue proposed distributing stimulus matrices over multiple `.RData` files. The actual
cause was a preallocated `zeros(img_size, img_size, n_trials)` array — ~1.5 GB at the
defaults — that existed in the parent environment, so `foreach` exported a full copy to
*every* worker, each of which wrote one slice and discarded the rest. Each iteration now
allocates only its own trial's noise. Chunking was never needed.

---

## Packaging, CI and tooling

### `styler` and `lintr` are deliberately not in the pre-commit config
They would reformat nearly every file in one sweep and destroy `git blame`. Do it as one
clearly-labelled commit or not at all. The pre-commit hooks are minimal and
language-agnostic on purpose.

### `fail_ci_if_error: false` on the coverage workflow
Chosen over getting a Codecov token or deleting the workflow. **This does not make coverage
reporting work** — no token means no badge, no per-PR comments, and `codecov.yml`'s
thresholds stay inert. What it fixes is narrower and more useful: a red `main` now means the
package is broken. The `token:` input is left wired so adding the secret later restores
reporting.

### `.Rbuildignore` and `.gitignore` are not interchangeable
`R CMD build` works from the *working directory*, not from git, so a file that is git-ignored
but present on disk still ships in the tarball unless `.Rbuildignore` also excludes it.
Untracking a file is never a substitute for `.Rbuildignore`ing it.

Any new top-level file needs an `.Rbuildignore` entry or `R CMD check` raises a
"non-standard file/directory found at top level" NOTE. The file holds `^`-anchored
**regexes**, not globs.

`.Rbuildignore` itself was listed in `.gitignore` from 2016 until 2026 — an unintentional
RStudio-template leftover, so it never shipped to CI or other contributors. Do not re-add it.

### Check `git diff --stat main...HEAD` before opening a PR
One branch carried ~75,000 lines of committed `R CMD check` artifacts across three commits,
unnoticed because only the diffs of edited files were ever read. `pre-commit.ci` caught it.
Both patterns are now in `.gitignore` and `.Rbuildignore`.

### CI workflows only trigger on PRs targeting `main`
A stacked PR based on another branch gets pre-commit and nothing else — no `R CMD check`.
Retargeting an existing PR does **not** re-fire the workflows; close and reopen it.

---

## Releases and git

### Trunk-based with tags; no `develop` branch
The standard R-package layout. CRAN has no concept of a develop branch, this is a
single-maintainer package, and it would add a permanent second merge direction for no gain.
Feature branches → PR → squash onto `main`; releases are tags.

### `main` carries a `.9000` development version, and that is safe because of the tag
Right after a release `DESCRIPTION` goes to `<released>.9000`; the release commit drops it to
the clean number. **Build the CRAN tarball from the tag, never from `main` HEAD.** This is
what makes the suffix safe: `Version contains large components` is only a blocker if the
*submitted tarball* carries it, and a tarball built from a release tag never does. Do not
abandon the development-version convention to dodge that NOTE.

Tagging was missing entirely until 2026-07-27, and it cost real clarity: `main` had moved two
PRs past the 1.1.0 awaiting CRAN submission with nothing recording which tree that was.
`v1.1.0` is tagged retroactively at `a3904e8`.

### Squash merges to `main`
One commit per PR keeps history readable and makes reverting a whole change straightforward;
here a PR is typically one fix plus its test and its `NEWS.md` entry. The consequence is that
per-commit detail disappears, so **measurements, rejected alternatives and reproducibility
impact must live in the squash message or `NEWS.md`.**

This convention also absorbed the 75,000-line artifact mistake above: `main` only ever saw
the final tree, so no history rewrite was needed. It is repo convention, not enforced by
GitHub settings.

Deleting merged branches is part of it. Note `gh pr merge --delete-branch` already removes
the remote branch, and `git push origin --delete a b` is **atomic** — one stale name fails
the whole command and takes the other branch down with it.

### `NEWS.md` entries go largest-impact first
Changes to numeric output or return values, then behaviour changes, then bug fixes that only
ever produced errors, then message-only fixes. Someone who stops reading after three bullets
should have read the three that could change their results.

### `ChangeLog` gets a pointer entry, not a duplicate
`NEWS.md` supersedes it, but someone opening `ChangeLog` first would otherwise conclude the
last release is whatever it last recorded. Two changelogs both claiming authority is worse
than one that defers.

---

## Documentation

### The Medium walkthrough moved into a vignette, and the post stays up
The post is the tutorial users are sent to; it 403s to `R CMD check`, and outside the repo it
cannot execute at build time, so it goes stale silently. The premise proved itself within
minutes of checking: two of its lines no longer worked (`autoscale(cis, saveasjpegs = TRUE)`,
now `save_as_pngs`; and `install_github(..., ref = "development")`, a branch that does not
exist).

Keeping the post published is deliberate — nine years of inbound links. This is a move of the
canonical copy, not a takedown.

### Figures must be checked by looking at them
Three problems in the walkthrough vignette were caught only by viewing the output:

- **`image()` auto-stretches its input across the palette**, so every linear rescaling of the
  same image renders *identically*. This made a four-way scaling comparison meaningless until
  `zlim = c(0, 1)` was pinned.
- **`zmapmethod = "quick"` z-scores across the pixels of one image**, so its values are
  relative to that image's own spatial structure, not a null distribution. At small sizes they
  span only ±1.65 and the **default `threshold = 3` returns a blank map** — which reads as "no
  signal" but means "wrong ruler".
- A white-noise synthetic base drowned every figure; a crude Gaussian-blob face is legible and
  still obviously synthetic.

### Base images in tests and vignettes are always synthetic
`png::writePNG(matrix(runif(n*n), n, n), ...)` or a Gaussian blob — never a real photograph.
Avoids licensing and consent questions entirely.

### The `.Rdata` anatomy belongs in `README.md`
It is the only link between the two halves of the package, and nothing about a stimulus set
is recoverable without it. The field-by-field table there was written by inspecting a real
generated file rather than by reading the `save()` call — which is how `trial` was identified
as a leftover loop counter carrying no information.

---

## Corrections — claims that did not survive checking

Kept because each was about to be written somewhere a reader could act on it, and the pattern
is more useful than any individual case.

- **"29x faster."** The benchmark precomputed the weighted patch array *outside* the timed
  call, so it measured only the step that changed. End-to-end it is **6.1x** — exactly what
  the original contributor had reported. **Benchmark the whole call, not the step you
  changed.** Correcting it took a code comment, `NEWS.md`, a PR title and body, and two
  places on a public thread.
- **"Medium blocks non-browser user agents."** It blocks by *network origin* (datacenter
  IPs). `curl -A "Mozilla/5.0"` still 403s. This was one command from being disproved and was
  about to go into `cran-comments.md` for a CRAN reviewer.
- **"The manual ERROR/WARNING resolved itself."** That run passed `--no-manual`, which
  *skips* the check rather than passing it. **Never treat a `--no-manual` run as evidence the
  manual builds.**
- **"Dependencies cut to 14."** `DESCRIPTION` said 15. The number was inferred (27 − 13)
  rather than counted, in a commit whose purpose was fixing inaccurate docs. **Numbers copied
  from neighbouring prose are not verified numbers.**
- **"Every P1 item is fixed."** Two P1 items were open. Reworded to "every P1 *code* item".
- **`NEWS.md` saying an issue is "addressed" is not evidence it can be closed.** Two issues
  that looked closed by 1.1.0 were still reproducible when actually run: a progress bar
  emitting zero characters, and old `.Rdata` files still failing because an `exists()` guard
  was added for two sibling variables but not the third, three lines apart.
- **A PR's diff is not its state.** One contributor's diff showed broken code while the
  review thread already contained a correct version, posted years earlier and never pushed.
  Read the comments before summarising a PR.
- **Check CI before merging, not after.** One merge went ahead with a check still
  `in_progress`, on the strength of "it's all green" plus three of four passing.
- **`BACKLOG.md`'s own summary table drifted out of step with its item sections** — five
  items marked done in their sections were still listed as outstanding in the table a cold
  session reads first. Update both, or the table becomes a trap.

---

## Open decisions — these need Ron, and should not be settled unilaterally

- **`plotZmap(mask = ...)` is dead code** (`BACKLOG.md` item 23). The mask is validated,
  dimension-checked and booleanised, then never applied. Documented behaviour does nothing.
  Fixing it changes rendered output for anyone currently passing `mask`, so it is a behaviour
  change, not a plain bug fix — the `autoscale()`/`$combined` lesson applies. Fix under a
  "Behaviour change" heading, or deprecate the argument.
- **CRAN resubmission** (`BACKLOG.md` item 1). Archived 2021-06-08 for an undeliverable
  maintainer address; the code was never the problem and the address now works. Every
  mechanical blocker is cleared. win-builder and R-hub are deliberately unrun because they
  upload the tarball to third parties and email the maintainer — two `<!-- TODO -->` markers
  in `cran-comments.md` wait on them. **Ron submits personally**, since CRAN mails the
  maintainer address to confirm, which is the entire point of the resubmission.
- **The announcement post** (`BACKLOG.md` item 21). Drafts in `notes/`. Held until the CRAN
  outcome is known so it can say where the package lives. **Do not soften the reproducibility
  section** when editing: it states both that default results are unchanged *and* that two
  fixes genuinely change infoVal. "Nothing changed" would be untrue and less useful to the
  researcher it is written for.
