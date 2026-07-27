# Decisions

Why `rcicr` is the way it is: the measurement that ruled an option out, the alternative that
looked obvious and was wrong, the thing that looks like a bug and is not. `NEWS.md` holds what
changed for users, `BACKLOG.md` what is left, `CLAUDE.md` and `CONTRIBUTING.md` the
conventions — none of them has room for *why*.

**Add an entry when a decision was not obvious**: a plausible alternative rejected, something
surprising measured, something deliberately not fixed. Entries are grouped by theme, not by
date, and are edited in place when they stop being true.

> This was a chronological session log (`.session-log.md`) until 2026-07-27. The original
> narrative, with dates and intermediate states, is in git history up to `887aea4`.

---

## The constraint that shapes everything

**Researchers re-run old analysis scripts years later and publish what comes out.** Every
other decision here is downstream of that.

The rules it implies — no silent changes to call syntax, argument meanings or numeric output;
the `.Rdata` file is append-only; anything that does change numbers gets a `NEWS.md`
"Reproducibility impact" entry — are spelled out in `CONTRIBUTING.md` and enforced by the
golden master in `tests/testthat/test-regression-baseline.R`. **If it goes red, the change
alters researchers' results.** Document it; never edit the numbers to match.

---

## Numerics and the random number stream

### `purrr::rbernoulli()` was replaced with `runif()`, not `rbinom()`
`BACKLOG.md` recommended `stats::rbinom(n, 1, p)`. **That advice was wrong.** `rbernoulli(n,
p)` is internally `runif(n) > (1 - p)`, and `rbinom` draws from the stream differently —
verified across 150 seed/probability combinations. Swapping it in would have silently changed
every reference distribution, and therefore every InfoVal, computed from a given seed. The
`runif` form is bit-identical to the old behaviour.

**The rule this stands for: check the random *stream*, not just the distribution.** Two
functions with the same marginal distribution are not interchangeable in a seeded pipeline.

### `rowMeans(x, dims = 2)` was adopted despite not being bit-identical
The patch-averaging step in `generateNoiseImage()` moved from `apply(..., 1:2, mean)`, about
6x faster end-to-end. The two sum in a different order and so differ by roughly 1 ULP (~1e-19
on pixel values of order 0.01). Adopted because an **independent oracle** — the average
written as an explicit triple loop, using neither function — put both forms ~5.6e-17 away
across noise types, scales and seeds. Neither is "more correct", and at the golden master's
own configuration the results came out bit-identical. Unlike the `rbinom` case above this is
floating-point summation order, not a changed random stream; `NEWS.md` says so explicitly.

The subtlety that bit once: **`rowMeans()` on a 3-D array defaults to `dims = 1`**, collapsing
dims 2 *and* 3, and `array()` then silently recycles the short result. The version originally
submitted did exactly that — measured max deviation 0.21 against data with SD ~0.01, not the
~1e-17 it claimed.

### `set.seed()` in `generateStimuli2IFC()` is load-bearing far beyond stimulus generation
`generateReferenceDistribution2IFC()` rebuilds the stimuli through it, so that
`set.seed(seed)` lands *before* the simulation loop's `runif()` draws — which is what makes
InfoVal reproducible from a stimulus file alone, independent of ambient RNG state and of
`ncores`. This was **emergent, not designed**: moving one line would have silently changed
every InfoVal ever computed, without touching `computeInfoVal2IFC()`. It is now a documented
guarantee in `?generateReferenceDistribution2IFC`, pinned by a test, with a comment on the
`set.seed()` call saying what depends on it.

### `response_seed`, not `seed`
`seed` is an *object in the `.Rdata` file* and `generateReferenceDistribution2IFC()` re-saves
its own frame, so an argument of that name would overwrite the stimulus seed and write it
back, corrupting the record of how the stimuli were generated.

It seeds the **responses**, applied after the stimulus rebuild rather than forwarded into it:
forwarding would rebuild a different stimulus set, so the null would describe stimuli the
participants never saw. `NULL` issues no `set.seed()` call at all, so the default path is
byte-identical rather than merely equivalent-looking. `computeInfoVal2IFC(response_seed=)`
additionally **forces regeneration and never caches** — without the former it would be
silently ignored on every file after the first, since the generator writes `reference_norms`
into the file; and caching a one-off Monte Carlo check would redefine what every later InfoVal
from that stimulus set means.

### `load()` assigns into the calling frame — check every new argument against saved names
An object in an `.Rdata` file silently overwrites a function argument of the same name.
`generateReferenceDistribution2IFC()` re-saved its whole frame, so the files it wrote
contained `rdata` and `ncores`; a second call then ignored the caller's `ncores` and wrote
back to the path recorded during the first. Fixed at the source, by excluding the function's
own arguments from the save, *and* defensively on read for files older versions already wrote.
This is also why `response_seed` is named as it is.

### The InfoVal formula is already correct — do not "fix" it
It matches the erratum to Schmitz et al. (2019): the Euclidean norm, with *k* supplied by R's
`mad()` (`constant = 1.4826`). Fixed years ago, now pinned by a regression test.

---

## Things that look like bugs and are not

### `autoscale()` leaves `$combined` untouched — intentional
Flagged as a bug and changed; **Ron overruled, and the change was reverted.** `$combined` must
stay as the caller supplied it, so a combination made before autoscaling survives the call and
existing scripts keep plotting the same image; `$scaled` carries the autoscaled result. The
user-facing problem was never the code but that nothing said which field to look at — fixed in
`?autoscale` and the vignette, plus a test asserting `$combined` comes back identical and
warning against "fixing" it.

The trap worth knowing: after `batchGenerateCI()` (which scales `'none'` first) `$combined` is
an overlay of *unscaled* noise and looks almost blank. Build the overlay as `(ci$scaled +
ci$base) / 2`, which is what `save_as_pngs = TRUE` writes.

**What kept this cheap to undo:** the change was filed under its own "Behaviour change"
heading rather than slipped in with unambiguous bug fixes. When a fix is debatable, say so.

### `computeInfoVal2IFC`'s test oracle mirrors the implementation — deliberate
It recomputes `(norm - median)/mad`, the same expression as the implementation, so it pins the
*implementation* rather than the published definition. Left alone: the formula was
hand-checked against the Schmitz erratum and the golden master pins the resulting number, so
the risk is covered from two other directions. Recorded so nobody mistakes it for the
independent check it resembles — the genuinely independent oracle is in
`test-generateNoiseImage.R`.

### `_R_CHECK_LIMIT_CORES_` handling caps cores under check only
`default_ncores()` returns 2 when that variable is set and `detectCores() - 1` otherwise. Only
the checker ever sets it, so CRAN's two-core policy is met and behaviour at the console is
unchanged. The trade-off — `detectCores()` over-subscribes under cgroup limits — is accepted
deliberately, not unnoticed.

### `generator_version` in old `.Rdata` files is not trustworthy
It was a hardcoded `'0.4.0'` from 2016 until 1.1.0.9000, so every file written in that range
claims to be 0.4.0; `p$generator_version` held the real value all along. Anything reading the
field must treat `'0.4.0'` as "unknown, somewhere in that range", accept both a character
string and a `package_version`, and compare with `numeric_version()` semantics — as text,
`'0.10.0' < '0.4.0'` is `TRUE`. The `pre_0.3.0` compatibility path does **not** key off this
field; it detects the old `sinusoids`/`sinIdx` layout structurally. Just as well.

---

## Testing

### Mutation testing: `cp` backups, and check the mutant applied
`CONTRIBUTING.md` covers the basic check (`git stash push -- R/`, not a plain `git stash`).
Beyond it: prefer `cp` backups in a scratchpad over git, because `git checkout <file>`
discards unstaged work and has destroyed an in-progress implementation here. Guard each
mutation with a `grep -q MUTANT` check — one that silently failed to apply looks exactly like
a surviving mutant.

### Vacuous assertions have shipped here twice, and both looked fine on the page
Two `batchGenerateCI*` tests titled "computes one CI per group" asserted only length, names
and `dim`: **grouping was never checked**, and a bug feeding all trials to every group passed
green. The other compared `p`/`stimuli_params` in an `.Rdata` file across a call that runs
with `save_rdata = FALSE` and never writes them back, so it could not have failed. Mutation
testing caught both; reading did not. Hence every grouping and threshold test now also asserts
that the *wrong* answer differs (a positional split ≠ a by-column split).

### To test a rendered image, render onto a uniform background
Then "drew nothing" is `expect_length(unique(as.vector(png)), 1)` — a comparison rather than
an eyeball. `plotZmap`'s tests were three-quarters `file.exists()` before this; a function
writing a uniformly blank PNG passed them all.

### The recovery test uses a permutation null, not a parametric one
`test-recovery.R` gives a simulated observer a known template and asserts `generateCI()`
recovers it, scored by `cor(vec(CI), vec(template))`. The null permutes response labels across
trials, which preserves the noise images and the 1/-1 balance and destroys only the
response-to-stimulus pairing.

**Why not a t-test on Pearson r:** it would use df = n_pixels − 2 = 1022 and be wildly
anticonservative. The basis has 60 patches, so a 32×32 CI has effective dimensionality ~60,
not 1024, and neighbouring pixels are autocorrelated by construction — measured null
correlations reach |r| = 0.454. Every null CI is built from the same basis, so the null
inherits that exact autocorrelation; that is what makes it correct, not merely convenient.

Sensitivity established by mutation: correct pairing 0.71, sign flip −0.74, responses reversed
−0.18, shifted one trial 0.24. **Documented limitation:** it cannot catch a transformation
applied consistently inside `generateNoiseImage()`, because the template is built through the
same function and the error cancels. The oracle test covers that.

### `NOT_CRAN` makes `skip_on_cran()` unverifiable the easy way
Both `devtools::test()` and `testthat::test_local()` set `NOT_CRAN=true` themselves, so
neither can be used to check that a skip fires.

---

## Performance and parallelism

### `ncores == 1` runs in-process instead of building a one-worker cluster
`startBackend()` in `zzz.R` registers `doSEQ` when `ncores < 2`, so the same `%dopar%` loops
run in the current process and **no loop body changed**. The test suite went from 140s to 4s;
under `R CMD check`, from `[8s/126s]` to `[8s/37s]` — eight seconds of CPU against 126 elapsed
was worker startup, not computation.

**Safe because neither parallel loop draws random numbers** — parameters are drawn under
`set.seed()` before the loop — so there is no per-worker stream to diverge. Verified, not
asserted: `ncores = 1` and `ncores = 2` give bit-identical stimulus parameters, noise basis
and CIs, pinned by `test-parallel-equivalence.R`, which also documents what would make this
unsafe.

### Parallelism stays on `parallel` + `doParallel`/`foreach`
Not `future` or `snowfall`. Issue #66 asked specifically for snowfall as a single-core
fallback; the outcome it wanted is what `registerDoSEQ()` provides, without a dependency.

### Memory: issue #12 was solved by deletion, not by chunking
The issue proposed spreading stimulus matrices over multiple `.RData` files. The actual cause
was a preallocated `zeros(img_size, img_size, n_trials)` array — ~1.5 GB at the defaults —
living in the parent environment, so `foreach` exported a full copy to *every* worker, each of
which wrote one slice and discarded the rest. Each iteration now allocates only its own
trial's noise.

---

## Packaging, CI and tooling

### `styler` and `lintr` are deliberately not in the pre-commit config
A one-sweep reformat would destroy `git blame`. If the package is ever styled it goes in as
one clearly-labelled commit of its own, never as a side effect of other work — which is why
the hooks that do run are minimal and language-agnostic.

### `fail_ci_if_error: false` on the coverage workflow
Chosen over getting a Codecov token or deleting the workflow. **This does not make coverage
reporting work** — no token means no badge, no per-PR comments, and `codecov.yml`'s thresholds
stay inert. What it buys is narrower and more useful: a red `main` now means the package is
broken. The `token:` input is left wired so adding the secret later restores reporting.

### `.Rbuildignore` was itself git-ignored for a decade
The mechanics are in `CLAUDE.md` and `CONTRIBUTING.md`. What is easy to reintroduce:
`.Rbuildignore` was listed in `.gitignore` from 2016 until 2026, an unintentional
RStudio-template leftover, so it never shipped to CI or to other contributors at all. Do not
re-add it.

### The 75,000-line artifact commit
One branch carried ~75,000 lines of committed `R CMD check` artifacts across three commits,
unnoticed because only the diffs of edited files were ever read — hence the `git diff --stat
main...HEAD` habit in `CONTRIBUTING.md`. `pre-commit.ci` caught it, and both patterns are now
in `.gitignore` and `.Rbuildignore`.

### CI workflows only trigger on PRs targeting `main`
A stacked PR based on another branch gets pre-commit and nothing else — no `R CMD check`.
Retargeting an existing PR does **not** re-fire the workflows; close and reopen it.

---

## Releases and git

`CLAUDE.md` states the conventions and their rationale — trunk-based with tags and no
`develop` branch, `.9000` on `main` between releases, squash merges, delete merged branches,
`NEWS.md` ordered largest-impact first. What those entries do not say:

- **The `.9000` suffix is safe only because the tarball is built from the tag.** `Version
  contains large components` is a blocker if and only if the *submitted* tarball carries it.
  Do not abandon the development-version convention to dodge that NOTE.
- **Untagged releases cost real clarity.** `main` had moved two PRs past the 1.1.0 awaiting
  CRAN submission with nothing recording which tree that was; `v1.1.0` is tagged retroactively
  at `a3904e8`. That is what the tagging convention buys.
- **Squashing is what made the 75,000-line artifact mistake above cheap** — `main` only ever
  saw the final tree, so no history rewrite was needed.
- `git push origin --delete a b` is **atomic**: one stale name fails the whole command and
  takes the other branch down with it. (`gh pr merge --delete-branch` has already removed the
  remote branch anyway.)

### `ChangeLog` gets a pointer entry, not a duplicate
`NEWS.md` supersedes it, but someone opening `ChangeLog` first would otherwise conclude the
last release is whatever it last recorded. Two changelogs both claiming authority is worse
than one that defers.

---

## Documentation

### The Medium walkthrough moved into a vignette, and the post stays up
The post is the tutorial users are sent to; it 403s to `R CMD check`, and outside the repo it
cannot execute at build time, so it goes stale silently — a premise that proved itself within
minutes of checking, since two of its lines no longer worked: the `saveasjpegs` argument to
`autoscale()`, now `save_as_pngs`, and `install_github(..., ref = "development")`, a branch
that does not exist. The post stays published for nine years of inbound links — this moves the
canonical copy, it is not a takedown.

### Figures must be checked by looking at them
Three problems in the walkthrough vignette were caught only by viewing the output:

- **`image()` auto-stretches its input across the palette**, so every linear rescaling of the
  same image renders *identically* — which made a four-way scaling comparison meaningless
  until `zlim = c(0, 1)` was pinned.
- **`zmapmethod = "quick"` z-scores across the pixels of one image**, so its values are
  relative to that image's own spatial structure, not a null distribution. At small sizes they
  span only ±1.65 and the **default `threshold = 3` returns a blank map** — which reads as "no
  signal" but means "wrong ruler".
- A white-noise synthetic base drowned every figure; a crude Gaussian blob is legible and
  still obviously synthetic.

### Base images in tests and vignettes are always synthetic
A `runif()` matrix or a Gaussian blob, never a real photograph — it avoids licensing and
consent questions entirely, which is worth more than a realistic-looking figure.

### The `.Rdata` anatomy belongs in `README.md`
It is the only link between the two halves of the package, and nothing about a stimulus set is
recoverable without it. The field-by-field table there was written by inspecting a real
generated file rather than by reading the `save()` call — which is how `trial` was identified
as a leftover loop counter carrying no information.

---

## Corrections — claims that did not survive checking

Kept because each was about to be written somewhere a reader could act on it, and the pattern
is more useful than any individual case. **Drop an entry when its pattern is already carried
by another one, or when the mistake has become impossible to repeat — not merely when the
particular fact has been corrected.**

- **"29x faster."** The benchmark precomputed the weighted patch array *outside* the timed
  call, so it measured only the step that changed; end-to-end it is **6.1x**, exactly what the
  original contributor had reported. **Benchmark the whole call, not the step you changed.**
  Correcting it took a code comment, `NEWS.md`, a PR title and body, and two places on a
  public thread.
- **"The manual ERROR/WARNING resolved itself."** That run passed `--no-manual`, which *skips*
  the check rather than passing it. **Never treat a `--no-manual` run as evidence the manual
  builds.**
- **"Dependencies cut to 14."** `DESCRIPTION` said 15. The number was inferred (27 − 13)
  rather than counted, in a commit whose purpose was fixing inaccurate docs. **Numbers copied
  from neighbouring prose are not verified numbers.**
- **`NEWS.md` saying an issue is "addressed" is not evidence it can be closed.** Two issues
  that looked closed by 1.1.0 were still reproducible when actually run: a progress bar
  emitting zero characters, and old `.Rdata` files still failing because an `exists()` guard
  was added for two sibling variables but not the third, three lines apart.
- **A PR's diff is not its state.** One contributor's diff showed broken code while the review
  thread already held a correct version, posted years earlier and never pushed. Read the
  comments before summarising a PR.
- **Check CI before merging, not after.** One merge went ahead with a check still
  `in_progress`, on the strength of "it's all green" plus three of four passing.
- **`BACKLOG.md`'s own summary table drifted out of step with its item sections** — five items
  marked done in their sections were still listed as outstanding in the table a cold session
  reads first. Update both, or the table becomes a trap.

---

## Open decisions — these need Ron, and should not be settled unilaterally

- **`plotZmap(mask = ...)` is dead code** (`BACKLOG.md` item 23). The mask is validated,
  dimension-checked and booleanised, then never applied. Fixing it changes rendered output for
  anyone currently passing `mask`, so it is a behaviour change rather than a plain bug fix —
  the `autoscale()`/`$combined` lesson applies. Fix under a "Behaviour change" heading, or
  deprecate the argument.
- **CRAN resubmission** (`BACKLOG.md` item 1). Archived 2021-06-08 for an undeliverable
  maintainer address; the code was never the problem and the address now works. Every
  mechanical blocker is cleared. win-builder and R-hub are deliberately unrun because they
  upload the tarball to third parties and email the maintainer — two `<!-- TODO -->` markers
  in `cran-comments.md` wait on them. **Ron submits personally**, since CRAN mails the
  maintainer address to confirm, which is the entire point of the resubmission.
- **The announcement post** (`BACKLOG.md` item 21). Drafts in `notes/`, held until the CRAN
  outcome is known so it can say where the package lives. **Do not soften the reproducibility
  section** when editing: it states both that default results are unchanged *and* that two
  fixes genuinely change infoVal. "Nothing changed" would be untrue and less useful to the
  researcher it is written for.
