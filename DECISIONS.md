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
other decision here is downstream of that. The rules it implies, and the golden master that
enforces them, are in `CONTRIBUTING.md`.

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
Four choices, each forced by something specific:

- **The name.** `seed` is an *object in the `.Rdata` file*, and
  `generateReferenceDistribution2IFC()` re-saves its own frame, so an argument of that name
  would overwrite the stimulus seed and write it back — corrupting the record of how the
  stimuli were generated.
- **It seeds the responses, applied after the stimulus rebuild**, not forwarded into it.
  Forwarding would rebuild a different stimulus set, so the null would describe stimuli the
  participants never saw.
- **`NULL` issues no `set.seed()` call at all**, so the default path is byte-identical rather
  than merely equivalent-looking.
- **`computeInfoVal2IFC(response_seed=)` forces regeneration and never caches.** Without the
  first it would be silently ignored on every file after the first, since the generator writes
  `reference_norms` into the file; caching a one-off Monte Carlo check would redefine what
  every later InfoVal from that stimulus set means.

### 4096 → 4092 parameters: why old stimulus files cannot be regenerated from their seed
rcicr 0.3.0 (2015-01-23) cut the per-trial random draws from 4096 to 4092. 4092 is the real
patch count — 6 orientations × 2 phases × `sum(4^0..4^4)` — while 4096 was a round `2^12`
over-allocation, so four contrasts were drawn per trial that no patch index ever referred to.
`ChangeLog` calls them "redundant" and says the change "does not affect anything else".

**That last claim is true for analysis and false for regeneration**, which is not obvious and
matters. Analysing a pre-0.3.0 file reads its *stored* parameters, and dropping four unused
columns changes nothing. But `runif()` is a stream: at the same seed, trial 1 gets identical
values either way and **every later trial is shifted by four draws** — verified, trial 1
identical, trials 2 and 3 not. So a pre-0.3.0 stimulus set can be re-analysed exactly and
cannot be re-created from its seed.

The same release also fixed `sinIdx` counting from 0 rather than 1, which is what the
`pre_0.3.0` flag exists for. The two are independent and both live in that boundary.

The truncation this left behind in `generateCI()` had a dead branch for eleven years: the
single-trial path tested for a length of 4092 and truncated to 4092, so it could never fire
on the 4096-length input it was written for (`BACKLOG.md` item 28). **A backward-compatibility
path that nothing exercises is indistinguishable from one that works** — this one had no test
until 2026-07-28, and neither did the `sinusoids`/`sinIdx` path, which was also broken.

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

A third instance was caught *before* shipping, when the z-map golden master was written on
2026-07-28, and it shows the shape to watch for. `zmapmethod = "quick"` ends in `scale()`, so
its z-map has sum 0 and sd 1 **by construction** — the obvious summary statistics to pin, and
both worthless as value checks. Mutating `sigma` from 3 to 4, a real change to the output, left
sum and sd bit-identical while `sum(abs())`, `min`, `max` and every individual cell moved. The
rule: **before pinning a summary statistic, ask what the transformation guarantees about it.**
Anything a normalisation fixes carries no information about the values that went in. They are
still asserted, labelled as a check that the standardisation happened, alongside statistics
that actually vary.

### To test a rendered image, render onto a uniform background
Then "drew nothing" is "the image is one flat value" — a comparison rather than an eyeball.
`plotZmap`'s tests were three-quarters `file.exists()` before this; a function writing a
uniformly blank PNG passed them all.

Count distinct values over the **colour channels only**, never the raw array. Whether a PNG
device writes an alpha channel is a property of the graphics backend, not of what was drawn:
cairo (Linux, Windows) writes RGB, macOS quartz writes RGBA. An opaque alpha plane is a solid
block of 1s, so it adds a second distinct value to the array while nothing has been painted.
The original `expect_length(unique(as.vector(png)), 1)` therefore measured the backend as much
as the drawing, and `test-plotZmap.R` now drops any alpha plane before counting.

This is not theoretical: it is the only assertion that failed when the suite first ran on
macOS (R-hub R-devel, 2026-07-28) — 220 passed, that one reported 2 distinct values where it
wanted 1. It was also the only assertion in the suite counting distinct values rather than
comparing two renders, which is why nothing else moved; every image-to-image comparison stays
green because both sides gain the same plane. Prefer the comparison form for new tests.

The replacement keeps the full detection power — verified against synthetic 1-, 2-, 3- and
4-channel images, painted and unpainted: every unpainted layout collapses to one value, every
painted one still yields two.

**The absolute pixel value is device-dependent too, and must not be asserted.** The first
attempt at this fix also pinned the flat value to the background grey, on the reasoning that
"uniform" alone would admit a `plotZmap()` flooding the image with a colour of its own. That
assertion passed on cairo and failed on macOS, which renders `bgimage` 0.5 to roughly 0.573
where cairo gives 0.502 — colour management, and precisely the same mistake the alpha fix was
correcting, made one line below it. The check is now an *ordering*: the same render over a
darker background must come out darker. Ordering survives any monotone transfer function,
which is all a colour pipeline can sanely be.

The general form, and the reason this cost two CI rounds rather than one: **when a test reads
pixels back from a graphics device, every absolute property of those pixels belongs to the
device.** Channel count and value are both device-dependent; only relationships between
renders — this differs from that, this is darker than that — are portable.

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

### The release gate runs the old code; the golden master only re-runs ours
`test-regression-baseline.R` pins values *this repository computed for itself*. That makes it
self-referential in one specific way: it can only catch drift away from the moment the numbers
were written down. Had a P0 fix already changed results before the baseline was recorded, the
baseline would have pinned the changed values and passed green forever.
`tools/compare-release-output.R` closes that gap by installing the reference commit into a
temporary library and running both versions over the same battery — it is the only thing here
that executes the old code. The two are complements, not substitutes: the golden master is
cheap enough to run on every commit, the gate costs two package installs and minutes of
compute, so it runs `--quick` on PRs and in full at release.

The battery is snapshotted into the temp directory before either side runs, so editing the
working copy mid-run cannot leave the two sides comparing different things — which happened
here once, and produces a "difference" that is purely an artefact.

### The v1.0.1 reference is pinned; the previous release is a *second* run, not a replacement
The obvious move once a release is green is to make it the new reference. It is wrong. Each
release would then be compared only against its predecessor, and a tree could walk away from
the published numbers one tolerated epsilon at a time, every step of the walk "identical to the
last release". The literature was produced with v1.0.1, so that comparison is the one that
protects it, and it stays pinned at `v1.0.1` (tagged retroactively at `b6ab269`, so the
default reads as a version rather than a bare SHA).

The second run (`--ref=v1.1.0`) answers a different and also useful question — did anything
break since the last release — and reaches further, because v1.0.1 *crashes* on calls that
v1.1.0 returns numbers for. Hence `EXPECTED` entries name the reference they apply to: a
deviation from v1.0.1 is not a deviation from v1.1.0, and an entry that fired for one would be
reported as stale by the other.

### The battery stops where the reference version crashes
Measured 2026-07-28 against v1.0.1 on R 4.3.3: it can produce a z-map **only** at 512px with
`zmapdecoration = TRUE`. Undecorated it dies in `if (bgimage != '')` ("the condition has length
> 1"), and at 64 and 128px in `plot.new()` ("figure margins too large" — 64px decorated fails
earlier still, on the `plt` graphical parameter). `mask` is unusable for the same class of
reason. All are fixed here, and none of them can be gated against v1.0.1: **a fix that turns a
crash into a number has no old number to compare against.** Those paths are covered by the test
suite, and by the gate only from v1.1.0 onward — the `SINCE` table in `tools/compare-harness.R`
records which extras need which reference.

### Tolerances: 8 ULP scaled to the values, plus an 8-bit pixel count
A flat `.Machine$double.eps` is the right bar for a classification image (values around 0.01)
and far too tight for a z-map, whose values run to several units and whose one-ULP steps are
therefore ~4× larger. So the tolerance is `8 * eps * max(1, max(abs(reference)))`: a few ULP
*of the values involved*. Anything that could only come from a changed random stream, algorithm
or file format — patch indices, drawn stimulus parameters, the base image, the stimulus PNGs —
is required to be bit-identical instead, with no tolerance at all.

**Why 8 and not 1.** These are not single passes over the data. A z-map is a convolution over
262,144 pixels followed by a standardisation over the same 262,144 values, so a 3.5e-18
difference in the classification image arrives as **1.33e-15** — measured at 512px, comparing
the current tree against v1.0.1. One ULP of the largest value (9.0e-16) rejects that; eight
accepts it and still sits three orders of magnitude below anything observable.

Widening a tolerance to make a run pass is exactly the move `CONTRIBUTING.md` warns against, so
the distinction matters: what protects this comparison is not the numeric tolerance but the two
**exact** checks beside it — 0 of N pixels may differ once quantised to 8 bits, and the NA
pattern must match cell for cell. The z-map sigma bug moved 1,282 cells across the threshold
and was caught by the NA check; every numeric tolerance considered here would have passed it.

Image-like outputs additionally have to survive quantisation: 0 of N pixels may differ once
rounded to 8 bits. That is the check that answers the question a researcher actually has, which
is whether the PNG they publish changes, and it is why a 1.11e-16 difference in the CI is
reported as a pass rather than argued about.

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

### If the package is ever run through `styler`, it goes in as a commit of its own
Never as a side effect of other work. (Why the hooks stay minimal and language-agnostic is in
`CONTRIBUTING.md`.)

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

### R-hub runs on `workflow_dispatch` only, never on push
The R-hub v2 workflow is the stock file `rhub::rhub_setup()` writes, kept unmodified so it can
be refreshed from upstream. It is left trigger-on-demand because R-hub answers a question that
only arises at release time — does this build under CRAN's own compiler flags and on the
odd platforms nobody develops on — and running it per-PR would spend a matrix of platforms on
a question nobody asked yet, competing with the ~11-minute reproducibility gate that *is*
required on every code PR.

This entry used to justify that by saying the everyday answer was "already covered by
`R-CMD-check.yaml` on release and devel". That was wrong, and expensively so: the matrix
varied only the *R version* and pinned `runs-on: ubuntu-latest`, so what it actually covered
was one platform twice. The hole stayed invisible until the first R-hub dispatch, on
2026-07-28 while preparing the CRAN resubmission, failed a `plotZmap()` test on macOS that
had been green on Linux since `b7cb6d9` — see the alpha-channel entry under Testing. CRAN
builds on macOS, so a platform-specific failure would have landed as a submission ERROR.
Both platforms CRAN gates on are now in the everyday matrix on release — `macos-latest`, then
`windows-latest` — and R-hub no longer stands in for per-platform coverage at all.

Windows was added on different grounds from macOS, and the distinction is the useful part.
macOS had never been run; Windows had, and passed, on that same R-hub dispatch and on
win-builder twice. So Windows was never an untested hole — it was a *late* one, checked only
at release, which is the same lateness the macOS row exists to fix. Adding it also puts the
pinned z-map values in `test-regression-baseline.R` in front of a third platform and a
different BLAS.

The lesson generalises past this workflow: "already covered by X" is a claim about what X
runs, and it is worth reading X to check rather than repeating.

Adding to that matrix is safe; renaming in it is not. `ubuntu-latest (release)` and
`ubuntu-latest (devel)` are required status checks matched **by name**, so a rename makes the
required check never report — which GitHub reads as pending forever, blocking every PR. Same
trap as the `paths:` filter on the reproducibility gate.

The cost of `workflow_dispatch`-only is that the workflow is invisible until it reaches the
**default branch** — GitHub offers no "Run workflow" button for a file that exists only on a
feature branch. So it merges to `main` ahead of the release it is meant to check, not
alongside it.

### The 75,000-line artifact commit
One branch carried ~75,000 lines of committed `R CMD check` artifacts across three commits,
unnoticed because only the diffs of edited files were ever read — hence the `git diff --stat
main...HEAD` habit in `CONTRIBUTING.md`. `pre-commit.ci` caught it, and both patterns are now
in `.gitignore` and `.Rbuildignore`.

---

## Releases and git

The conventions and their rationale are in `CLAUDE.md` — trunk-based with tags and no
`develop` branch, `.9000` on `main` between releases, tarballs built from the tag, squash
merges, delete merged branches, `NEWS.md` ordered largest-impact first. Two things learned
the hard way that those entries do not cover:

- **Squashing is what made the 75,000-line artifact mistake above cheap** — `main` only ever
  saw the final tree, so no history rewrite was needed.
- `git push origin --delete a b` is **atomic**: one stale name fails the whole command and
  takes the other branch down with it. (`gh pr merge --delete-branch` has already removed the
  remote branch anyway.)

### External checks run on the release branch; the tag and `.9000` still precede acceptance
At 1.2.1 the order was: merge the release PR, tag, publish the GitHub release, bump to
`.9000` — and only *then* run win-builder and R-hub. So the check results landed on a
development tree, and the release commit and its tag recorded no evidence that the tree they
name had passed anything. `usethis:::release_checklist()` puts `check_win_devel()` and
"Update `cran-comments.md`" under **Prepare for release**, before the version is finalised,
and leaves the tag (`use_github_release()`) and the `.9000` reopen (`use_dev_version()`)
until after "Wait for CRAN… Accepted". The checks now run in step 2 of
`CONTRIBUTING.md` → "Releasing", against the release branch.

**The circularity that appears to force the old order is not real.** It looks as though
results cannot be recorded in the commit they describe, because writing them changes the
tree and so invalidates the tarball that was checked. But `cran-comments.md` is
`.Rbuildignore`d and verifiably absent from the built tarball, so editing it produces the
same tarball byte for byte. It is not a package artifact at all — it is the text pasted into
the "Optional comment" field of the CRAN submission form.

**What is deliberately *not* copied from usethis: tagging and reopening `.9000` still happen
before CRAN accepts.** Two reasons. The `.9000` suffix is what selects `--quick` over the
full ~20-minute battery, so holding `main` at a clean version across the weeks CRAN can take
would run the full gate on every unrelated PR in that window. And a rejection means shipping
X.Y.Z+1 anyway, leaving the tag naming a tree that was never on CRAN — already true of
`v1.2.0`, so the repository tolerates that outcome cheaply.

This variant is in one respect tighter than usethis's: that checklist runs win-builder while
`DESCRIPTION` still carries `.9000` and never re-checks after `use_version()`, so the
recorded results are for a different version string than the one submitted. Checking the
release branch, where the version is already clean, avoids that.

Not applied retroactively to 1.2.1. `cran-comments.md` never being in the tarball is exactly
what makes a 1.2.2 unnecessary: the already-checked `rcicr_1.2.1.tar.gz` is what should be
submitted, and re-cutting a release to carry updated comments would produce an identical
tarball at the cost of a fresh round of external checks.

### GitHub releases carry notes only — the built tarball is not attached
Asked and settled at the 1.2.1 release. A release page already offers "Source code (tar.gz)",
which GitHub generates from the tag, and that is **not** the same artifact `R CMD build`
produces: the git archive has no built vignettes in `inst/doc/` and does carry every
`.Rbuildignore`d development file. Attaching the real tarball would put two different
"source" downloads on one page, which is a support question waiting to happen, and a second
artifact that can drift from the tag it claims to be.

The convention this follows is r-lib's and the tidyverse's: CRAN hosts the tarball, and
`remotes::install_github('rdotsch/rcicr@vX.Y.Z')` covers everyone else. That rcicr is
currently *off* CRAN is the strongest argument the other way, and it was considered — it is
temporary, and it does not outweigh publishing an artifact the project would then have to
keep consistent at every release.

What this deliberately gives up: the tarball is **not byte-reproducible from the tag**,
because `R CMD build` stamps `Packaged: <timestamp>; <user>` into `DESCRIPTION`. So the exact
bytes sent to CRAN are not archived anywhere and can only be rebuilt, into something
marginally different. That is accepted: what has to be reproducible is the *tree*, which the
tag pins exactly, and which is what both the release gate and `R CMD check` actually operate
on. Nothing about a result depends on the timestamp in a tarball header.

### `ChangeLog` gets a pointer entry, not a duplicate
`NEWS.md` supersedes it, but someone opening `ChangeLog` first would otherwise conclude the
last release is whatever it last recorded. Two changelogs both claiming authority is worse
than one that defers.

---

## Documentation

### The Medium walkthrough moved into a vignette, and the post stays up
A tutorial outside the repo cannot execute at build time, so it goes stale silently. The
premise proved itself during the port: two of the post's lines no longer worked —
`autoscale()`'s `saveasjpegs` argument, now `save_as_pngs`, and `install_github(..., ref =
"development")`, a branch that does not exist. The post stays published for nine years of
inbound links; this moves the canonical copy, it is not a takedown.

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

### When the docs and the code disagreed about `mask`, the code was the contract
Both `?plotZmap` and `?generateCI` said a *matrix* masks where the cell is `1`/`TRUE` while a
*PNG* masks where it is black (`0`) — two opposite conventions in one sentence. `applyMask()`
has always done `mask_matrix == 0` unconditionally, so the matrix half was simply wrong.

**The documentation was corrected to the code, not the reverse**, because users build masks
against observed behaviour: someone whose masked CI came out right has a mask encoding `0` =
masked, whatever the page said. Changing the code to match the prose would have silently
inverted every existing mask — the one outcome with real cost, since a mask covers a face and
an inverted one looks plausible.

**This nearly shipped backwards.** The first version of the `plotZmap` fix believed the docs
and branched on provenance, masking where `0` for a PNG and where `TRUE` for a matrix. That
would have made the same mask remove complementary halves in the two functions. It was caught
by Ron asking whether the two forms should really be opposite, and the answer was in the
sibling implementation rather than in any documentation. **When two documented conventions
conflict, run the one that has been executing for a decade.**

### `plotZmap(mask = ...)` was applied rather than deprecated
It had been documented since 2016 and never worked: commit `18e07cb` landed the import half as
*"add mask import ... (todo: applying the mask)"*, the todo was never picked up, and the
feature was implemented for `generateCI()` two weeks later instead. The pre-refactor monolith
at `d93cb2b^` confirms nothing was lost in the 2017 file split.

Filed as a behaviour change, but the risk that framing guards against could not materialise:
nothing has ever been masked, so no published z-map depends on the current output. GitHub code
search finds no caller of `plotZmap()` outside this repository — weak evidence on its own,
since analysis scripts usually live in OSF supplements, but enough alongside the rest.
Deprecating would have removed a documented feature on the grounds that it was never built.

**The generalisable part:** "this changes rendered output" is not by itself an argument for
leaving something broken. Ask who is currently *relying* on the broken output. Here the answer
was nobody, and the entry sat open a day longer than it needed to.

### Base images in tests and vignettes are always synthetic
Avoiding licensing and consent questions entirely is worth more than a realistic-looking figure.

### The `.Rdata` anatomy belongs in `README.md`
It is the only link between the two halves of the package, and nothing about a stimulus set is
recoverable without it. The field-by-field table there was written by inspecting a real
generated file rather than by reading the `save()` call — which is how `trial` was identified
as a leftover loop counter carrying no information.
