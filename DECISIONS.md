# Decisions

Why `rcicr` is the way it is: the measurement that ruled an option out, the alternative that
looked obvious and was wrong, the thing that looks like a bug and is not. `NEWS.md` holds what
changed for users, the issue tracker what is left, `AGENTS.md` and `CONTRIBUTING.md` the
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
An earlier note recommended `stats::rbinom(n, 1, p)`. **That advice was wrong.** `rbernoulli(n,
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

**The 0-based path in `generateNoiseImage()` looks like it corrupts the whole image and does
not.** When `min(patchIdx) == 0`, `params[p$patchIdx]` drops every 0-indexed cell (R drops
rather than `NA`s a 0 subscript), returning a vector too short for the patch array, which
`array()` then recycles — the setup for a general misalignment. It is harmless because the
same counter-from-0 leaves the *last* patch layer never written, so `patchIdx == 0` iff
`patches == 0` (verified at every `nscales`), and those zero-index cells sit contiguously at
the end of the column-major order. The recycled values land only where the patch is
identically zero and are multiplied away. Measured: the warning-branch output equals the
honest "one patch not shown" result exactly (max abs diff 0) across 36 size/nscales/seed
combinations, so the warning is accurate, not an understatement. Do **not** "fix" the recycle
by 1-offsetting the index — that would change which sinusoid is dropped, altering the CI of
every genuinely pre-0.3.0 file. `test-generateNoiseImage.R` pins the equality so the masking
cannot be broken silently.

The truncation this left behind in `generateCI()` had a dead branch for eleven years: the
single-trial path tested for a length of 4092 and truncated to 4092, so it could never fire
on the 4096-length input it was written for. **A backward-compatibility
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

### `base_face_files` validation rejects two inputs that used to run
Duplicate names, and a list element that is not a single file name, now stop the call. Both
ran before, so this is a considered exception to the constraint above rather than a plain bug
fix — and it is allowed because neither ever produced what the call asked for. `list(face =
'a.png', face = 'b.png')` looks each name up by string, so it wrote **one** set of stimuli,
from `a.png`, and nothing from `b.png`; verified against the pre-fix code, four files for two
trials. A script that "worked" this way was generating stimuli for a base image it never
named. Every other new check rejects input that already failed, just further downstream —
inside a parallel worker, as `attempt to select less than one element in get1index`.

**Readability is checked by attempting the read, not by `file.access()`.** The obvious version
— `file.access(filename, 4)` — is documented as unreliable on Windows and on networked
filesystems, so it can reject a file that reads fine. `png::readPNG()` / `jpeg::readJPEG()`
inside `tryCatch()` cannot be wrong about it, and keeping the reader's own message means the
distinction between "not a PNG" and "permission denied" survives instead of being flattened
into one guess.

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

### `return_as_dataframe = TRUE` returns one noise image per trial, not per trial × base image
`generateStimuli2IFC()`'s early `return()` (`R/generateStimuli2IFC.R:231`) sits *inside* the
per-base-image loop at `:193`, so with several base images it fires on the first and the rest
never run. Under the default
`use_same_parameters = TRUE` that is correct — every base image shares one parameter set, so one
noise image per trial is all there is — and the returned frame has one column per trial, so it
could not represent the alternative anyway. Documented on `@param return_as_dataframe` rather
than changed.

InfoVal is unaffected, checked rather than assumed:
`generateReferenceDistribution2IFC()` is the only in-package caller, it never passes
`use_same_parameters`, and the first base image's parameters come from the same leading block of
the RNG stream either way — measured identical, max absolute difference 0. Widening the frame to
trial × base image would change the return shape, so it needs a **new argument**, never a
redefinition.

### Repopulating `ref_lookup` costs four measurements — and the two halves stand or fall together
`AGENTS.md` covers what the table is (not a cache, empty since 2018, every lookup misses). What
belongs here is the way out, because either half done alone is worse than the status quo.

Repopulating means measuring four numbers: `median(reference_norms)` and `mad(reference_norms)`
under the current formula for seed 1, 512px, 10000 iterations at 100/300/500/1000 trials, one
`generateReferenceDistribution2IFC()` run each. The alternative is deleting the ~55 lines of
matching and prompt machinery. **Do not do half of either** — delete the machinery while
intending to re-measure and the feature becomes unrecoverable rather than merely dormant.

---

## Testing

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

### `computeInfoVal2IFC`'s test oracle mirrors the implementation, and is kept anyway
`test-computeInfoVal2IFC.R` recomputes `(norm(ci, "f") - median(reference_norms)) /
mad(reference_norms)` — the implementation's own expression — so it pins the *implementation*,
not the published definition. A formula wrong in both places passes.

Left as it is: the formula was hand-checked against the erratum and the golden master pins the
resulting number, so the risk is covered from two other directions. Recorded so nobody mistakes
this test for the independent check it resembles — the genuinely independent oracle in the suite
is in `test-generateNoiseImage.R`. If it is ever revisited, the replacement is a worked example
from the paper or a hand-computed 2×2 CI with a known reference vector, not a tidier
restatement of the same expression.

### Pixel assertions have measured the graphics device twice
The practice this produced — uniform backgrounds, colour channels only, comparisons rather
than pinned values — is in `CONTRIBUTING.md`. Here is what it cost to learn, twice in the same
fix. (`plotZmap`'s tests were three-quarters `file.exists()` before any of it; a function
writing a uniformly blank PNG passed them all.)

Count distinct values over the **colour channels only**, never the raw array. Whether a PNG
device writes an alpha channel is a property of the graphics backend, not of what was drawn:
cairo (Linux, Windows) writes RGB, macOS quartz writes RGBA. An opaque alpha plane is a solid
block of 1s, so it adds a second distinct value to the array while nothing has been painted.
The original `expect_length(unique(as.vector(png)), 1)` therefore measured the backend as much
as the drawing, and `test-plotZmap.R` now drops any alpha plane before counting.

It was the *only* assertion in the suite counting distinct values rather than comparing two
renders, which is why it was also the only one that failed when the suite first ran on macOS:
every image-to-image comparison stays green because both sides gain the same plane. Prefer
the comparison form for new tests.

**The absolute pixel value is device-dependent too, and must not be asserted.** The first
attempt at this fix pinned the flat value to the background grey — and failed on macOS, which
renders `bgimage` 0.5 at ~0.573 where cairo gives 0.502. That is colour management, and
precisely the same mistake the alpha fix was correcting, made one line below it. The check is
now an *ordering*: the same render over a darker background must come out darker, which
survives any monotone transfer function.

**The general form:** when a test reads pixels back from a graphics device, every absolute
property of those pixels belongs to the device — channel count and value alike. Only
relationships between renders are portable.

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

### The legacy `.Rdata` fixtures are committed, not generated when the tests run
Both directions of the compatibility promise now have a test. The gate above checks that this
version still computes what old versions computed; `test-legacy-rdata.R` checks that it can
still *read* what they wrote, which no other test can, because every other fixture is built by
the current generator.

Generating those files needs the old version installed — each one builds a cluster whose
workers call `library(rcicr)`, so sourcing its R files is not enough, and v1.0.1 additionally
needs `raster`, dropped in #186. Doing that at test time would put a package install and a
network round trip inside the suite. Instead `tools/make-legacy-rdata.R` installs each tag into
a throwaway library once and the resulting files are committed: 205 KB for 1.0.1 and 45 KB for
1.1.0 at 32px — 1.0.1 is the larger because it is generated at `nscales = 5`, the historical
default, so its missing-`nscales` fallback lands on the right noise basis. The safeguard runs
in every CI job, on every platform, with no network. Regenerating is a deliberate act, as with
the `released-formals` fixtures — a red test here means this version can no longer read a file
a researcher already has, which is the failure, not the fixture.

Writing them also settled the `generator_version` question with data rather than documentation:
the 1.0.1 and 1.1.0 files really do carry a top-level `'0.4.0'` while `p$generator_version`
holds `'1.0.1'` and `'1.1.0'`.

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

### Write paths are required arguments, not defaults of `tempdir()`
CRAN's review of 1.2.1 asked us to "omit any default path in writing functions". Two answers
would satisfy the letter of that: default the paths to `tempdir()`, or remove the defaults.
We removed them.

`tempdir()` was rejected because it is the more dangerous of the two for this package's
users. A researcher whose script relied on `./cis` would keep running, silently writing
classification images into a directory that is deleted when the session ends — the failure
surfaces days later as missing output, with nothing to connect it to an upgrade. Removing
the default turns the same script into an immediate error naming the argument to supply. A
breaking change that announces itself beats a silent relocation of somebody's results.

It also cost nothing to verify: the release gate reports `max|d| = 0` across 135 checks
against v1.2.1, because `tools/compare-harness.R` already passed every path explicitly.

### `captureArgs()` skips required-and-absent arguments, but never defaulted ones
The `load()` guard snapshots a function's arguments and restores them after reading an
`.Rdata` file. Once paths became required, `mget(names(formals()))` started aborting: it
forces the promise, and a wrapper forwarding its own missing argument — `batchGenerateCI()`
passing `targetpath = targetpath` — makes that promise a missing symbol.

The fix must not over-correct. `missing()` reports a *defaulted* argument missing too when
the caller did not supply it, and skipping those would silently reopen the hazard the guard
exists for: a default is the value the function goes on to use, and is exactly as
replaceable by a field in the `.Rdata` as one passed explicitly. Dropping them removed the
`step` guard in `computeCumulativeCICorrelation()` and the test caught it. The predicate is
therefore *required* (no default in `formals()`) **and** absent — not absent alone.

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

### Dependabot watches the actions, not the R packages
`.github/dependabot.yml` covers `github-actions` and nothing else. That is not an oversight,
and the asymmetry is the point.

**The R dependencies cannot be watched, and pretending otherwise is worse than the gap.**
Dependabot has no CRAN ecosystem, and **CRAN publishes no security advisory database** — there
is no `npm audit` equivalent for R. The nearest thing is Sonatype's OSS Index via `oysteR`,
which as of 2026-07-29 returns **401 to anonymous requests** and so needs an account to run at
all. Weighed and rejected: a scanner whose CRAN coverage is thin reports "clean" because its
database is empty, not because the tree is safe, and a green badge that means nothing is worse
than a known gap. The realistic dependency failure — an import getting archived — is already
caught by `R CMD check` on four platforms on every push.

Where the actual CVE surface lives is worth stating, because it is not in R code at all: 29 of
the 53 packages in the recursive tree need compilation, and their vulnerabilities belong to the
C libraries they bind to — **libpng** (`png`), **libjpeg** (`jpeg`), **GDAL/GEOS/PROJ**
(`terra` ← `raster`). Those are patched by the user's operating system, not by CRAN, and
nothing this package does can affect them. Dropping `raster` (issue #186) is the only
lever that shrinks it.

**The actions genuinely needed watching.** Twelve of the thirteen in use were on floating
major tags, so whoever controls those repositories could change what runs in CI at any time —
the one supply-chain surface here with real incident history, and the one Dependabot supports.
Updates are **grouped into a single PR**, because `.github/workflows/` is deliberately not on
the gate's inert allowlist, so ungrouped bumps would each spend a full CI round.

### R-hub runs on `workflow_dispatch` only, never on push
The R-hub v2 workflow is the stock file `rhub::rhub_setup()` writes, kept unmodified so it can
be refreshed from upstream. It is left trigger-on-demand because R-hub answers a question that
only arises at release time — does this build under CRAN's own compiler flags and on the
odd platforms nobody develops on — and running it per-PR would spend a matrix of platforms on
a question nobody asked yet, competing with the ~11-minute reproducibility gate that *is*
required on every code PR.

**R-hub does not stand in for per-platform coverage, and once wrongly appeared to.** This
entry used to justify the per-PR gap by saying the everyday answer was "already covered by
`R-CMD-check.yaml` on release and devel". That matrix varied only the *R version* against a
pinned `ubuntu-latest` — one platform twice. The hole stayed invisible until the first R-hub
dispatch failed a `plotZmap()` test on macOS that had been green on Linux for months (see
the alpha-channel entry under Testing); CRAN builds on macOS, so it would have landed as a
submission ERROR. **The lesson generalises: "already covered by X" is a claim about what X
runs, and worth reading X to check rather than repeating.**

Both platforms CRAN gates on are now in the everyday matrix. Windows was added on different
grounds from macOS, and the distinction is the useful part: macOS had never been run, whereas
Windows had passed on R-hub and win-builder — it was not an untested hole but a *late* one,
checked only at release.

The cost of `workflow_dispatch`-only is that the workflow is invisible until it reaches the
**default branch** — GitHub offers no "Run workflow" button for a file that exists only on a
feature branch. So it merges to `main` ahead of the release it is meant to check.

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

- **A step in `ubuntu-latest (release)`, not a new job or workflow.** Required checks are
  matched by name, and adding a name to the ruleset is a repo-settings write agents here
  cannot make — a new job would report without ever blocking. It runs last so a failure
  cannot cost us the check results.
- **Not a pre-commit hook.** R hooks were kept out of `.pre-commit-config.yaml` for their own
  reasons, but the deciding one here is that the check job already has R and every dependency
  installed, where a hook would install them per run — and an optional local hook is not a
  gate.
- **roxygen2 is pinned to `RoxygenNote`, not `any::`.** `.Rd` output varies between generator
  versions, so `any::` lets a CRAN release of roxygen2 redden a required check on formatting
  alone — an external event breaking the gate, in an environment that cannot re-run a job. The
  step asserts the installed version matches `DESCRIPTION` so the pin cannot silently desync;
  bumping the generator means bumping both.

Nothing rebuilds the *site* locally or in a hook — the pkgdown workflow builds it on every PR
and deploys on push, which is what catches a vignette that no longer knits.

### `CITATION.cff` is generated by `cffr` and compared, not hand-written and field-checked
The first version of this was hand-written, with a script comparing each field to a source of
truth: title and authors to `inst/CITATION`, licence, URLs and addresses to `DESCRIPTION`,
version and release date to `NEWS.md`. It was rejected *after* being built, because reviewing
it turned up five defects in four rounds and every one was the same shape — a field the
comparison did not reach. Unordered URL membership let `repository-code` and `url` swap; author
comparison by name left the address unchecked, on a package CRAN archived for an undeliverable
address; a missing source address compared as `""` against `"(absent)"`, which would have
failed permanently for a future author with no email. A hand-written comparator is a list of
the fields someone remembered, and the failure mode is silence.

`cffr::cff_create()` derives the whole file, so the comparison is structural — regenerate,
compare every key, and nothing depends on remembering a field. It also validates against the
CFF schema, which is what makes a malformed file fail loudly; GitHub declines to render the
"Cite this repository" button on a file it cannot parse and says nothing.

Two measurements shaped how it is generated. `dependencies = TRUE`, the default, emits 380
lines describing every dependency's authors, with years read from whichever versions happen to
be installed — regenerating on a different machine produces a different file, so it is off.
And `preferred-citation`'s year comes from `inst/CITATION`, which reads the clock when
`DESCRIPTION` carries no publication date, so that one field changes on 1 January with nothing
else: the comparison excludes it, or a required check would go red on a calendar boundary.
`gh_keywords = FALSE` keeps generation off the network.

The dependency and the second version pin — the two costs that argued for hand-writing it —
are real but bought something. `version` now tracks `DESCRIPTION`, including the `.9000`
suffix between releases, where the hand-written file deliberately named the last release; that
is the one thing the switch gave up.

---

## Releases and git

The conventions and their rationale are in `AGENTS.md` — trunk-based with tags and no
`develop` branch, `.9000` on `main` between releases, tarballs built from the tag, squash
merges, delete merged branches, `NEWS.md` ordered largest-impact first. Two things learned
the hard way that those entries do not cover:

- **Squash merging is what makes a bad commit cheap.** One branch here carried ~75,000 lines
  of committed `R CMD check` artifacts across three commits; `main` only ever saw the final
  tree, so no history rewrite was needed. (Both patterns are now in `.gitignore` and
  `.Rbuildignore`, and the `git diff --stat main...HEAD` habit is in `CONTRIBUTING.md`.)
- `git push origin --delete a b` is **atomic**: one stale name fails the whole command and
  takes the other branch down with it. (`gh pr merge --delete-branch` has already removed the
  remote branch anyway.)

### External checks run on the release branch; the tag and `.9000` still precede acceptance
The procedure is in `CONTRIBUTING.md` → "Releasing", step 2. What is not obvious is why it
can work at all.

**The circularity that appears to force checks after the tag is not real.** It looks as
though results cannot be recorded in the commit they describe, because writing them changes
the tree and so invalidates the tarball that was checked. But `cran-comments.md` is
`.Rbuildignore`d and verifiably absent from the built tarball, so editing it produces the
same tarball byte for byte. It is not a package artifact at all — it is the text pasted into
the "Optional comment" field of the CRAN submission form. At 1.2.1 the checks ran *after* the
tag and the `.9000` reopen, so the results landed on a development tree and the tag recorded
no evidence that the tree it names had passed anything.

**Deliberately *not* copied from `usethis:::release_checklist()`: tagging and reopening
`.9000` still happen before CRAN accepts.** The `.9000` suffix is what selects `--quick` over
the full ~20-minute battery, so holding `main` at a clean version across the weeks CRAN can
take would run the full gate on every unrelated PR in that window. And a rejection means
shipping X.Y.Z+1 anyway, leaving the tag naming a tree that was never on CRAN — already true
of `v1.2.0` and `v1.2.2`, so the repository tolerates that cheaply.

This variant is in one respect tighter than usethis's, which runs win-builder while
`DESCRIPTION` still carries `.9000` and never re-checks after `use_version()` — recording
results for a different version string than the one submitted.

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

### Answer CRAN from the review text, and never send them a question a sweep can answer
CRAN's review of 1.2.1 listed two `.Rd` files under "some code lines in examples are
commented out". Working from a summary that had kept one of them, 1.2.2 replied that we could
not find the commented line and **asked the reviewer which line she meant** — while the fix
for the file she had named correctly sat in the same commit, made under a different point.
Two further drafts repeated it.

Three rules came out of that, and the first two are in `AGENTS.md`: log every reply verbatim
in `notes/cran-review-<version>.md` and answer from the file; re-verify every claim in
`cran-comments.md` against the tree rather than carrying it forward. The third is a matter of
tone: **a question to the reviewer costs a round trip measured in weeks, and a sweep costs
minutes.** Sweep everything the point could apply to, then report what was found.

The same instinct applies in reverse — **do not explain a note the reviewer does not have.**
`cran-comments.md` carried a paragraph answering a `medium.com` 403 that appears only on our
local check, on no external check for 1.2.3, and not in CRAN's own pretest of 1.2.1. It was
removed; the link stays in `README.md`, because `README.md` ships and removing it would
invalidate the built tarball and force all five external checks to re-run for something no
CRAN-side check has ever raised (issue #192).

### Claims must survive on someone else's machine
`cran-comments.md` twice stated a total runtime for the example set — "about nine seconds",
then "about fifteen" — measured here. Both win-builder runs report 18s, so the reviewer's own
log would have contradicted the letter arguing for the package. The rule is in `AGENTS.md`;
what it is *not* is a ban on numbers. A ratio survives, and so does a comparison to a fixed
bar such as CRAN's five-second per-example limit; `NEWS.md`'s "about 6x faster, 1.66s to
0.28s per call" earns its place because a user feels that difference and the ratio holds
anywhere. A bare absolute describes our hardware to a reader who has their own.

### A work tracker carries no project state
`BACKLOG.md` once held a "state as of" block plus a "Last updated" narrative — some 78 lines
describing where the release stood. Every fact in them lived somewhere with a better claim:
`NEWS.md`, `cran-comments.md`, `CONTRIBUTING.md` → "Releasing", and the open PRs. **It drifted
four times, each within a day of a release**, which is what a duplicated fact does rather than
what a careless author does; the fourth fix was to write it more carefully, and it was wrong
within the hour. The rule outlived the file and applies to the issue tracker: a *hold
condition* belongs on the issue it holds; the project's current position belongs nowhere in
the tracker at all.

Same reasoning as deleting the hand-maintained `Version:`/`Date:` table from
`man/rcicr-package.Rd` in 1.2.3: the way to keep a fact current is to stop writing it twice.

### The backlog moved to GitHub Issues, after the stale issues were triaged
Done 2026-08-08: 21 issues opened as #174–#194 under `P0`–`P3` labels, and `BACKLOG.md`
deleted. The sweep came first and was the hard prerequisite.

**The tracker is treated as a working surface, not a curated public product surface** — so
internal maintenance work belongs in it alongside user-visible bugs, as in r-lib repos. The
alternative considered was keeping the tracker for things a researcher would recognise as
affecting them and holding chores elsewhere; rejected because two backlogs is the same
duplication problem one layer up, and because the drift this file produced came precisely
from status living somewhere a human has to remember to update.

**The ordering was not incidental.** All 25 open issues dated from 2016–2017, and opening 21
new ones into that would have buried the triage under the migration. It also matters which half of
this carries the value: **the sweep alone fixes the thing that actually costs something** —
a tracker that reads as abandoned in 2017 while the package is being resubmitted to CRAN.
The migration is the smaller remaining gain, which is why it is sequenced second rather than
bundled.

What deliberately does **not** become an issue: the "already correct — do not re-fix" entries.
Those are decisions and belong here, and four moved into this file rather than onto the
tracker — the InfoVal formula, `return_as_dataframe`'s one-image-per-trial shape, the emptied
`ref_lookup` rows, and the infoVal test oracle that mirrors its implementation. The formula one
was *already* duplicated here, which is what made the misfiling visible in the first place.

**Known cost, accepted:** changes to the plan stop going through a reviewed PR, because issue
edits are not reviewable. For a single maintainer that is close to zero.

### `CLAUDE.md` is a stub that imports `AGENTS.md`, not a symlink
Claude Code reads `CLAUDE.md` and does not read `AGENTS.md`, so renaming the conventions file
in #166 silently stopped it loading in every agent session until 2026-08-08 — including the
rules on `NEWS.md` headings, `.Rbuildignore` and `test-fixed-bugs.R`, each of which costs a
release when broken.

The documented fixes are a `CLAUDE.md` containing `@AGENTS.md`, or `ln -s AGENTS.md
CLAUDE.md`. **The symlink was rejected**: it requires Administrator privileges or Developer
Mode on Windows, and this package has Windows contributors and a Windows CI runner. The
import is a plain file that git and every platform treat identically, and it leaves room for
Claude-specific instructions below it. `AGENTS.md` stays the single source of truth.

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

### Three vignette figures were wrong in ways only viewing them showed
The rule this stands behind is in `CONTRIBUTING.md`. What it caught, none of it visible to an
assertion:

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
It had been documented since 2016 and never worked — the import half landed with an explicit
"todo: applying the mask" that was never picked up. Filed as a behaviour change, but the risk
that framing guards against could not materialise: nothing has ever been masked, so no
published z-map depends on the current output. Deprecating would have removed a documented
feature on the grounds that it was never built.

**The generalisable part:** "this changes rendered output" is not by itself an argument for
leaving something broken. Ask who is currently *relying* on the broken output. Here the answer
was nobody, and the entry sat open a day longer than it needed to.

### The `.Rdata` anatomy belongs in `README.md`
It is the only link between the two halves of the package, and nothing about a stimulus set is
recoverable without it. The field-by-field table there was written by inspecting a real
generated file rather than by reading the `save()` call — which is how `trial` was identified
as a leftover loop counter carrying no information.
