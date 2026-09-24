# Decisions

Why `rcicr` behaves as it does: the measurement that ruled an option out, the alternative that looked obvious and was wrong, the thing that looks like a bug and is not.

**This file is about the package, not about the repository**: a decision here would still matter if `rcicr` were maintained somewhere else entirely. `AGENTS.md` → "Which file a thing goes in" says what that excludes, when an entry is worth adding, and where other material goes. Entries are grouped by theme and edited in place.

**Keep this file under 5200 words**; over budget, something comes out before something goes in. Write the decision and the evidence, not the route to it: an entry earns its length from a measurement or a rejected alternative.

> Until 2026-07-27 this was a chronological session log (`.session-log.md`). The original narrative, with dates and intermediate states, is in git history up to `887aea4`.

---

## The constraint that shapes everything

**Researchers re-run old analysis scripts years later and publish what comes out.** Every other decision here follows from that. The rules it implies, and the golden master that enforces them, are in `CONTRIBUTING.md`.

---

## Numerics and the random number stream

### `purrr::rbernoulli()` was replaced with `runif()`, not `rbinom()`
`purrr::rbernoulli(n, p)` computes `runif(n) > (1 - p)`. `rbinom()` consumes the seeded stream differently (verified across 150 seed and probability combinations) and would change every derived InfoVal. The replacement therefore uses the bit-identical `runif()` expression.

### A base image's alpha channel is discarded, not composited
Greyscale conversion drops alpha and uses the stored colour channels. Compositing would invent a background value and make it part of the contract; rejecting transparency would break files the package has always accepted. `rcicr` renders opaque stimuli, so alpha means nothing downstream. A cut-out may therefore reveal colour that a viewer hid under transparency.

### `rowMeans(x, dims = 2)` was adopted despite not being bit-identical
Patch averaging in `generateNoiseImage()` moved from `apply(..., 1:2, mean)`, about 6x faster end to end. The two sum in a different order and so differ by about 1 ULP (~1e-19 on pixel values around 0.01). Adopted because an **independent oracle**, the average as an explicit triple loop using neither function, put both within ~5.6e-17 across noise types, scales and seeds. At the golden master's configuration they are bit-identical.

The trap: **`rowMeans()` on a 3-D array defaults to `dims = 1`**, collapsing dimensions 2 *and* 3, and `array()` then silently recycles the short result. The first version submitted did exactly that, with a measured maximum deviation of 0.21 on data with an SD of ~0.01, not the ~1e-17 it claimed.

### The stimulus seed's stream is load-bearing well beyond stimulus generation
`generateStimuli2IFC()` seeds with the stimulus seed and draws one `runif()` value per parameter per trial. The simulated reference responses continue that stream, so the null depends on the file and the session's `RNGkind()`, not on the current random state or `ncores`. Files do not record the RNG kind, and changing it changes the response draws but not the saved noise; [#315](https://github.com/rdotsch/rcicr/issues/315) tracks this. `?generateReferenceDistribution2IFC` documents the historical stream, and tests pin it.

The stimuli are no longer regenerated to reach that point in the stream ([#301](https://github.com/rdotsch/rcicr/issues/301)), so `seedResponseStream()` **replays** their draws. Do not "simplify" that loop away: it looks inert and is the whole guarantee. Its draw count is the saved matrix's width *after* `selectStimulusParams()`; a pre-0.3.0 file's raw 4096 columns would shift the stream, as below.

### `response_seed`, not `seed`
Four choices, each forced by something specific:

- **The name.** `seed` is an *object in the `.Rdata` file*, and `generateReferenceDistribution2IFC()` re-saves its own frame. An argument called `seed` would overwrite the stimulus seed and write it back, corrupting the record of how the stimuli were made.
- **It seeds the responses**, replacing the replayed stimulus stream. Forwarding it into stimulus generation would describe stimuli participants never saw.
- **`NULL` replays the stimulus stream instead of seeding afresh**, so the default path stays byte-identical, not merely equivalent-looking.
- **`computeInfoVal2IFC(response_seed=)` forces a new reference and never stores it.** Without forcing, the seed would be silently ignored on every call after the first, since the generator writes `reference_norms` into the file. Storing a one-off Monte Carlo check would redefine what every later InfoVal from that stimulus set means.

### 4096 → 4092 parameters: why old stimulus files cannot be regenerated from their seed
rcicr 0.3.0 (2015-01-23) cut the random draws per trial from 4096 to 4092. 4092 is the real patch count, 6 orientations × 2 phases × `sum(4^0..4^4)`; 4096 was a round `2^12` over-allocation, so four contrasts per trial were drawn that no patch index used. `ChangeLog` says the change "does not affect anything else".

**That is true for analysis and false for regeneration.** Analysing a pre-0.3.0 file reads its *stored* parameters, and four unused columns change nothing. But `runif()` is a stream: at the same seed, trial 1 gets identical values either way and **every later trial is shifted by four draws** (verified: trial 1 identical, trials 2 and 3 not). A pre-0.3.0 stimulus set can be re-analysed exactly, but not re-created from its seed.

The same release fixed `sinIdx` counting from 0 instead of 1, an independent change that the `pre_0.3.0` flag handles.

**The 0-based path in `generateNoiseImage()` looks as if it corrupts the whole image, and does not.** When `min(patchIdx) == 0`, `params[p$patchIdx]` drops every cell indexed 0 (R drops a 0 subscript rather than returning `NA`). The result is too short for the patch array, and `array()` recycles it, which would normally misalign everything. It is harmless because the same counting from 0 leaves the *last* patch layer unwritten, so every cell with `patchIdx == 0` also has `patches == 0`. Those cells come last in column-major order, so the recycled values land only where the patch is zero and are multiplied away. Measured: identical to the honest "one patch not shown" result (maximum absolute difference 0) across 36 combinations of size, `nscales` and seed.

A stimulus file stores the patch array its own generator wrote, and the genuine pre-0.3.0 generator is the same `co = 0` / `idx = 0` loop that the `pre_0.3.0` flag runs. That was verified against the R-Forge source (`git show 7d0d9e6:pkg/R/rcicr.R`): 0.3.0 only flipped the default and added the flag. Do **not** "fix" the recycling by offsetting the index by one: that changes which sinusoid is dropped and alters the CI of every genuine pre-0.3.0 file. `test-generateNoiseImage.R` pins both properties.

**A backward-compatibility path that nothing exercises cannot be told apart from one that works.** The truncation left one broken in `generateCI()` for eleven years, untested until 2026-07-28; so was the `sinusoids`/`sinIdx` path.

### `load()` assigns into the calling frame — check every new argument against saved names
An object in an `.Rdata` file silently overwrites a function argument of the same name. `generateReferenceDistribution2IFC()` re-saved its whole frame, so the files it wrote contained `rdata` and `ncores`. A second call then ignored the caller's `ncores` and wrote back to the path recorded by the first. It is fixed at the source, by leaving the function's own arguments out of the save, *and* defensively on read, for files older versions already wrote.

### The InfoVal formula is already correct — do not "fix" it
It is equations 2 and 3 of Brinkman et al. (2019, *Behavior Research Methods* 51, 2059-2073): `(norm - median) / (k * MAD)`, with `k = 1.4826`. Schmitz, Rougier and Yzerbyt (2019, <https://doi.org/10.3758/s13428-019-01295-1>) reported two miscomputations in rcicr: the one norm instead of the Euclidean norm, and a missing `k`. Their erratum (Schmitz et al., 2020, *Behavior Research Methods* 52, 1800-1801, <https://doi.org/10.3758/s13428-020-01367-7>) withdraws the second, because R's `mad()` already applies `constant = 1.4826`. The code uses the Euclidean norm (`norm(x, "f")`) and plain `mad()`. **Do not add a `k`**: it would be applied twice. A regression test pins the result.

---

## Things that look like bugs and are not

### `autoscale()` leaves `$combined` untouched — intentional
At Ron's direction, `$combined` stays as the caller passed it, and `$scaled` holds the autoscaled result. After `batchGenerateCI()`, which scales with `'none'` before autoscaling, `$combined` overlays unscaled noise and can look almost blank. `(ci$scaled + ci$base) / 2` gives what `save_as_pngs = TRUE` writes. Changing `$combined` would change the images existing scripts plot.

### A CI with no range renders neutral, not NaN — `matched` differently again
Exactly cancelling responses give an all-zero CI. That is a result, so the range-based methods fill it with a neutral value instead of dividing by a zero range into NaN ([#303](https://github.com/rdotsch/rcicr/issues/303)). `0.5` is a **choice, not a limit**. The constant is derived from the CI, so along `ci = t * x` the scaled result does not depend on `t`: an all-positive pattern gives 1, an all-negative one 0. `0.5` is what `constant` scaling already returns here, and what `autoscale()` gives a zero CI next to one with signal.

`matched` uses the midpoint of the **base image's** range instead, because it renders into that range and `0.5` can fall outside it: with a base spanning `[0, 0.3]`, a no-signal CI would be brighter than any pixel of the base. With the default contrast maximization the base spans `[0, 1]`, so the two agree.

The triggers differ on purpose. `matched` needs a range, so a uniform *non-zero* CI triggers it too; `independent` only an exactly zero one. A fully masked CI has no values rather than no range, and is left alone. Raising an error was rejected: one cancelling participant would abort a whole `save_individual_cis = TRUE` call.

### `base_face_files` validation rejects two inputs that used to run
Duplicate names, and a list element that is not a single file name, now stop the call. Both used to run, so this is a considered exception to the constraint above, allowed because neither ever produced what the call asked for. `list(face = 'a.png', face = 'b.png')` looks each name up as a string, so it wrote **one** stimulus set, from `a.png` (verified against the pre-fix code: four files for two trials). A script that "worked" this way was generating stimuli for a base image it never named. Every other new check rejects input that already failed further on, inside a parallel worker.

**Readability is checked by attempting the read, not with `file.access()`**, which is documented as unreliable on Windows and network filesystems and can reject a file that reads fine. Keeping the reader's own error message also preserves the difference between "not a PNG" and "permission denied".

### `_R_CHECK_LIMIT_CORES_` handling caps cores under check only
`default_ncores()` returns 2 when that variable is set, and `detectCores() - 1` otherwise. Only the checker sets it, so CRAN's two-core policy is met and behaviour at the console is unchanged. The known downside, that `detectCores()` over-subscribes under cgroup limits, is accepted.

### `generator_version` in old `.Rdata` files is not trustworthy
It was hardcoded as `'0.4.0'` from 2016 until 1.1.0.9000, so every file written in that range claims to be 0.4.0; `p$generator_version` held the real value all along. Code reading the field must treat `'0.4.0'` as "unknown, somewhere in that range", accept both a character string and a `package_version`, and compare with `numeric_version()`: as text, `'0.10.0' < '0.4.0'` is `TRUE`. The `pre_0.3.0` compatibility path does **not** depend on this field; it detects the old `sinusoids`/`sinIdx` layout from the file's structure.

### `return_as_dataframe = TRUE` returns one noise image per trial, not per trial × base image
The frame has one column per trial, so it cannot hold trial × base image and carries the first base image's noise. Under the default `use_same_parameters = TRUE`, that is every base image's noise. This is documented on `@param return_as_dataframe` rather than changed: widening it would change the return shape, which needs a **new argument**, never a redefinition. Stimuli are written for every base image either way; that they once were not was [#302](https://github.com/rdotsch/rcicr/issues/302), a bug, not this decision.

### Independent base images are scored against their own reference
With `use_same_parameters = FALSE`, base images after the first were scored against the first one's null ([#299](https://github.com/rdotsch/rcicr/issues/299)), at a cost measured in [`analyses/infoval-reference-impact.md`](analyses/infoval-reference-impact.md). Independent-base references now take a `baseimage` label and are stored per base. Shared-parameter numbers, and the first base's numbers from 0.3.0 on, are unchanged (maximum absolute difference 0), because the first base's parameters come from the same leading block of the random stream. A *pre-0.3.0* independent file is the exception: its trials cannot be rebuilt from the seed (see "4096 → 4092" above), so its first base moves too.

### `computeCumulativeCICorrelation()` does not aggregate repeated stimuli, and its curve ends at 1 by construction
`generateCI()` averages the responses to each unique stimulus before building its CI (`aggregateResponses()` in `R/ci-inputs.R`). `computeCumulativeCICorrelation()` does not: it takes trials in presentation order, and averaging repeats would discard the order a cumulative curve is about.

Without a `targetci`, the final CI is built from the same unaveraged trials as the curve, so the curve **ends at exactly 1**: self-consistency, not convergence. That holds only when the evaluated trials reach the last one (measured: six responses at `step = 2` stop at the fifth and end at 0.967; at `step = 3`, at 0.938). It does not hold at all when responses cancel exactly: that CI is zero everywhere, and a correlation with a constant is `NA` at *every* point.

That final CI equals `generateCI()`'s only when every stimulus was shown equally often, and then they are bit-identical. Unequal counts weight the data differently, each trial equally here and each unique stimulus equally there: counts 3/1 correlate at 0.845, counts 4/2/1/1 at 0.773.

Pinned by a test rather than changed: averaging the final CI would change the numeric output for anyone calling without `targetci`, *and* stop the curve ending at 1. Whether `generateCI()`'s own weighting suits unbalanced designs is a separate issue.

### A cached reference is trusted only on positive evidence
Old references were built by reconstructing the stimuli, using settings that could be missing and the session's RNG kind. So a default-stream reference without a provenance marker is rebuilt once; a seeded one is kept, as a deliberate choice. One resolver handles both the shared and the per-base layout. An automatic rebuild keeps the stored iteration count and the caller's random state, and a read-only file is scored in memory.

Older rcicr versions keep unknown fields while replacing the norms, so the marker is tied to a full copy of the norms, checked with `identical()`. That covers every value with no hashing, formatting or arithmetic. It costs eight bytes per double (about 80 KB for 10,000 norms), small next to the noise basis. Checking a sample could miss a replacement, and a hashing package would add a compiled dependency.

A changed reference is reported with a message, not a warning, so that under `warn = 2` reporting the correction cannot stop the corrected value being returned. The cost: `warnings()` does not collect it and `suppressMessages()` hides it, which `NEWS.md` states. Iteration counts the caller chose still get the reliability warning.

### Repopulating `ref_lookup` costs four measurements — and the two halves stand or fall together
`AGENTS.md` says what the table is; this entry is the way out, because doing either half alone is worse than the present state. Repopulating means four numbers: `median(reference_norms)` and `mad(reference_norms)` under the current formula, for seed 1, 512px and 10000 iterations, at 100, 300, 500 and 1000 trials. The alternative is deleting the ~55 lines of matching and prompt code. **Do not do half of either**: delete the code while meaning to re-measure, and the feature becomes unrecoverable rather than dormant.

---

## Testing

### Vacuous assertions have shipped here twice, and both looked fine on the page
Two `batchGenerateCI*` tests asserted only length, names and `dim`, so **grouping was never checked**, and a bug feeding all trials to every group passed. Another compared `.Rdata` fields that its `save_rdata = FALSE` call never writes. Mutation testing caught both; reading the tests did not. Every grouping and threshold test now also asserts that the *wrong* answer differs.

**Before pinning a summary statistic, ask what the transformation already guarantees.** A `zmapmethod = "quick"` z-map ends in `scale()`, so its mean of 0 and SD of 1 hold regardless: changing `sigma` left both bit-identical while every cell moved. They are still asserted, as a check that the standardization happened, alongside statistics that do vary.

### `computeInfoVal2IFC`'s test oracle mirrors the implementation, and is kept anyway
`test-computeInfoVal2IFC.R` recomputes `(norm(ci, "f") - median(reference_norms)) / mad(reference_norms)`, the implementation's own expression. So it pins the *implementation*, not the published definition: a formula wrong in both places passes.

It stays: the formula was checked against the papers (see "The InfoVal formula" above), and the golden master pins the resulting number. The suite's genuinely independent oracle is in `test-generateNoiseImage.R`. A real replacement is a hand-computed 2×2 CI with a known reference vector.

### Pixel assertions have measured the graphics device twice
The practice is in `CONTRIBUTING.md`; both measurements behind it came from one fix.

**The channel count belongs to the backend.** cairo (Linux, Windows) writes RGB where macOS quartz writes RGBA, so counting distinct values measured the backend as much as the drawing; it was the only assertion that failed when the suite first ran on macOS. **So does the absolute value**: quartz renders a 0.5 background at ~0.573 where cairo gives 0.502. The check is now an *ordering* (the same render over a darker background must come out darker), which survives any monotone transformation.

### The recovery test uses a permutation null, not a parametric one
`test-recovery.R` gives a simulated observer a known template and asserts that `generateCI()` recovers it, scored by `cor(vec(CI), vec(template))`. The null permutes response labels across trials. That keeps the noise images and the balance of 1 and -1, and destroys only the pairing of responses with stimuli.

**Why not a t-test on Pearson's r:** it would use df = n_pixels − 2 = 1022 and be far too liberal. The basis has 60 patches, so a 32×32 CI has an effective dimensionality of about 60, not 1024, and neighbouring pixels are correlated by construction; measured null correlations reach |r| = 0.454. Every null CI is built from the same basis, so the null has exactly that autocorrelation.

Sensitivity, established by mutation: correct pairing 0.71, sign flipped −0.74, responses reversed −0.18, shifted by one trial 0.24. **A known limitation:** it cannot catch an error applied consistently inside `generateNoiseImage()`, because the template is built through the same function and the error cancels. The oracle test covers that.

### The release gate runs the old code; the golden master only re-runs ours
`test-regression-baseline.R` pins values *this repository computed for itself*, so it catches only drift after the numbers were recorded. Had a P0 fix changed results before the baseline was recorded, the baseline would have pinned the changed values and passed forever. `tools/compare-release-output.R` closes that gap: it installs the reference commit into a temporary library and runs both versions over the same battery, the only place here that runs the old code. The two complement each other. The golden master runs on every commit; the gate costs two package installs and minutes of compute, so it runs `--quick` on PRs and in full at release.

The battery is copied into the temporary directory before either side runs, so editing the working copy mid-run cannot make the two sides compare different things.

### The legacy `.Rdata` fixtures are committed, not generated when the tests run
Both directions of the compatibility promise have a test. The gate checks that this version still computes what old versions computed; `test-legacy-rdata.R` checks that it can still *read* what they wrote. No other test can, because every other fixture is made by the current generator.

Generating the fixtures needs the old version *installed*: each builds a cluster whose workers call `library(rcicr)`, so sourcing its R files is not enough, and v1.0.1 also needs `raster`, dropped in #186. Doing that at test time would put a package install and a network round trip inside the suite. Instead, `tools/make-legacy-rdata.R` installs each tag into a throwaway library once, and the files are committed, so the check runs in every CI job, on every platform, with no network. A red test here means this version can no longer read a file a researcher already has; do not regenerate the fixture to fix it.

Each fixture is generated at its era's **defaults** (`nscales = 5`, `sigma = 25` for 1.0.1), the situation a returning researcher is actually in. The gabor fixture is the only one whose saved basis is not sinusoidal. It was written to exercise the fallbacks for missing fields, which [#301](https://github.com/rdotsch/rcicr/issues/301) removed by reading the saved basis instead of rebuilding it; it now checks that basis being read back.

### The v1.0.1 reference is pinned; the previous release is a *second* run, not a replacement
The obvious move once a release is green is to make it the new reference. It is wrong. Each release would then be compared only with its predecessor, and the code could walk away from the published numbers one tolerated epsilon at a time, every step "identical to the last release". The literature was produced with v1.0.1, so that comparison is the one that protects it, and it stays pinned at `v1.0.1` (tagged retroactively at `b6ab269`, so the default reads as a version rather than a bare SHA).

The second run, against the previous release, answers a different and also useful question: did anything break since then? It also reaches further, because v1.0.1 *crashes* on calls that later versions return numbers for. That is why `EXPECTED` entries name the reference they apply to: a deviation from v1.0.1 is not a deviation from the previous release, and an entry that fired for one would be reported as stale by the other.

### The battery stops where the reference version crashes
Measured on 2026-07-28 against v1.0.1 on R 4.3.3: it can produce a z-map **only** at 512px with `zmapdecoration = TRUE`. Undecorated, it dies in `if (bgimage != '')` ("the condition has length > 1"); at 64 and 128px it dies in `plot.new()` ("figure margins too large"), and at 64px decorated even earlier, on the `plt` graphical parameter. `mask` fails for the same kind of reason. All are fixed here, and none can be compared with v1.0.1: **a fix that turns a crash into a number has no old number to compare with.** The test suite covers those paths, and the gate does from v1.1.0 on; the `SINCE` table in `tools/compare-harness.R` records which extras need which reference.

### Tolerances: 8 ULP scaled to the values, plus an 8-bit pixel count
A flat `.Machine$double.eps` suits a classification image (values around 0.01) but is far too tight for a z-map, whose values reach several units and whose ULP steps are ~4× larger. So the tolerance is `8 * eps * max(1, max(abs(reference)))`: a few ULP *of the values involved*. Anything that could only change through a different random stream, algorithm or file format (patch indices, drawn parameters, the base image, the stimulus PNGs) must be bit-identical instead.

**Why 8 and not 1.** These are not single passes over the data. A z-map is a convolution over 262,144 pixels followed by a standardization over the same values, so a 3.5e-18 difference in the CI arrives as **1.33e-15** (measured at 512px against v1.0.1). One ULP of the largest value (9.0e-16) rejects that; eight accepts it and is still three orders of magnitude below anything visible.

What protects this comparison is not the tolerance but the two **exact** checks beside it: after quantizing to 8 bits, 0 of N pixels may differ, and the NA pattern must match cell for cell. The z-map sigma bug moved 1,282 cells across the threshold and was caught by the NA check; every numeric tolerance considered here would have let it through. The 8-bit check also answers a researcher's actual question: does the PNG they publish change?

---

## Performance and parallelism

### `ncores == 1` runs in-process instead of building a one-worker cluster
`startBackend()` in `parallel.R` registers `doSEQ` when `ncores < 2`, so the same `%dopar%` loops run in the current process and **no loop body changed**. The test suite went from 140s to 4s; under `R CMD check`, from `[8s/126s]` to `[8s/37s]`. Eight seconds of CPU against 126 elapsed was worker startup, not computation.

**This is safe because neither parallel loop draws random numbers**, so there is no per-worker random stream to diverge. Verified: `test-parallel-equivalence.R` pins `ncores = 1` and `ncores = 2` to bit-identical output and spells out what would make this unsafe again.

### Parallelism stays on `parallel` + `foreach`
Not `future`, `snowfall` or MPI. Issue #66 asked for snowfall as a single-core fallback; `registerDoSEQ()` provides that without a dependency. Issue #63 asked for MPI and `doSNOW`; MPI needs the script launched under `mpirun` instead of a plain R session, a usability cost no speedup here justifies.

`doSNOW` replaced `doParallel` for #178, because only it honours `.options.snow`. Its `progress` callback runs in the **parent**, whereas a progress bar ticked inside a `%dopar%` body only moves a worker's private copy. `registerDoSEQ()` ignores the callback (measured), so the serial path keeps its in-body ticks behind an `is.null(cl)` guard. `doParallel` had no other use, so it was removed.

### Memory: issue #12 was solved by deletion, not by chunking
The issue proposed spreading stimulus matrices over several `.RData` files. The real cause was a preallocated `zeros(img_size, img_size, n_trials)` array, about 1.5 GB at the defaults, living in the parent environment. `foreach` therefore exported a full copy to *every* worker, each of which wrote one slice and discarded the rest. Each iteration now allocates only its own trial's noise.

---

## Arguments and internal guards

### Trial alignment errors stop computation

Recycled participant IDs and malformed stimulus indices can produce plausible CIs from the wrong trials (#294, #300). A warning would still return those results, so these calls raise an error instead. Factor and character IDs need explicit conversion based on the experiment records: factor codes and matrix row names are not stimulus numbers. Validation comes before averaging, because averaging can drop missing IDs. The all-NA convention for "no participants" stays; what a partly missing participant vector should mean is a separate question.

### Write paths are required arguments, not defaults of `tempdir()`
CRAN's review of 1.2.1 asked us to "omit any default path in writing functions". A `tempdir()` default would satisfy the letter of that, and was **rejected as the more dangerous option**. A script relying on `./cis` would keep running while silently writing its classification images into a directory deleted at the end of the session, surfacing days later as missing output with nothing linking it to an upgrade. Without a default, the same script stops immediately with an error naming the argument.

It cost nothing to verify: the gate reports `max|d| = 0` across 135 checks, because `tools/compare-harness.R` already passed every path explicitly.

### A stimulus file is saved in place behind a backup, not replaced by a rename
Renaming a copy over the original (#333) was **rejected**: it changes the file's owner, group and ACLs, which base R cannot restore. Where it can, `saveRdataSafely()` keeps a verified backup until `save()` completes. Barring concurrent saves, rcicr never deletes or overwrites a `.rcicr-backup` it cannot prove it made; restoring needs a check first.

### `captureArgs()` skips required-and-absent arguments, but never defaulted ones
The `load()` guard copies a function's arguments and restores them after reading an `.Rdata` file. Once paths became required, `mget(names(formals()))` started failing: it forces each promise, and a wrapper passing on its own missing argument (`batchGenerateCI()` passing `targetpath = targetpath`) turns that promise into a missing symbol. The guard is used by `computeInfoVal2IFC()` and `computeCumulativeCICorrelation()`; `generateCI()` loads into a private environment instead, so it has no frame to guard.

`missing()` also reports a *defaulted* argument as missing when the caller did not supply it, and skipping those would reopen the hazard: a default is the value the function uses, and an `.Rdata` field can replace it just as easily. Skipping them removed the `step` guard in `computeCumulativeCICorrelation()`, and a test caught it. So the condition is *required* (no default in `formals()`) **and** absent.

## Documentation

### The Medium walkthrough moved into a vignette, and the post stays up
A tutorial outside the repository cannot run at build time, so it goes stale unnoticed, as this one had: `saveasjpegs` is now `save_as_pngs`, and `install_github(..., ref = "development")` names a branch that does not exist. The post stays up for its inbound links; the vignette is the canonical copy.

### Three vignette figures were wrong in ways only viewing them showed
Viewing caught what assertions missed:

- `image()` stretched each rescaling to look identical, until `zlim = c(0, 1)` fixed the palette.
- The quick z-map's scores, relative to the image's own pixels, spanned only ±1.65 at small sizes, so the default threshold returned a misleading blank map.
- A white-noise base image obscured every figure; a Gaussian blob stays legible and is still synthetic.

### When the docs and the code disagreed about `mask`, the code was the contract
Both `?plotZmap` and `?generateCI` said a *matrix* masks where the cell is `1`/`TRUE`, while a *PNG* masks where it is black (`0`): two opposite conventions in one sentence. `applyMask()` has always done `mask_matrix == 0` regardless, so the matrix half was simply wrong.

**The documentation was corrected to match the code, not the reverse**, because users build masks from observed behaviour: anyone whose masked CI came out right has a mask with `0` = masked, whatever the page said. Changing the code to match the prose would have silently inverted every existing mask, the one outcome with a real cost, since a mask covers a face and an inverted one looks plausible.

**The first fix got this backwards.** It believed the docs for `plotZmap` and branched on the mask's origin, masking `0` for a PNG and `TRUE` for a matrix, so the same mask would have removed opposite halves in the two functions. The answer was in the sibling implementation, not in any documentation. **When two documented conventions conflict, run the one that has been executing for a decade.**

### `plotZmap(mask = ...)` was applied rather than deprecated
It had been documented since 2016 and never worked: the import half landed with an explicit "todo: applying the mask" that nobody picked up. It is filed as a behaviour change, but the risk that label guards against could not occur: nothing had ever been masked, so no published z-map depends on the old output. Deprecating it would have removed a documented feature because it was never built.

"This changes rendered output" is not on its own a reason to leave something broken; ask who *relies* on the broken output.

### The `.Rdata` anatomy belongs in `README.md`
The file is the only link between the two halves of the package, and nothing about a stimulus set can be recovered without it. The field-by-field table there was written by inspecting a real generated file rather than by reading the `save()` call, which is how `trial` was identified as a leftover loop counter carrying no information.
