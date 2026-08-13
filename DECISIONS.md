# Decisions

Why `rcicr` behaves as it does: the measurement that ruled an option
out, the alternative that looked obvious and was wrong, the thing that
looks like a bug and is not.

**This file is about the package, not about the repository** — a
decision that would still matter if `rcicr` were maintained somewhere
else entirely. What that excludes, when an entry is worth adding, and
where the other material goes are in `AGENTS.md` → “Which file a thing
goes in”; entries here are grouped by theme and edited in place.

**Keep this file under 5200 words.** It is read to answer “why is this
like this”, and a file long enough to skim is one whose answer is never
found. Over budget, something comes out before something goes in. Write
the decision and the evidence, not the route taken to it: an entry earns
its length from a measurement or a rejected alternative, never from
narrating how it was reached.

> This was a chronological session log (`.session-log.md`) until
> 2026-07-27. The original narrative, with dates and intermediate
> states, is in git history up to `887aea4`.

------------------------------------------------------------------------

## The constraint that shapes everything

**Researchers re-run old analysis scripts years later and publish what
comes out.** Every other decision here is downstream of that. The rules
it implies, and the golden master that enforces them, are in
`CONTRIBUTING.md`.

------------------------------------------------------------------------

## Numerics and the random number stream

### `purrr::rbernoulli()` was replaced with `runif()`, not `rbinom()`

An earlier note recommended `stats::rbinom(n, 1, p)`. **That advice was
wrong.** `rbernoulli(n, p)` is internally `runif(n) > (1 - p)`, and
`rbinom` draws from the stream differently — verified across 150
seed/probability combinations. Swapping it in would have silently
changed every reference distribution, and therefore every InfoVal,
computed from a given seed. The `runif` form is bit-identical to the old
behaviour.

**The rule this stands for: check the random *stream*, not just the
distribution.** Two functions with the same marginal distribution are
not interchangeable in a seeded pipeline.

### `rowMeans(x, dims = 2)` was adopted despite not being bit-identical

The patch-averaging step in
[`generateNoiseImage()`](https://rdotsch.github.io/rcicr/reference/generateNoiseImage.md)
moved from `apply(..., 1:2, mean)`, about 6x faster end-to-end. The two
sum in a different order and so differ by roughly 1 ULP (~1e-19 on pixel
values of order 0.01). Adopted because an **independent oracle** — the
average written as an explicit triple loop, using neither function — put
both forms ~5.6e-17 away across noise types, scales and seeds. Neither
is “more correct”, and at the golden master’s own configuration the
results came out bit-identical. Unlike the `rbinom` case above this is
floating-point summation order, not a changed random stream; `NEWS.md`
says so explicitly.

The subtlety that bit once:
**[`rowMeans()`](https://rdrr.io/r/base/colSums.html) on a 3-D array
defaults to `dims = 1`**, collapsing dims 2 *and* 3, and
[`array()`](https://rdrr.io/r/base/array.html) then silently recycles
the short result. The version originally submitted did exactly that —
measured max deviation 0.21 against data with SD ~0.01, not the ~1e-17
it claimed.

### `set.seed()` in `generateStimuli2IFC()` is load-bearing far beyond stimulus generation

[`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
rebuilds the stimuli through it, so that `set.seed(seed)` lands *before*
the simulation loop’s [`runif()`](https://rdrr.io/r/stats/Uniform.html)
draws — which is what makes InfoVal reproducible from a stimulus file
alone, independent of ambient RNG state and of `ncores`. This was
**emergent, not designed**: moving one line would have silently changed
every InfoVal ever computed, without touching
[`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md).
It is now a documented guarantee in
[`?generateReferenceDistribution2IFC`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md),
pinned by a test, with a comment on the
[`set.seed()`](https://rdrr.io/r/base/Random.html) call saying what
depends on it.

### `response_seed`, not `seed`

Four choices, each forced by something specific:

- **The name.** `seed` is an *object in the `.Rdata` file*, and
  [`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
  re-saves its own frame, so an argument of that name would overwrite
  the stimulus seed and write it back — corrupting the record of how the
  stimuli were generated.
- **It seeds the responses, applied after the stimulus rebuild**, not
  forwarded into it. Forwarding would rebuild a different stimulus set,
  so the null would describe stimuli the participants never saw.
- **`NULL` issues no [`set.seed()`](https://rdrr.io/r/base/Random.html)
  call at all**, so the default path is byte-identical rather than
  merely equivalent-looking.
- **`computeInfoVal2IFC(response_seed=)` forces regeneration and never
  caches.** Without the first it would be silently ignored on every file
  after the first, since the generator writes `reference_norms` into the
  file; caching a one-off Monte Carlo check would redefine what every
  later InfoVal from that stimulus set means.

### 4096 → 4092 parameters: why old stimulus files cannot be regenerated from their seed

rcicr 0.3.0 (2015-01-23) cut the per-trial random draws from 4096 to
4092. 4092 is the real patch count — 6 orientations × 2 phases ×
`sum(4^0..4^4)` — while 4096 was a round `2^12` over-allocation, so four
contrasts were drawn per trial that no patch index ever referred to.
`ChangeLog` calls them “redundant” and says the change “does not affect
anything else”.

**That last claim is true for analysis and false for regeneration**,
which is not obvious and matters. Analysing a pre-0.3.0 file reads its
*stored* parameters, and dropping four unused columns changes nothing.
But [`runif()`](https://rdrr.io/r/stats/Uniform.html) is a stream: at
the same seed, trial 1 gets identical values either way and **every
later trial is shifted by four draws** — verified, trial 1 identical,
trials 2 and 3 not. So a pre-0.3.0 stimulus set can be re-analysed
exactly and cannot be re-created from its seed.

The same release also fixed `sinIdx` counting from 0 rather than 1,
which is what the `pre_0.3.0` flag exists for. The two are independent
and both live in that boundary.

**The 0-based path in
[`generateNoiseImage()`](https://rdotsch.github.io/rcicr/reference/generateNoiseImage.md)
looks like it corrupts the whole image and does not.** When
`min(patchIdx) == 0`, `params[p$patchIdx]` drops every 0-indexed cell (R
drops rather than `NA`s a 0 subscript), returning a vector too short for
the patch array, which [`array()`](https://rdrr.io/r/base/array.html)
then recycles — the setup for a general misalignment. It is harmless
because the same counter-from-0 leaves the *last* patch layer never
written, so every `patchIdx == 0` cell is also a `patches == 0` cell —
the one-way implication masking needs, not an equivalence. Those cells
sit contiguously at the end of the column-major order, so recycled
values land only where the patch is identically zero and are multiplied
away. Measured: identical to the honest “one patch not shown” result
(max abs diff 0) across 36 size/nscales/seed combinations, so the
warning is accurate rather than an understatement.

Not a self-fulfilling reconstruction: a stimulus file stores the patch
array its own generator wrote, and the genuine pre-0.3.0 generator is
the identical `co = 0` / `idx = 0` loop the `pre_0.3.0` flag runs —
verified against the R-Forge source (`git show 7d0d9e6:pkg/R/rcicr.R`),
where 0.3.0 only flipped the default and added the flag. Do **not**
“fix” the recycle by 1-offsetting the index: that changes which sinusoid
is dropped, altering the CI of every genuinely pre-0.3.0 file.
`test-generateNoiseImage.R` pins both properties.

The truncation left a dead branch in
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
for eleven years: the single-trial path tested for a length of 4092 and
truncated to 4092, so it could never fire on the 4096-length input it
was written for. **A backward-compatibility path that nothing exercises
is indistinguishable from one that works** — it had no test until
2026-07-28, and neither did the `sinusoids`/`sinIdx` path, also broken.

### `load()` assigns into the calling frame — check every new argument against saved names

An object in an `.Rdata` file silently overwrites a function argument of
the same name.
[`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
re-saved its whole frame, so the files it wrote contained `rdata` and
`ncores`; a second call then ignored the caller’s `ncores` and wrote
back to the path recorded during the first. Fixed at the source, by
excluding the function’s own arguments from the save, *and* defensively
on read for files older versions already wrote. This is also why
`response_seed` is named as it is.

### The InfoVal formula is already correct — do not “fix” it

It matches the erratum to Schmitz et al. (2019): the Euclidean norm,
with *k* supplied by R’s [`mad()`](https://rdrr.io/r/stats/mad.html)
(`constant = 1.4826`). Fixed years ago, now pinned by a regression test.

------------------------------------------------------------------------

## Things that look like bugs and are not

### `autoscale()` leaves `$combined` untouched — intentional

Flagged as a bug and changed; **Ron overruled, and the change was
reverted.** `$combined` must stay as the caller supplied it, so a
combination made before autoscaling survives the call and existing
scripts keep plotting the same image; `$scaled` carries the autoscaled
result. The user-facing problem was never the code but that nothing said
which field to look at.

The trap worth knowing: after
[`batchGenerateCI()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
(which scales `'none'` first) `$combined` is an overlay of *unscaled*
noise and looks almost blank. Build it as `(ci$scaled + ci$base) / 2`,
which is what `save_as_pngs = TRUE` writes.

What kept this cheap to undo: it was filed under its own “Behaviour
change” heading rather than slipped in with unambiguous bug fixes.
**When a fix is debatable, say so.**

### `base_face_files` validation rejects two inputs that used to run

Duplicate names, and a list element that is not a single file name, now
stop the call. Both ran before, so this is a considered exception to the
constraint above — allowed because neither ever produced what the call
asked for. `list(face = 'a.png', face = 'b.png')` looks each name up by
string, so it wrote **one** stimulus set, from `a.png`: verified against
the pre-fix code, four files for two trials. A script that “worked” this
way was generating stimuli for a base image it never named. Every other
new check rejects input that already failed further downstream, inside a
parallel worker.

**Readability is checked by attempting the read, not by
[`file.access()`](https://rdrr.io/r/base/file.access.html)**, which is
documented as unreliable on Windows and networked filesystems and can
reject a file that reads fine. Keeping the reader’s own message also
preserves the difference between “not a PNG” and “permission denied”.

### `_R_CHECK_LIMIT_CORES_` handling caps cores under check only

`default_ncores()` returns 2 when that variable is set and
`detectCores() - 1` otherwise. Only the checker ever sets it, so CRAN’s
two-core policy is met and behaviour at the console is unchanged. The
trade-off — `detectCores()` over-subscribes under cgroup limits — is
accepted deliberately, not unnoticed.

### `generator_version` in old `.Rdata` files is not trustworthy

It was a hardcoded `'0.4.0'` from 2016 until 1.1.0.9000, so every file
written in that range claims to be 0.4.0; `p$generator_version` held the
real value all along. Anything reading the field must treat `'0.4.0'` as
“unknown, somewhere in that range”, accept both a character string and a
`package_version`, and compare with
[`numeric_version()`](https://rdrr.io/r/base/numeric_version.html)
semantics — as text, `'0.10.0' < '0.4.0'` is `TRUE`. The `pre_0.3.0`
compatibility path does **not** key off this field; it detects the old
`sinusoids`/`sinIdx` layout structurally. Just as well.

### `return_as_dataframe = TRUE` returns one noise image per trial, not per trial × base image

[`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)’s
early [`return()`](https://rdrr.io/r/base/function.html) sits *inside*
the per-base-image loop, so with several base images it fires on the
first and the rest never run. Under the default
`use_same_parameters = TRUE` that is correct — every base image shares
one parameter set — and the returned frame has one column per trial, so
it could not represent the alternative anyway. Documented on
`@param return_as_dataframe` rather than changed.

InfoVal is unaffected, checked rather than assumed:
[`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
is the only in-package caller, never passes `use_same_parameters`, and
the first base image’s parameters come from the same leading block of
the RNG stream either way — measured identical, max absolute difference
0. Widening the frame would change the return shape, so it needs a **new
argument**, never a redefinition.

### `computeCumulativeCICorrelation()` does not aggregate repeated stimuli, and its curve ends at 1 by construction

[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
averages the responses to each unique stimulus before building its CI
(`R/generateCI.R:184-191`);
[`computeCumulativeCICorrelation()`](https://rdotsch.github.io/rcicr/reference/computeCumulativeCICorrelation.md)
does not, and walks trials in presentation order. Deliberate —
collapsing repeats would discard exactly the order a cumulative curve is
about.

With no `targetci` the final CI is built from the same un-aggregated
trials as the curve, so the curve **ends at exactly 1** —
self-consistency, not convergence. Three conditions, all measured: it
holds only where the evaluated trials reach the last one (six responses
at `step = 2` stop at the fifth and end at 0.967, at `step = 3` at
0.938); and not at all when responses cancel exactly, since that CI is
uniformly zero and correlating against a constant gives `NA` at *every*
point.

That self-computed final CI equals
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)’s
only under equal repeat counts, where the two are bit-identical. Unequal
counts weight the data differently — each trial equally here, each
unique stimulus equally there: counts 3/1 correlate at 0.845, counts
4/2/1/1 at 0.773.

Documented and pinned by a test rather than changed. Aggregating the
self-computed final CI would move numeric output for anyone calling
without `targetci` *and* stop the curve ending at 1 — a worse default
than the one being fixed. Whether
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)’s
own weighting is right for unbalanced designs is filed separately.

### Repopulating `ref_lookup` costs four measurements — and the two halves stand or fall together

`AGENTS.md` covers what the table is (not a cache, empty since 2018,
every lookup misses). What belongs here is the way out, because either
half done alone is worse than the status quo.

Repopulating means measuring four numbers: `median(reference_norms)` and
`mad(reference_norms)` under the current formula for seed 1, 512px,
10000 iterations at 100/300/500/1000 trials, one
[`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
run each. The alternative is deleting the ~55 lines of matching and
prompt machinery. **Do not do half of either** — delete the machinery
while intending to re-measure and the feature becomes unrecoverable
rather than merely dormant.

------------------------------------------------------------------------

## Testing

### Vacuous assertions have shipped here twice, and both looked fine on the page

Two `batchGenerateCI*` tests titled “computes one CI per group” asserted
only length, names and `dim`: **grouping was never checked**, and a bug
feeding all trials to every group passed green. The other compared
`p`/`stimuli_params` in an `.Rdata` file across a call that runs with
`save_rdata = FALSE` and never writes them back, so it could not have
failed. Mutation testing caught both; reading did not. Hence every
grouping and threshold test now also asserts that the *wrong* answer
differs (a positional split ≠ a by-column split).

A third instance was caught *before* shipping, when the z-map golden
master was written on 2026-07-28, and it shows the shape to watch for.
`zmapmethod = "quick"` ends in
[`scale()`](https://rdrr.io/r/base/scale.html), so its z-map has sum 0
and sd 1 **by construction** — the obvious summary statistics to pin,
and both worthless as value checks. Mutating `sigma` from 3 to 4, a real
change to the output, left sum and sd bit-identical while `sum(abs())`,
`min`, `max` and every individual cell moved. The rule: **before pinning
a summary statistic, ask what the transformation guarantees about it.**
Anything a normalisation fixes carries no information about the values
that went in. They are still asserted, labelled as a check that the
standardisation happened, alongside statistics that actually vary.

### `computeInfoVal2IFC`’s test oracle mirrors the implementation, and is kept anyway

`test-computeInfoVal2IFC.R` recomputes
`(norm(ci, "f") - median(reference_norms)) / mad(reference_norms)` — the
implementation’s own expression — so it pins the *implementation*, not
the published definition. A formula wrong in both places passes.

Left as it is: the formula was hand-checked against the erratum and the
golden master pins the resulting number, so the risk is covered from two
other directions. Recorded so nobody mistakes this test for the
independent check it resembles — the genuinely independent oracle in the
suite is in `test-generateNoiseImage.R`. If it is ever revisited, the
replacement is a worked example from the paper or a hand-computed 2×2 CI
with a known reference vector, not a tidier restatement of the same
expression.

### Pixel assertions have measured the graphics device twice

The practice — uniform backgrounds, colour channels only, comparisons
rather than pinned values — is in `CONTRIBUTING.md`. Here are the two
measurements behind it, made one line apart in the same fix.

**Channel count belongs to the backend.** cairo (Linux, Windows) writes
RGB, macOS quartz writes RGBA, and an opaque alpha plane is a solid
block of 1s — so it adds a distinct value while nothing has been
painted. `expect_length(unique(as.vector(png)), 1)` measured the backend
as much as the drawing. It was the only assertion in the suite counting
distinct values rather than comparing two renders, and so the only one
that failed when the suite first ran on macOS: an image-to-image
comparison stays green because both sides gain the same plane.

**So does the absolute value.** The first attempt at that fix pinned the
flat value to the background grey and failed on macOS, which renders
`bgimage` 0.5 at ~0.573 where cairo gives 0.502 — colour management, the
same mistake being corrected one line above. The check is now an
*ordering*: the same render over a darker background must come out
darker, which survives any monotone transfer function.

### The recovery test uses a permutation null, not a parametric one

`test-recovery.R` gives a simulated observer a known template and
asserts
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
recovers it, scored by `cor(vec(CI), vec(template))`. The null permutes
response labels across trials, which preserves the noise images and the
1/-1 balance and destroys only the response-to-stimulus pairing.

**Why not a t-test on Pearson r:** it would use df = n_pixels − 2 = 1022
and be wildly anticonservative. The basis has 60 patches, so a 32×32 CI
has effective dimensionality ~60, not 1024, and neighbouring pixels are
autocorrelated by construction — measured null correlations reach \|r\|
= 0.454. Every null CI is built from the same basis, so the null
inherits that exact autocorrelation; that is what makes it correct, not
merely convenient.

Sensitivity established by mutation: correct pairing 0.71, sign flip
−0.74, responses reversed −0.18, shifted one trial 0.24. **Documented
limitation:** it cannot catch a transformation applied consistently
inside
[`generateNoiseImage()`](https://rdotsch.github.io/rcicr/reference/generateNoiseImage.md),
because the template is built through the same function and the error
cancels. The oracle test covers that.

### The release gate runs the old code; the golden master only re-runs ours

`test-regression-baseline.R` pins values *this repository computed for
itself*. That makes it self-referential in one specific way: it can only
catch drift away from the moment the numbers were written down. Had a P0
fix already changed results before the baseline was recorded, the
baseline would have pinned the changed values and passed green forever.
`tools/compare-release-output.R` closes that gap by installing the
reference commit into a temporary library and running both versions over
the same battery — it is the only thing here that executes the old code.
The two are complements, not substitutes: the golden master is cheap
enough to run on every commit, the gate costs two package installs and
minutes of compute, so it runs `--quick` on PRs and in full at release.

The battery is snapshotted into the temp directory before either side
runs, so editing the working copy mid-run cannot leave the two sides
comparing different things — which happened here once, and produces a
“difference” that is purely an artefact.

### The legacy `.Rdata` fixtures are committed, not generated when the tests run

Both directions of the compatibility promise now have a test. The gate
above checks that this version still computes what old versions
computed; `test-legacy-rdata.R` checks that it can still *read* what
they wrote, which no other test can, because every other fixture is
built by the current generator.

Generating them needs the old version *installed* — each builds a
cluster whose workers call
[`library(rcicr)`](https://rdotsch.github.io/rcicr/), so sourcing its R
files is not enough, and v1.0.1 additionally needs `raster`, dropped in
\#186. At test time that would put a package install and a network round
trip inside the suite. Instead `tools/make-legacy-rdata.R` installs each
tag into a throwaway library once and the files are committed, so the
safeguard runs in every CI job, on every platform, with no network.
Regenerating is a deliberate act: a red test here means this version can
no longer read a file a researcher already has, which is the failure,
not the fixture.

Each is generated at the era’s **defaults** — `nscales = 5`,
`sigma = 25` for 1.0.1 — so its missing-field fallbacks land on the
right noise basis, which is the situation a returning researcher is
actually in. Any other value would enshrine a wrong null. The gabor row
exists because `sigma` reaches the basis through
[`generateGabor()`](https://rdotsch.github.io/rcicr/reference/generateGabor.md)
alone: sinusoidal norms are identical at `sigma` 25 and 10, so no
sinusoidal fixture can exercise that fallback at all.

Writing them also settled the `generator_version` question with data:
the 1.0.1 and 1.1.0 files really do carry a top-level `'0.4.0'` while
`p$generator_version` holds the truth.

### The v1.0.1 reference is pinned; the previous release is a *second* run, not a replacement

The obvious move once a release is green is to make it the new
reference. It is wrong. Each release would then be compared only against
its predecessor, and a tree could walk away from the published numbers
one tolerated epsilon at a time, every step of the walk “identical to
the last release”. The literature was produced with v1.0.1, so that
comparison is the one that protects it, and it stays pinned at `v1.0.1`
(tagged retroactively at `b6ab269`, so the default reads as a version
rather than a bare SHA).

The second run (`--ref=v1.1.0`) answers a different and also useful
question — did anything break since the last release — and reaches
further, because v1.0.1 *crashes* on calls that v1.1.0 returns numbers
for. Hence `EXPECTED` entries name the reference they apply to: a
deviation from v1.0.1 is not a deviation from v1.1.0, and an entry that
fired for one would be reported as stale by the other.

### The battery stops where the reference version crashes

Measured 2026-07-28 against v1.0.1 on R 4.3.3: it can produce a z-map
**only** at 512px with `zmapdecoration = TRUE`. Undecorated it dies in
`if (bgimage != '')` (“the condition has length \> 1”), and at 64 and
128px in [`plot.new()`](https://rdrr.io/r/graphics/frame.html) (“figure
margins too large” — 64px decorated fails earlier still, on the `plt`
graphical parameter). `mask` is unusable for the same class of reason.
All are fixed here, and none of them can be gated against v1.0.1: **a
fix that turns a crash into a number has no old number to compare
against.** Those paths are covered by the test suite, and by the gate
only from v1.1.0 onward — the `SINCE` table in `tools/compare-harness.R`
records which extras need which reference.

### Tolerances: 8 ULP scaled to the values, plus an 8-bit pixel count

A flat `.Machine$double.eps` is right for a classification image (values
around 0.01) and far too tight for a z-map, whose values run to several
units and whose one-ULP steps are ~4× larger. So the tolerance is
`8 * eps * max(1, max(abs(reference)))`: a few ULP *of the values
involved*. Anything that could only come from a changed random stream,
algorithm or file format — patch indices, drawn parameters, the base
image, the stimulus PNGs — must be bit-identical instead.

**Why 8 and not 1.** These are not single passes over the data. A z-map
is a convolution over 262,144 pixels then a standardisation over the
same values, so a 3.5e-18 difference in the CI arrives as **1.33e-15** —
measured at 512px against v1.0.1. One ULP of the largest value (9.0e-16)
rejects that; eight accepts it and still sits three orders of magnitude
below anything observable.

Widening a tolerance to make a run pass is the move `CONTRIBUTING.md`
warns against, so the distinction matters: what protects this comparison
is not the tolerance but the two **exact** checks beside it — 0 of N
pixels may differ once quantised to 8 bits, and the NA pattern must
match cell for cell. The z-map sigma bug moved 1,282 cells across the
threshold and was caught by the NA check; every numeric tolerance
considered here would have passed it. The 8-bit check is also the one
answering the question a researcher actually has — does the PNG I
publish change — which is why a 1.11e-16 difference in the CI is a pass
rather than an argument.

------------------------------------------------------------------------

## Performance and parallelism

### `ncores == 1` runs in-process instead of building a one-worker cluster

`startBackend()` in `parallel.R` registers `doSEQ` when `ncores < 2`, so
the same `%dopar%` loops run in the current process and **no loop body
changed**. The test suite went from 140s to 4s; under `R CMD check`,
from `[8s/126s]` to `[8s/37s]` — eight seconds of CPU against 126
elapsed was worker startup, not computation.

**Safe because neither parallel loop draws random numbers**, so there is
no per-worker stream to diverge. Verified, not asserted:
`test-parallel-equivalence.R` pins `ncores = 1` and `ncores = 2` to
bit-identical output and spells out what would make this unsafe again.

### Parallelism stays on `parallel` + `foreach`

Not `future`, `snowfall` or MPI. Issue \#66 asked specifically for
snowfall as a single-core fallback; the outcome it wanted is what
`registerDoSEQ()` provides, without a dependency. Issue \#63 asked for
MPI and `doSNOW`; MPI needs the script launched under `mpirun` rather
than a plain R session, a usability cost no speedup here justifies.

`doSNOW` replaced `doParallel` for \#178: only it honours
`.options.snow`, whose `progress` callback runs in the **parent**, so a
bar ticked inside a `%dopar%` body only moves a worker’s private copy.
`registerDoSEQ()` ignores it (measured), so the serial path keeps
in-body ticks behind an `is.null(cl)` guard. A swap, not an addition —
`doParallel` had no other use — though `snow` is new to the graph, not
transitive via `doParallel`.

### Memory: issue \#12 was solved by deletion, not by chunking

The issue proposed spreading stimulus matrices over multiple `.RData`
files. The actual cause was a preallocated
`zeros(img_size, img_size, n_trials)` array — ~1.5 GB at the defaults —
living in the parent environment, so `foreach` exported a full copy to
*every* worker, each of which wrote one slice and discarded the rest.
Each iteration now allocates only its own trial’s noise.

------------------------------------------------------------------------

## Arguments and internal guards

### Write paths are required arguments, not defaults of `tempdir()`

CRAN’s review of 1.2.1 asked us to “omit any default path in writing
functions”. Defaulting to
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) would satisfy the
letter and was **rejected as the more dangerous option**: a script
relying on `./cis` would keep running, silently writing classification
images into a directory deleted when the session ends, surfacing days
later as missing output with nothing connecting it to an upgrade.
Removing the default turns the same script into an immediate error
naming the argument. A breaking change that announces itself beats a
silent relocation of someone’s results.

Cost nothing to verify: the gate reports `max|d| = 0` across 135 checks,
because `tools/compare-harness.R` already passed every path explicitly.

### `captureArgs()` skips required-and-absent arguments, but never defaulted ones

The [`load()`](https://rdrr.io/r/base/load.html) guard snapshots a
function’s arguments and restores them after reading an `.Rdata` file.
Once paths became required, `mget(names(formals()))` started aborting:
it forces the promise, and a wrapper forwarding its own missing argument
—
[`batchGenerateCI()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
passing `targetpath = targetpath` — makes that promise a missing symbol.

The fix must not over-correct.
[`missing()`](https://rdrr.io/r/base/missing.html) reports a *defaulted*
argument missing too when the caller did not supply it, and skipping
those would silently reopen the hazard the guard exists for: a default
is the value the function goes on to use, and is exactly as replaceable
by a field in the `.Rdata` as one passed explicitly. Dropping them
removed the `step` guard in
[`computeCumulativeCICorrelation()`](https://rdotsch.github.io/rcicr/reference/computeCumulativeCICorrelation.md)
and the test caught it. The predicate is therefore *required* (no
default in [`formals()`](https://rdrr.io/r/base/formals.html)) **and**
absent — not absent alone.

## Documentation

### The Medium walkthrough moved into a vignette, and the post stays up

A tutorial outside the repo cannot execute at build time, so it goes
stale silently. The premise proved itself during the port: two of the
post’s lines no longer worked —
[`autoscale()`](https://rdotsch.github.io/rcicr/reference/autoscale.md)’s
`saveasjpegs` argument, now `save_as_pngs`, and
`install_github(..., ref = "development")`, a branch that does not
exist. The post stays published for nine years of inbound links; this
moves the canonical copy, it is not a takedown.

### Three vignette figures were wrong in ways only viewing them showed

The rule this stands behind is in `CONTRIBUTING.md`. What it caught,
none of it visible to an assertion:

- **[`image()`](https://rdrr.io/r/graphics/image.html) auto-stretches
  its input across the palette**, so every linear rescaling of the same
  image renders *identically* — which made a four-way scaling comparison
  meaningless until `zlim = c(0, 1)` was pinned.
- **`zmapmethod = "quick"` z-scores across the pixels of one image**, so
  its values are relative to that image’s own spatial structure, not a
  null distribution. At small sizes they span only ±1.65 and the
  **default `threshold = 3` returns a blank map** — which reads as “no
  signal” but means “wrong ruler”.
- A white-noise synthetic base drowned every figure; a crude Gaussian
  blob is legible and still obviously synthetic.

### When the docs and the code disagreed about `mask`, the code was the contract

Both
[`?plotZmap`](https://rdotsch.github.io/rcicr/reference/plotZmap.md) and
[`?generateCI`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
said a *matrix* masks where the cell is `1`/`TRUE` while a *PNG* masks
where it is black (`0`) — two opposite conventions in one sentence.
`applyMask()` has always done `mask_matrix == 0` unconditionally, so the
matrix half was simply wrong.

**The documentation was corrected to the code, not the reverse**,
because users build masks against observed behaviour: someone whose
masked CI came out right has a mask encoding `0` = masked, whatever the
page said. Changing the code to match the prose would have silently
inverted every existing mask — the one outcome with real cost, since a
mask covers a face and an inverted one looks plausible.

**This nearly shipped backwards.** The first version of the `plotZmap`
fix believed the docs and branched on provenance, masking where `0` for
a PNG and where `TRUE` for a matrix. That would have made the same mask
remove complementary halves in the two functions. It was caught by Ron
asking whether the two forms should really be opposite, and the answer
was in the sibling implementation rather than in any documentation.
**When two documented conventions conflict, run the one that has been
executing for a decade.**

### `plotZmap(mask = ...)` was applied rather than deprecated

It had been documented since 2016 and never worked — the import half
landed with an explicit “todo: applying the mask” that was never picked
up. Filed as a behaviour change, but the risk that framing guards
against could not materialise: nothing has ever been masked, so no
published z-map depends on the current output. Deprecating would have
removed a documented feature on the grounds that it was never built.

**The generalisable part:** “this changes rendered output” is not by
itself an argument for leaving something broken. Ask who is currently
*relying* on the broken output. Here the answer was nobody, and the
entry sat open a day longer than it needed to.

### The `.Rdata` anatomy belongs in `README.md`

It is the only link between the two halves of the package, and nothing
about a stimulus set is recoverable without it. The field-by-field table
there was written by inspecting a real generated file rather than by
reading the [`save()`](https://rdrr.io/r/base/save.html) call — which is
how `trial` was identified as a leftover loop counter carrying no
information.
