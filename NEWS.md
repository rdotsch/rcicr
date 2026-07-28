# rcicr (development version)

## Reproducibility impact — read this if you made z-maps with 1.1.0

- **`generateCI(zmap = TRUE)` smoothed z-maps with the wrong sigma, and ignored the
  `sigma` you passed.** `generateCI()` reads the stimulus set with `load()`, which assigns
  straight into the function's own frame, and 1.1.0 started storing the *noise* `sigma` in
  the `.Rdata` file — the same name as the z-map blur argument. The saved value therefore
  replaced the argument on every call: z-maps were blurred with 25 (the default noise sigma,
  or whatever you generated your stimuli with) instead of the documented default of 3, and
  passing `sigma` explicitly did nothing at all.

  **Who is affected.** Only z-maps, and only from stimulus sets generated with 1.1.0 —
  `.Rdata` files written by 1.0.1 and earlier have no `sigma` field, so their z-maps were
  always correct. The classification images themselves, their scaling, InfoVal and every
  saved number are untouched; `sigma` is used for nothing else. If you produced a z-map from
  a 1.1.0 stimulus set, regenerate it: the blur it received was roughly eight times wider
  than intended, which spreads and weakens exactly the localised signal a z-map exists to
  show.

  Every argument is now kept across the `load()`, so a field added to the `.Rdata` later
  cannot quietly capture another one. Found by `tools/compare-release-output.R`, the release
  gate that compares this tree's output against the last published version — this is the
  first bug it caught.

## Behaviour change

- **`plotZmap(mask = ...)` now actually masks the z-map.** The argument has been documented
  since 2016 — "if a cell evaluates to TRUE, the corresponding zmap pixel will be masked" —
  and until now it did nothing at all: the mask was read from its PNG or matrix, checked
  against the z-map's dimensions, validated as binary, converted to boolean, and then
  discarded before plotting. A correct mask produced an unmasked z-map, with no error and no
  warning. Masked cells are now dropped from the z-map exactly as sub-threshold cells are.

  **Who is affected.** Only direct calls to `plotZmap(mask = ...)`. `generateCI()` does not
  pass `mask` to `plotZmap()` — it masks the classification image itself, via a separate and
  working code path — so z-maps produced through the normal pipeline are unchanged, and no
  stored numbers change anywhere. If you have been passing a mask and your z-maps looked
  unmasked, that is why; they will now come out masked, and the earlier images were wrong
  about which regions carry signal.

  A second bug is fixed alongside: the conversion to boolean set every cell to `FALSE`
  whatever you passed, so even once applied the mask would have masked nothing.

- **The `mask` convention is documented correctly for the first time**, in both
  `?plotZmap` and `?generateCI`. Both said a *matrix* masks where the value is `1`/`TRUE`
  while a *PNG* masks where it is black (`0`) — two opposite conventions in one sentence.
  `generateCI()`'s implementation has always masked where the value is `0`, for a matrix and
  a PNG alike, so the matrix half of the documentation was simply wrong. **The code is
  unchanged and the documentation now matches it**, since existing masks were built against
  the behaviour, not the prose. `plotZmap()` follows the same single convention, so one mask
  can be passed to both functions — which is now asserted by a test rather than assumed.

  If you built a mask by reading `?generateCI` rather than by looking at your output, check
  it: `0`/black/`FALSE` is the region that gets masked away.

## New features

- `generateReferenceDistribution2IFC()` and `computeInfoVal2IFC()` gained a
  `response_seed` argument, so the null distribution InfoVal is scored against can be
  varied deliberately. Until now there was no way to draw a second, independent null from
  the same stimuli — which meant you could not check how much Monte Carlo error your choice
  of `iter` was leaving in your InfoVal. `response_seed` seeds the simulated responses only;
  the stimuli, and so the noise basis the null is built on, are untouched.

  **Existing calls are unaffected.** The default (`NULL`) issues no `set.seed()` call at
  all, so the reference distribution is byte-identical to what earlier versions produced.
  Verified against norms generated before the change, not merely assumed.

  In `computeInfoVal2IFC()`, passing `response_seed` forces the reference distribution to be
  regenerated even when the `.Rdata` file already holds one, and the result is deliberately
  *not* written back — a one-off check of the Monte Carlo error cannot silently become the
  number every later analysis of that stimulus set reports.

- `generateReferenceDistribution2IFC()` gained `save_rdata` (default `TRUE`, i.e. unchanged)
  and now returns the reference distribution invisibly instead of returning nothing, so the
  norms are reachable when you ask it not to write them to the `.Rdata` file.

- The `.Rdata` file gained a `reference_norms_seed` field recording the `response_seed` the
  stored `reference_norms` were generated with (`NULL` for the default). Purely additive;
  files written by earlier versions simply lack it. A stimulus set carrying a deliberately
  varied null is no longer indistinguishable from one carrying the default.

## Bug fixes

- `autoscale()` works on masked classification images. `generateCI(mask = ...)` sets masked
  pixels to `NA` by design, and `autoscale()` took a bare `range()` over them, so the scaling
  constant became `NA` and the call died with `missing value where TRUE/FALSE needed`. The
  scaling constant is now computed from the unmasked pixels, exactly as `generateCI()`'s own
  scaling has always done, and masked pixels stay masked in the result. A CI that is
  *entirely* `NA` now raises an error naming the CI instead of failing the same opaque way.

- `generateCI()` accepts a pre-0.3.0 `.Rdata` file when computing a CI from a **single**
  trial. rcicr 0.3.0 stopped drawing four random contrasts per trial that no patch index
  ever referred to (4096 → 4092), and `generateCI()` has truncated older files ever since —
  but the single-trial branch tested for a length of 4092 and then truncated to 4092, a
  no-op that could never fire on the 4096-parameter input it existed for. Such a call failed
  with `Stimulus generation aborted: number of parameters doesn't equal number of patches!`.
  The multi-trial path was always correct and is unchanged.

- `computeInfoVal2IFC()` and `generateReferenceDistribution2IFC()` work on `.Rdata` files
  written before `noise_type` was saved (#94). Such a file failed outright with
  `object 'noise_type' not found`, and the workaround on record was to load the file and
  assign the variable by hand. It now falls back to `sinusoid` with a loud warning, matching
  how `nscales` is handled — a warning rather than a silent default, because guessing wrong
  means the null is built on a different *kind* of noise than participants saw, and the
  resulting InfoVal would be wrong. Files written by 1.1.0 or later already store the field
  and are unaffected.

- `generateStimuli2IFC(return_as_dataframe = TRUE)` shows its progress bar (#82). The
  `return` handing back each trial's noise exits the entire loop body, so it jumped past the
  progress-bar update and the bar sat at zero for the whole run — on the slowest path there
  is, since `generateReferenceDistribution2IFC()` takes it for every InfoVal.

- The `.Rdata` file written by `generateStimuli2IFC()` now records the rcicr version that
  actually wrote it (#29). `generator_version` was a hardcoded `'0.4.0'` string from 2016
  onwards, so every file produced by 0.4.0 through 1.1.0 claims to come from 0.4.0 —
  useless for the provenance the field exists for, and it disagreed with
  `p$generator_version`, which held the real version all along.

  No result changes: nothing in the package has ever read this field. If your own code
  does, note two things. Existing files cannot be trusted to say what wrote them, so treat
  a stored `'0.4.0'` as "unknown, somewhere between 0.4.0 and 1.1.0" rather than as a
  version. And the field is now a `package_version` object rather than a character string,
  so compare with `utils::packageVersion()` or `numeric_version()`, never as text —
  `'0.10.0' < '0.4.0'` is `TRUE` when compared as strings.

## Documentation

- `README.md` now describes the package's architecture and, more usefully, the **anatomy of
  the `.Rdata` file** field by field. That file is the only link between stimulus generation
  and analysis — nothing about a stimulus set is recoverable without it — and until now its
  contents were documented nowhere a user would look.

- A `CONTRIBUTING.md` sets out how to work on the package, leading with the constraint that
  makes it unusual: researchers re-run old analysis scripts years later and publish the
  results, so numeric output must not change silently.

- `?generateStimuli2IFC` documents a restriction on `return_as_dataframe = TRUE`: the frame
  holds one noise image per trial, so it is meaningful only under the default
  `use_same_parameters = TRUE`. With `use_same_parameters = FALSE` and more than one base
  image, only the first base image's noise comes back — the frame's shape cannot represent
  trial × base image. Behaviour is unchanged, and the files written to disk were never
  affected; the restriction simply was not stated.

- `?generateReferenceDistribution2IFC` now documents as a **guarantee** what was previously
  only true by accident: with the default `response_seed`, the reference distribution — and
  therefore InfoVal — is reproducible from the stimulus `.Rdata` file alone, independent of
  the calling session's random number state and of `ncores`. This held before, but nobody
  had chosen it: it is a consequence of `generateStimuli2IFC()`'s internal `set.seed()`
  landing before the simulation draws. That call now carries a comment saying what depends
  on it, so it is not moved casually.

# rcicr 1.1.0 (2026-07-27)

> First release since 1.0.1, and the version submitted to CRAN to reinstate the
> package after its 2021 archival. The minor bump rather than 1.0.2 is deliberate:
> some of the changes below alter behaviour rather than only repairing it. The
> public API is unchanged — no function, argument or argument meaning was removed
> or redefined — so a 2.0.0 is not warranted. Anyone re-running an old analysis
> script should read the "Reproducibility impact" section below.

## Bug fixes

- `generateCI(zmap = TRUE, zmapdecoration = FALSE)` works at all. The undecorated
  branch of `plotZmap()` tested its background image with `if (bgimage != '')`, a
  condition of length `img_size^2`, which R >= 4.2 treats as an error rather than
  silently taking the first element. `generateCI()` always supplies a background
  image, so this path could never run. Same root cause as the `mask` bug above; the
  decorated branch already used `identical()`.
- `plotZmap(decoration = FALSE)` works on small images. `plot.new()` was called
  before the margins were reset to zero, and it rejects a device too small to hold
  the default margins, so any z-map below roughly 100 pixels failed with
  `figure margins too large`. Rendered output at usual sizes is unchanged.
- `generateCI(zmaptargetpath = ...)` is honoured. The argument was documented and
  accepted but never forwarded to `plotZmap()`, so z-maps were always written to
  `./zmaps` relative to the working directory regardless of what was requested.
- `generateStimuli2IFC()` now validates that the base image matches `img_size`, with an
  error naming the file and both sizes, instead of failing inside a parallel worker with
  `non-conformable arrays` (#124).
- The `mask` argument of `generateCI()` works again. It was unusable on R >= 4.2, was
  restricted to 512px stimuli by a hardcoded size check, and rejected greyscale-as-RGB
  PNG masks even when it had successfully converted them.
- `generateCI()`, `batchGenerateCI()` and `computeCumulativeCICorrelation()` accept
  tibble columns, not only `data.frame` columns (#70, #123). `generateCI()` also checks
  that `stimuli` and `responses` are the same length up front.
- `generateStimuli2IFC()` saves `nscales` and `sigma` to the `.Rdata` file, and
  `generateReferenceDistribution2IFC()` uses them (#81). See below.
- `computeInfoVal2IFC(force_gen_ref_dist = TRUE)` regenerates the reference distribution
  instead of silently reusing the cached one (#113).
- `generateNoiseImage()` accepts pre-0.3.3 `.Rdata` noise patterns using the
  `sinusoids`/`sinIdx` layout. This backward-compatibility path had never worked.
- `simulateNoiseIntensities()` runs at all, and honours its `img_size` argument.
- `generateReferenceDistribution2IFC()` gained an `ncores` argument instead of always
  using `detectCores() - 1`.
- `computeInfoVal2IFC()` no longer creates a stray file called `data` in your working
  directory. One `write()` call omitted its destination, and base R's `write()` defaults
  to `file = "data"` — so the "reference distribution has been saved" message was written
  to that file instead of the console.
- `generateCI()` accepts a greyscale-encoded RGB PNG as a `mask`. The conversion worked,
  but the error for genuinely non-greyscale images was raised unconditionally afterwards,
  so every such mask was rejected anyway.
- `generateReferenceDistribution2IFC()` no longer writes its own arguments into the
  `.Rdata` file it updates. `load()` assigns straight into the calling function's frame,
  so a stored `rdata` (or, newly, `ncores`) object overwrote the argument of the same name
  on the next call — meaning a second call ignored the `ncores` you passed and wrote back
  to the path recorded during the first call. `computeInfoVal2IFC()` additionally guards
  against this when reading `.Rdata` files written by older versions, which still contain
  the stale path.

## Performance and dependencies

- `generateStimuli2IFC()` no longer preallocates the full stimulus array before starting
  its parallel workers. At the default 512px / 770 trials that array was ~1.5 GB, and
  `foreach` copied it to every worker, which each wrote one slice into and discarded.
  Addresses issue #12.
- Parallel clusters are now stopped via `on.exit()`, so workers are released even when an
  error interrupts the loop. Fixes the "closing unused connections" warnings (issue #50).
- `generateNoiseImage()` is about **6x faster** — 1.66s to 0.28s per call at the default
  512px with `nscales = 5` — and allocates about 30% less memory. The per-pixel average
  across patch layers is now computed with `rowMeans(..., dims = 2)` instead of
  `apply(..., 1:2, mean)`; that step alone is ~31x faster, and what remains is building
  the weighted patch array, which is unavoidable. Because this function is called for
  every trial during stimulus generation and again for every CI and z-map, the saving
  compounds. Thanks to [@hvalev](https://github.com/hvalev), who diagnosed this and
  benchmarked it in #122. See "Reproducibility impact" below — the result is not
  *bit*-identical to the old one.
- **`ncores = 1` no longer builds a cluster.** The parallel loops in
  `generateStimuli2IFC()` and `generateCI()` previously started a one-worker cluster even
  when asked for a single core — a second R process, loading the package, receiving each
  iteration and sending it back. They now run in the current process. The test suite went
  from 140s to 4s, and `R CMD check` from `[8s/126s]` on tests (eight seconds of CPU
  against 126 elapsed, over 22 cluster spawns) to a fraction of that. Results are
  unaffected and this is verified rather than asserted: `ncores = 1` and `ncores = 2`
  produce bit-identical stimulus parameters, noise bases and classification images, now
  pinned by `tests/testthat/test-parallel-equivalence.R`. Neither loop draws random
  numbers — parameters are drawn under `set.seed()` before the loop — so there is no
  worker RNG stream to diverge. Addresses the remainder of issue #50's neighbourhood.
- `Imports` shrank from 27 packages to 15; none of the removed ones were used.
- Deprecated calls replaced: `dplyr::progress_estimated()`, `purrr::rbernoulli()`, and
  `citEntry()`/`personList()` in `inst/CITATION`. The `rbernoulli()` replacement was
  chosen to preserve the random number stream exactly — see below.

## Reproducibility impact — read this if you have published or in-progress results

This release fixes several long-standing bugs. Most of them turn a **crash into
working behaviour**, and therefore cannot change any result you already have —
if the old code errored, there was no number to change.

Only two fixes can change a value you may already have reported. Both are listed
first, with the exact conditions under which they apply.

A golden-master test (`tests/testthat/test-regression-baseline.R`) pins the
numeric output of the default pipeline — noise basis, classification image,
scaling, and infoVal — to the values produced *before* these fixes. It passes,
which is the evidence for the claim that default-configuration results are
unchanged.

### Results CAN differ

**1. infoVal, if you used a non-default `nscales` (issue #81)**

`generateStimuli2IFC()` did not save `nscales` or `sigma` into the `.Rdata` file.
`generateReferenceDistribution2IFC()` re-generates the stimulus set from that
file to build the infoVal null distribution, so with `nscales` missing it
silently fell back to the default `nscales = 5` — building the reference
distribution on a different noise basis than the stimuli participants actually
saw.

- **Affected:** anyone who called `generateStimuli2IFC(..., nscales = <not 5>)`
  and then reported an infoVal (or `sigma` with `noise_type = "gabor"`).
- **Not affected:** the default `nscales = 5`, where the fallback happened to
  coincide with the real value. This is the common case.
- **After the fix:** infoVal for affected stimulus sets will change, because it is
  now computed against the correct null distribution. The old value was wrong.
- **What to do:** if this applies to you, recompute infoVal. Old `.Rdata` files do
  not contain `nscales`, so the fix warns rather than guessing — regenerate the
  stimulus set, or pass the value explicitly.

**2. infoVal, if you passed `force_gen_ref_dist = TRUE` (issue #113)**

The flag was silently ignored whenever a `reference_norms` vector already existed
in the `.Rdata` file, so the reference distribution was reused rather than
recomputed.

- **Affected:** anyone who passed `force_gen_ref_dist = TRUE` expecting a fresh
  simulation. You received the cached distribution instead.
- **Not affected:** the default `force_gen_ref_dist = FALSE`, which was always
  meant to reuse the cache.
- **After the fix:** the flag regenerates the distribution, so infoVal will differ
  from a cached run. Because the reference distribution is simulated, it also
  varies slightly between regenerations — that is expected.

### Results cannot differ

These fixes only affect code paths that previously errored out, so no previously
obtained result changes:

- **Tibble inputs to `generateCI()`** (issues #70, #123) — previously failed with
  `arguments must have same length`.
- **The `mask` argument** — previously failed on R >= 4.2 with
  `the condition has length > 1`, and was restricted to 512px stimuli by a
  hardcoded size check.
- **Base images not matching `img_size`** (issue #124) — previously failed inside a
  parallel worker with `non-conformable arrays`.
- **Pre-0.3.3 `.Rdata` files** using the `sinusoids`/`sinIdx` layout — the
  backward-compatibility path was unreachable and always errored.
- **`simulateNoiseIntensities()`** — errored on every call.
- **Arguments overwritten by `load()`** — only reachable if you renamed or moved an
  `.Rdata` file between calls, in which case the old code either errored with
  `cannot open the connection` or wrote the reference distribution back to the file's
  former path. The `ncores` half affects speed only. No infoVal changes.

### Changed, but below any scale that can matter: the patch average

`generateNoiseImage()` now averages patch layers with `rowMeans(..., dims = 2)`
rather than `apply(..., 1:2, mean)`. These compute the same quantity but sum in a
different order, so they are **not bit-identical**: they differ by about one unit
in the last place, ~1e-19 in absolute terms on pixel values of order 0.01.

This is a different class of thing from the `rbernoulli` case below. That one would
have changed the random *stream* — a large, systematic divergence. This is
floating-point summation order, and it was checked against an independent oracle
(the average written as an explicit triple loop, using neither `apply()` nor
`rowMeans()`) across noise types, spatial scales and seeds: both the old and new
forms sit ~5.6e-17 from that oracle. Neither is "more correct" than the other.

- **Affected:** nothing you would report. No CI, z-map, infoVal or scaling decision
  changes at any precision a paper prints, and the difference is far below the
  noise you would get from re-running with a different seed.
- **At the configuration the golden master pins**, the results came out bit-identical.
- It is recorded here anyway, because the standard for this package is that any
  numeric change is written down rather than discovered later by someone re-running
  a five-year-old script.
- `tests/testthat/test-generateNoiseImage.R` now pins the new implementation against
  the original `apply()` form across 12 configurations, so it cannot drift further.

### Deliberately unchanged: the random number stream

`purrr::rbernoulli(n, p)` (deprecated) is internally `runif(n) > (1 - p)`. The obvious
modern replacement, `rbinom(n, 1, p)`, draws from the stream differently — verified across
150 seed/probability combinations that the two diverge. Swapping it in would have silently
changed every reference distribution, and therefore every infoVal, computed from a given
seed. The `runif` form is used instead, verified bit-identical to the old behaviour.

### Unchanged and verified

- **The infoVal formula itself.** rcicr already implements the corrected version
  from the erratum to Schmitz et al. (2019): the Euclidean norm, with the *k*
  constant supplied by R's `mad()` (`constant = 1.4826`). This was fixed years
  ago; it is now covered by a regression test so it cannot silently drift.
- **The noise basis** (`generateNoisePattern()`), the stimulus parameter draw, and
  all four scaling methods, for the default configuration.

---

*Older changes are recorded in `ChangeLog`.*
