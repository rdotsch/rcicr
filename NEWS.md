# rcicr 1.0.1.9000 (development version)

> The `.9000` suffix marks an unreleased development version, per the usual R
> convention. It will become a release number when these fixes ship — most
> likely **1.1.0** rather than 1.0.2, since some of them change behaviour
> (see below) rather than only repairing it. The public API is unchanged, so a
> 2.0.0 is not warranted.

## Bug fixes

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
