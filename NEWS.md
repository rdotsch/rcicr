# rcicr 1.0.1.9000 (development version)

> The `.9000` suffix marks an unreleased development version, per the usual R
> convention. It will become a release number when these fixes ship — most
> likely **1.1.0** rather than 1.0.2, since some of them change behaviour
> (see below) rather than only repairing it. The public API is unchanged, so a
> 2.0.0 is not warranted.

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

### Unchanged and verified

- **The infoVal formula itself.** rcicr already implements the corrected version
  from the erratum to Schmitz et al. (2019): the Euclidean norm, with the *k*
  constant supplied by R's `mad()` (`constant = 1.4826`). This was fixed years
  ago; it is now covered by a regression test so it cannot silently drift.
- **The noise basis** (`generateNoisePattern()`), the stimulus parameter draw, and
  all four scaling methods, for the default configuration.

---

*Older changes are recorded in `ChangeLog`.*
