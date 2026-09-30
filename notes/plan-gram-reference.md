# Plan: build the InfoVal reference through the stimulus Gram matrix (#354)

## Position

The reference distribution is the norms of CIs built from random responses. rcicr renders every
saved stimulus into a pixels-by-trials matrix and multiplies it by each response vector. The same
norms follow from the trials-by-trials Gram matrix, `norm(S r / n) = sqrt(t(r) G r) / n`, where `G`
is built from the saved parameters and the sparse basis without rendering any image.

The two routes sum in a different order, and none of the configurations measured gave bit-identical results. That relaxes an
exactness contract the repository holds today. `test-generateReferenceDistribution2IFC.R:255` pins
the reference norms bit-for-bit to the rendered arithmetic. The maintainer accepted the trade on
#354, on the measurements in `analyses/gram-reference-accuracy.Rmd`, committed on this branch and
kept after the merge.

## What the analysis measured

From the knitted `analyses/gram-reference-accuracy.md`, all with the reference BLAS and 10,000
reference draws: seven configurations with scored CIs, from 64px with 100 trials to 512px with 300
trials, sinusoid and gabor, plus the package default (512px, 770 trials), compared through its two
references alone.

- **Accuracy.** No configuration gave bit-identical references. The largest relative difference in
  a single norm was 4.8e-14 (4.6e-14 at the default). The largest InfoVal difference among 2,800
  CIs was 7.9e-13, at 512px with 300 trials.
- **Decisions.** The norms whose call at 1.96 differs between the routes span at most 5.6e-13 on
  the rendered route's InfoVal scale (2.3e-14 at the default), so only a CI within that distance of
  the cut-off can change its call. That distance is 5.6e10 times smaller than the estimated
  standard deviation of InfoVal at 1.96 across 10,000-draw references (0.031, from 40). No call
  changed among the 2,800 CIs, 228 of them within 0.5 of 1.96.
- **Speed, serial only.** Against the rendered route rendering serially, as the package does with
  `ncores = 1`: 13x or more faster at 128px, and over 100x at 512px. The package's default
  `ncores = detectCores() - 1` parallelises only the rendering. Its per-iteration loop, which
  dominates at 10,000 draws, is serial whatever `ncores` is. The analysis does not time the
  parallel default, so these ratios do not describe it.
- **Memory at the default.** Peak R heap was 3.3 GB for the rendered route and 1.6 GB for the Gram
  route, most of the latter the sparse basis and the loaded file. As component sizes, the rendered
  noise matrix is 1.5 GB, and the Gram matrix is 4.5 MB plus a 128 MB basis cross-product.

## Changes

1. **One helper builds the norms for all three paths:** the shared default in
   `generateReferenceDistribution2IFC()`, `generateBaseReference()`, and
   `generateSubsetReference()`. Today each has its own copy of the loop. Signature
   `referenceNorms(source, label, reference_stimuli, iter)`, called after `seedResponseStream()`,
   which is unchanged. It:
   - takes the saved parameters from `savedReferenceParams()`, which already truncates pre-0.3.0
     files' 4,096 columns to 4,092, and selects `reference_stimuli` rows if given;
   - builds the sparse basis `P` from `p` (or `s`), mirroring `generateNoiseImage()`:
     - the legacy `sinusoids`/`sinIdx` names are renamed;
     - 0-based `patchIdx` cells are dropped, since their `patches` are 0, the property
       `DECISIONS.md` → "4096 → 4092" measures;
     - parameters beyond `max(patchIdx)` are ignored, with `generateNoiseImage()`'s warning given
       once rather than once per trial;
   - computes `G = X t(P) P t(X)`;
   - draws responses in blocks as one `runif(n * k)`, which consumes the stream exactly as `k`
     sequential `runif(n)` calls do, so the RNG state afterwards is unchanged.

   The progress bar ticks per block.
2. **`ncores`** stays an argument, because scripts pass it. The reference no longer renders noise,
   so it has nothing to parallelise. Its documentation says so. No warning, because passing it is
   not an error.
3. **`Matrix` moves into `Imports`.** It is an R "recommended" package, already in rcicr's recursive
   dependencies through `spatstat.explore`, `spatstat.data`, `spatstat.random` and
   `spatstat.sparse`, so nothing new is installed.
4. **Stored references are untouched.** A cached reference is reused as stored, and its fingerprint
   check compares the stored norms with their stored copy, not with a recomputation. Only
   references computed fresh change: new files, `force_gen_ref_dist`, `response_seed`, and the
   automatic rebuild of an unmarked legacy reference.
5. **Docs.** A `NEWS.md` entry under "Reproducibility impact", sized from the knitted analysis, and
   a performance note. Its figures come from timing `generateReferenceDistribution2IFC()` itself on
   `main` and on the branch, at the default `ncores`, in the implementation PR; the analysis's
   serial ratios are not quoted there. `DECISIONS.md` generalises "`rowMeans(x, dims = 2)` was adopted
   despite not being bit-identical" to cover both, since they share a rationale: an independent
   oracle, and measured differences far below the Monte Carlo error InfoVal already carries. Both
   entries state the bound only for what was measured, the configurations in the analysis and
   the BLAS builds CI ran, and name what was not. The file is at 5,197 of 5,200 words, so the
   merged entry must fit by tightening the existing text, not adding to it.

## Risks

- **Most likely to fail: exact legacy parity in the basis.** The 0-based and `sinusoids` paths in
  `generateNoiseImage()` are exercised only by committed legacy fixtures (`test-legacy-rdata.R`).
  A mismatch there shows up as a large difference, not a rounding one. Mitigation: a test comparing
  the Gram reference with the rendered arithmetic on every legacy fixture, at the same tolerance as
  the rest.
- **The release gate.** InfoVal checks use an absolute tolerance of `1e-9`
  (`tools/compare-release-output.R:350`). The measured differences are three orders of magnitude
  below that, so no `EXPECTED` entry should be needed. That is a prediction; the gate run against
  `origin/main`, and the full dispatch, settle it.
- **Memory at 512px.** The Gram route peaked at 1.6 GB of R heap at the default in the analysis,
  mostly the sparse basis (15.7M non-zeros) and the loaded file. That is half the rendered route's
  3.3 GB, but not small. If the implementation cannot lower it, the PR states its measured peak.

## Verification

- The bit-identity test becomes a tolerance test, `expect_equal(tolerance = 1e-12)`, against the
  rendered arithmetic it already spells out. The RNG state afterwards and the saved fields stay
  `identical()`. The analysis measured only the reference BLAS. The four CI platforms (macOS,
  Windows, and Ubuntu release and devel) run this test, so they measure the tolerance under their
  own BLAS builds. If any exceeds `1e-12`, the PR reports the measured difference and sizes the
  tolerance and the `NEWS.md` entry from it, rather than loosening either unmeasured.
- A second parity test at 128px with 300 trials, a few seconds' work, gives those four BLAS builds
  a configuration of realistic size, not only the 32px fixture. The 512px default under a
  non-reference BLAS stays unmeasured, and the docs say so.
- The legacy fixtures are checked against the rendered arithmetic, as under Risks.
- The full suite passes, and `analyses/gram-reference-accuracy.Rmd` knits against the implementation.
  Its package check requires agreement with both routes to 1e-12.
- The quick gate against `origin/main` locally, and the full `reproducibility.yaml` dispatch in CI.
- Mutation: breaking the 0-based handling, or the block stream order, must fail a test.
