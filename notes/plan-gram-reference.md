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
  `ncores = 1`, in the committed knit: 25x to 34x at 64px with 100 trials, 23x at 128px with 300
  trials, 44x and 41x at 256px with 300 and 770 trials, 107x at 512px with 300 trials, and 147x at
  the 770-trial default. Between knits these ratios have varied by up to about 10%. The package's default
  `ncores = detectCores() - 1` parallelises the rendering but not the per-iteration loop. The
  analysis times each route whole, not the two phases, so how far the parallel default narrows
  these ratios is not measured.
- **Memory at the default.** Peak R heap was 3.3 GB for the rendered route and 1.6 GB for the Gram
  route. The analysis records only these totals, not what makes them up. As component sizes, the rendered
  noise matrix is 1.5 GB, and the Gram matrix is 4.5 MB plus a 128 MB basis cross-product.

## Changes

1. **One helper builds the norms for all three paths:** the shared default in
   `generateReferenceDistribution2IFC()`, `generateBaseReference()`, and
   `generateSubsetReference()`. Today each has its own copy of the loop. Signature
   `referenceNorms(source, label, reference_stimuli, iter, ncores)`, called after `seedResponseStream()`,
   which is unchanged. It:
   - takes the saved parameters from `savedReferenceParams()`, which already truncates pre-0.3.0
     files' 4,096 columns to 4,092, and selects `reference_stimuli` rows if given;
   - builds the sparse basis `P` from `p` (or `s`), mirroring `generateNoiseImage()`:
     - the legacy `sinusoids`/`sinIdx` names are renamed;
     - 0-based `patchIdx` cells are dropped, since their `patches` are 0, the property
       `DECISIONS.md` → "4096 → 4092" measures;
     - parameters beyond `max(patchIdx)` are ignored, with `generateNoiseImage()`'s warning given
       once rather than once per trial;
   - **routes on size.** When the selected trials outnumber the pixels, `G` (trials x trials) would
     be larger than the rendered noise, and so would each iteration's product. Such sets, like
     32px with 10,000 trials, keep today's rendered calculation unchanged and stay bit-identical.
     Every other set uses `G`. All eight configurations in the analysis have fewer trials than
     pixels, so all take the Gram route it measures;
   - builds `G` whichever way holds less: as `crossprod(S)` from the sparse-rendered noise
     `S = P t(X)` when pixels x trials is at most the squared parameter count, otherwise as
     `X t(P) P t(X)`, whose dense cross-product is 128 MB for five scales at any image size. Small
     stimulus sets therefore never pay for it. The analysis implements the same rule and records
     which build each configuration used;
   - draws responses in blocks as one `runif(n * k)`, which consumes the stream exactly as `k`
     sequential `runif(n)` calls do, so the RNG state afterwards is unchanged.

   The progress bar ticks per block.
2. **`ncores`** keeps its meaning where rendering remains. On the rendered fallback (more trials
   than pixels) it goes to `referenceNoise()` exactly as today, so that path stays parallel and
   bit-identical. The Gram path renders no image, so it has nothing to parallelise and ignores
   `ncores`. The documentation says which path uses it. There is no warning, because passing it is
   not an error.
3. **`Matrix` moves into `Imports`.** It is an R "recommended" package, already in rcicr's recursive
   dependencies through `spatstat.explore`, `spatstat.data`, `spatstat.random` and
   `spatstat.sparse`, so nothing new is installed.
4. **Stored references change only when rebuilt.** A marked (vouched) or seeded reference is reused
   as stored: its fingerprint check compares the stored norms with their stored copy, not with a
   recomputation. Everything computed fresh changes: new files, `force_gen_ref_dist`,
   `response_seed`, and the automatic rebuild of an unmarked default-stream reference left by an
   older rcicr, which `resolveReferenceNorms()` treats as stale. The `NEWS.md` entry names that
   migration explicitly, because it changes a number already stored in a researcher's file.
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
  half the rendered route's 3.3 GB, but not small. The analysis does not break that peak down. If
  the implementation cannot lower it, the PR measures what makes it up and states the peak.

## Verification

- A fixture with more trials than pixels (8px, one scale, 100 trials) takes the rendered route,
  and its reference stays `identical()` to the rendered arithmetic. A one-pixel mutation of the
  routing threshold must fail a test.
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
