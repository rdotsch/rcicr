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

From the knitted `analyses/gram-reference-accuracy.md`: seven configurations, from 64px with 100
trials to 512px with 300 trials, sinusoid and gabor, each with 10,000 reference draws.

- **Accuracy.** The largest relative difference in a single norm was 4.8e-14. The largest InfoVal
  difference among 2,800 CIs was 7.9e-13, at 512px.
- **Decisions.** At 1.96 the InfoVal difference was at most 5.6e-13, so only a CI within that
  distance of the cut-off can change its call. That distance is 5.6e10 times smaller than the Monte
  Carlo spread of InfoVal at 1.96 for a 10,000-draw reference (0.031). No call changed among the
  2,800 CIs, 228 of them within 0.5 of 1.96.
- **Speed.** 13x to 107x faster at 128px or more.
- **Memory.** At the 512px, 770-trial default, the rendered noise matrix is 1.5 GB; the Gram matrix
  is 4.5 MB, plus a 128 MB basis cross-product.

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
   a performance note. `DECISIONS.md` generalises "`rowMeans(x, dims = 2)` was adopted despite not
   being bit-identical" to cover both, since they share a rationale: an independent oracle, and
   differences many orders of magnitude below anything a researcher reports. The file is at 5,197
   of 5,200 words, so the merged entry must fit by tightening the existing text, not adding to it.

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
- **Memory at 512px.** Building the sparse basis peaked at 1.4 GB in a Gram-only run at the default, 512px, 770 trials and 10,000 iterations (15.7M non-zeros,
  plus the loaded file). That is below the rendered route, but not small. If the implementation
  cannot lower it, the PR states the measured peak.

## Verification

- The bit-identity test becomes a tolerance test, `expect_equal(tolerance = 1e-12)`, against the
  rendered arithmetic it already spells out. The RNG state afterwards and the saved fields stay
  `identical()`.
- The legacy fixtures are checked against the rendered arithmetic, as under Risks.
- The full suite passes, and `analyses/gram-reference-accuracy.Rmd` knits against the implementation.
  Its package check requires agreement with both routes to 1e-12.
- The quick gate against `origin/main` locally, and the full `reproducibility.yaml` dispatch in CI.
- Mutation: breaking the 0-based handling, or the block stream order, must fail a test.
