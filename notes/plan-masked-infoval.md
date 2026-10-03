# Plan: InfoVal for masked classification images (#374)

## The defect

`computeInfoVal2IFC()` and `batchComputeInfoVal2IFC()` return `NA` for a CI made with `mask`: `generateCI()` stores `NA` in masked pixels, and the Frobenius norm of a matrix holding `NA` is `NA`. Dropping the `NA`s alone would be wrong, because the reference norms are over the whole image.

## The definition

The InfoVal of a masked CI is computed exactly as for an unmasked one, with both sides restricted to the same pixels: the CI's norm over its unmasked pixels, against the norms of random-response CIs over those same pixels, from the same stimuli and the same replayed response stream. That is Brinkman et al.'s (2019) statistic on the region the researcher analysed. A DECISIONS.md entry records the definition and why dropping `NA`s against the whole-image reference was rejected.

## The change

1. **Which pixels.** The mask is read from the CI itself: the pixels of `target_ci$ci` that are `NA`. No new argument to the scoring functions, so the CI cannot be scored against a mask other than its own. A CI with no unmasked pixel stops with an error.
2. **The reference.** The restricted reference that #349 added for `reference_stimuli` is generalised to restrict pixels too:
   - `stimulusGram(params, p, pixels = NULL)`: with `pixels`, only those rows of each basis block enter the Gram matrix (both the sparse cross-product build and the rendered build). `referenceNoise()` subsets rows the same way for `reference_method = "images"`.
   - The response stream is the same replay as every other reference, so a masked reference is reproducible from the stimulus file alone.
3. **Storage.** Stored in `reference_norms_by_stimuli` (unreleased, added in this development cycle, so no released version reads it), with one more key, `mask`: `NULL` for unmasked, otherwise the run-length encoding of the masked pixels (`rle()` of the logical vector), which is exact under `identical()` and small for any mask drawn as a shape. An entry with `reference_stimuli = NULL` and a mask means "all stimuli, these pixels". Unmasked CIs never reach this path.
4. **`generateReferenceDistribution2IFC(mask = NA)`**, appended as the last formal and accepting what `generateCI(mask =)` accepts, so a seedless or read-only file can be given a stored, seeded masked reference, which is what the existing advice messages tell users to do. `mask` is removed from the frame before `load()`, as `reference_stimuli` is, and listed as an internal so the shared path never saves it into the file.
5. **Batch.** References are grouped by stimuli and mask together, so CIs sharing both share one simulation.

## Reproducibility

An unmasked CI takes exactly today's code path; no existing number moves. The gate run with `--ref=$(git rev-parse origin/main)` must show 0 deviations.

## Tests

- **Oracle:** for a small file, the masked Gram norms equal norms computed by brute force (render every stimulus, drop masked pixels, multiply by the same responses) to 1e-12, under both reference methods.
- A masked CI's InfoVal equals `(norm over kept pixels - median) / mad` of that reference; it is finite where it used to be `NA` (fails without the fix).
- Two different masks get two stored entries; the same mask is reused (simulation mocked to fail on the second call); an unmasked CI still uses `reference_norms`.
- `generateReferenceDistribution2IFC(mask =)` then `computeInfoVal2IFC()` reuses the stored reference; the `.Rdata` gains no `mask` object.
- Batch: mixed masked and unmasked CIs match the single-CI values exactly.

## Step most likely to fail

The pixel restriction inside `stimulusGram()`'s block loop: rows of a block are column-major pixels of that block of columns, so the kept-pixel index has to be cut per block. The brute-force oracle is what catches an off-by-one there.
