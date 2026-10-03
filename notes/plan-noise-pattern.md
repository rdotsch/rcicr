# Plan: compute each scale's patch size directly (#386)

Part of #384. `generateNoisePattern()` builds three `img_size` × `img_size` × `nscales` arrays (`matlab::meshgrid()`, about 10 MB each at 512px and 5 scales) and divides two of them to read one element per scale: `img_size / scale`.

## The change

- `size <- img_size / scale`; `mg` and `patchSize` go.
- `for (col in 1:scale)` and `for (row in 1:scale)` become `seq_len(scale)` (`scale >= 1`, so identical).

Verified on #386: `identical(patchSize[s, img_size], img_size / s)` for every valid `img_size` in 16, 64, 512 and `nscales` 1 to 5.

## Proving nothing changed

- **Full** release gate against `main` (`Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`): `0 expected deviations, 0 unexpected`.
- Full `testthat::test_local()` passing with **no existing test modified**.
- `lintr::lint_package()`: 0 lints. The gate compares `patchIdx` exactly and the patches within ULPs.

## Step most likely to fail

None of substance; the gate is the check.
