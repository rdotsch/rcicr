# Plan: standardise a masked quick z-map over its kept pixels (#429)

## The change

`computeZmapQuick()` (`R/zmap-compute.R`) sets the masked pixels of the blurred CI to `NA` before `scale()`:

```r
zmap <- as.matrix(blur(as.im(ci), sigma = sigma))
zmap[is.na(ci)] <- NA
zmap <- matrix(scale(as.vector(zmap)), img_size, img_size)
```

`generateCI()` applies the mask before the z-map, so `is.na(ci)` is exactly the mask. `scale()` already drops `NA` from both the mean and the SD (checked: `scale(c(1, 2, NA, 4))` equals `(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)`). The t-test z-map is untouched: it already returns `NA` where the CI is `NA`.

## What it changes, measured

The issue's reproduction (128px, pure noise, central quarter kept, `sigma = 3`), current code against the change:

| | current | with the change |
|---|---|---|
| unmasked z-map | | `identical()` to current |
| z ratio, masked / unmasked, median over the core | 1.98 | 1.04 |
| core pixels with abs(z) > 3 (unmasked: 0) | 183 of 1024 | 0 of 1024 |
| masked pixels that are not `NA` (before thresholding) | 12288 of 12288 | 0 |
| mean and SD over kept pixels | | 0 and 1 |

The remaining 1.04 is the mask-edge effect the issue names as separate: kept pixels next to the mask are blurred against the zeros `blur()` returns there, which lowers the kept region's SD a little. Out of scope here, as is the same effect at the image border, which unmasked z-maps already have.

Only masked quick z-maps change, and only `$zmap` and the z-map PNG. `$ci`, `$scaled`, `$combined` and every unmasked call stay bit-identical.

## Tests

In `tests/testthat/test-fixed-bugs.R`, through `generateCI(mask =, zmap = TRUE, threshold = 0)` on the existing fixtures:

1. Every masked pixel of `$zmap` is `NA`, and no kept pixel is.
2. Over the kept pixels, the mean is 0 and the SD is 1 (tolerance 1e-12).
3. The kept values equal the unmasked call's blurred CI standardised over the same pixels, and **differ** from the current all-pixel standardisation, so the test cannot pass on the old code. Checked by running it with the one-line change stashed.

The unmasked quick z-map is already pinned by `test-regression-baseline.R`; it must stay green unchanged.

## Gate and records

- No battery configuration combines `mask` with a z-map, so the gate against v1.0.1 and v1.6.0 sees no difference. `--ref="$(git rev-parse origin/main)"` must report 0 deviations, showing the change moves nothing the battery measures. Adding a mask + z-map configuration was considered and left out: v1.0.1 cannot run it (`DECISIONS.md`), and the tests above pin the behaviour.
- `NEWS.md`, under the development version: a "Reproducibility impact" entry. Masked quick z-maps from 1.1.0 to 1.6.1 have kept values too large by about 1 / sqrt(share of pixels kept), so they show more pixels above the threshold; recompute them. Masked pixels are now `NA`.
- `DECISIONS.md`: one short entry on why the standardisation uses the kept pixels only, and why the edge effect was left alone.

## Most likely to fail

Test 3's comparison with the unmasked call: `blur()` near the mask edge depends on how `as.im()` treats `NA` pixels, so the "expected" value has to be built from the same blurred masked CI, not from the unmasked one. If that makes the test tautological, it falls back to asserting the measured ratio range on a noise CI.
