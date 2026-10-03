# Plan: keep out-of-range pixels from wrapping in PNGs or crashing the z-map (#371, #373)

## The defect

Both issues are one cause: a value outside [0, 1] reaches an image writer.

- `png::writePNG()` (png 0.1.9) wraps values above 1 (1.004 is written 0, 1.1 is written 0.098) and clips values below 0. Stimulus PNGs exceed 1 where noise above 0.3 meets a bright base (67% of trials at `nscales = 1`, 0.5% at 4, none measured at the default 5). CI PNGs exceed it under `scaling = "constant"` with a small constant.
- `rasterImage()` in `plotZmap()` stops on any value outside [0, 1], so `generateCI(zmap = TRUE)` aborts under `scaling = "none"` (whenever a base pixel is 0 and the CI is negative there, and contrast maximization always makes one pixel 0) and under an out-of-range `"constant"`.

## The change

1. One internal helper, `clampUnit(x)`: `pmin(pmax(x, 0), 1)`, which keeps `NA` (masked pixels) as `NA`.
2. Apply it to what is **drawn or written**, never to what is **returned**:
   - `generateStimuli2IFC()`: both `writePNG()` calls.
   - `saveToImage()` (group and individual CI PNGs).
   - `autoscale(save_as_pngs = TRUE)`: in range by construction; clamped so all writers agree.
   - `plotZmap()`: both `rasterImage(bgimage, ...)` calls.
   `$combined`, `$scaled` and every number in the `.Rdata` file stay as they are.
3. The `"constant"` range warning says "clipping will occur"; after this it is true, so the wording stays.

## Reproducibility impact

Stimulus PNG bytes change only where a pixel was above 1, which was written as a dark pixel and is now white. CI PNGs change only under out-of-range `"constant"` scaling. No CI, InfoVal or `.Rdata` value changes. NEWS.md gets a "Reproducibility impact" entry giving the affected configurations with the measured rates, and noting that stimuli already shown cannot be corrected.

## Tests

- `clampUnit()`: values below 0, above 1, inside, `NA`.
- A stimulus set at `nscales = 1` on a white base: no written pixel is darker than the computed stimulus by more than quantisation where the computed value is at least 1 (fails without the fix: those pixels read back near 0).
- `generateCI(scaling = "constant", scaling_constant = 0.01)` on a bright base: every PNG pixel whose `combined` exceeds 1 reads back as 1.
- `generateCI(zmap = TRUE, scaling = "none")` for both decorations writes a z-map instead of erroring.

## Verification, and the step most likely to fail

The gate's stimulus PNG MD5s are all at `nscales = 5`, where no overflow was measured, so they should not move, but that is a sample. **I will run the gate with `--ref=$(git rev-parse origin/main)`**: any MD5 that moves is a pixel that was wrapping in a pinned config and needs an `EXPECTED` entry, not a quiet acceptance.
