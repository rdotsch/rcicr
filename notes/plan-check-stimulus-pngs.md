# Plan: `checkStimulusPNGs()`

## What it is

An exported function that checks a stimulus `.Rdata` file against a folder of stimulus PNGs, trial by trial:

```r
checkStimulusPNGs(rdata, png_dir)
```

It returns a data frame with one row per base label and trial: `base`, `trial`, `share` (the share of unclipped pixels that agree), `compared` (how many pixels were compared) and `missing` (`TRUE` when either PNG is absent).

`vignette("recipes")` → "When the `.Rdata` file is lost" then calls it, and its `matches()` helper and the copy of the rendering formula it holds go.

## Why a function, not the vignette's helper

- **The formula belongs to the package.** A stimulus is `renderStimulus(noise, base)`, `clampUnit((((noise + 0.3) / 0.6) + base) / 2)` (`R/generateStimuli2IFC.R:302`). The vignette's helper copies `0.3` and `0.6`. The function reads them from one place that `renderStimulus()` also uses, so the two cannot drift.
- **The edge cases are real.** All of these were found in review of #412, and each one made the helper report a wrong answer:
  - PNG names have to be built as `stimulusPngPath()` builds them. A non-whole seed such as `1.5` broke `%d`, and a label containing `_ori.png` broke substitution.
  - Every base label has to be checked. The first base draws the same parameters whether or not `use_same_parameters` is right.
  - RGB PNGs need one channel.
- **Label, seed and base labels come from the file.** The function needs only `rdata` and `png_dir`, so there is nothing to mistype.

## Behaviour

- **Comparison.** `ori - inv` equals the noise divided by 0.6 wherever neither image is clipped at 0 or 1, because both hold the same base. A pixel agrees when `abs(ori - inv - noise / 0.6) <= 1 / 255`.
- **Why 1/255.** Measured over 20 trials at 128 px, the largest residual on unclipped pixels was 0.5/255. That held for PNGs written by this version and by 1.0.1, at `nscales = 1` and `4`. Pixels that 1.0.1 wrote dark when they overflowed differ by about 256/255; they are 0.13% of pixels at `nscales = 1`. The vignette's 2/255 was looser than needed.
- **Parameters.** Read through `loadStimulusParams()`, so legacy `s` files and pre-0.3.0 4096-wide matrices are handled as `generateCI()` handles them. Noise is rendered with `generateNoiseImage()`.
- **Missing PNGs.** A row with `missing = TRUE` and `share = NA`, plus one warning naming the first few missing files (`firstFew()`). No PNG found at all stops with an error naming the pattern looked for.
- **No verdict.** It returns shares, not `TRUE`/`FALSE`. The measured gap between right and wrong settings (at least 98.6% against a median of 2% to 18%) is documented in `@return` and the recipe.
- `png_dir` is required, with no default (AGENTS.md → write paths; it is read-only, but a default folder would be a guess).

## Not changed

- `renderStimulus()` keeps its arithmetic. Only its constants move to a shared internal place. The golden master and the release gate confirm stimuli are unchanged.
- No `.Rdata` field, no numeric output, no existing argument.

## Tests (`tests/testthat/test-check-stimulus-pngs.R`)

- Right settings: every share is 1 on a synthetic set.
- Wrong `nscales`, and the wrong seed: every share is below 0.2.
- Two bases with `use_same_parameters = FALSE`, checked against a file regenerated with `TRUE`: the first base is all 1 and the second is below 0.2. The function reports the second base, so the mismatch shows.
- A seed of `1.5`, and a label containing `_ori.png`: PNGs are found.
- RGB PNGs: same shares as grey.
- One PNG deleted: that row is `missing`, with a warning. All deleted: an error.

## Documentation

- The recipe calls `checkStimulusPNGs()` in place of `matches()`.
- `NEWS.md` gets a new-features bullet.
- `_pkgdown.yml` lists the function under "Generating stimuli".
- `R/zzz.R` gets globals only if needed.

## Most likely to fail

Moving the rendering constants without changing a single stimulus pixel. The golden master and the gate's `stimulus_pngs` MD5s must stay identical. If the refactor risks that, the fallback is to leave `renderStimulus()` untouched, add the constants beside it, and add a test asserting that the two agree.
