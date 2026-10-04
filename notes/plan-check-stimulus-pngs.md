# Plan: `checkStimulusPNGs()`

## What it is

An exported function that checks a stimulus `.Rdata` file against a folder of stimulus PNGs, trial by trial:

```r
checkStimulusPNGs(rdata, png_dir, label = NULL, seed = NULL)
```

It returns a data frame with one row per base label and trial: `base`, `trial`, `share` (the share of unclipped pixels that agree), `compared` (how many pixels were compared) and `missing` (`TRUE` when either PNG is absent).

`vignette("recipes")` → "When the `.Rdata` file is lost" then calls it, and its `matches()` helper and the copy of the rendering formula it holds go.

## Why a function, not the vignette's helper

- **The formula belongs to the package.** A stimulus is `renderStimulus(noise, base)`, `clampUnit((((noise + 0.3) / 0.6) + base) / 2)` (`R/generateStimuli2IFC.R:302`). The vignette's helper copies `0.3` and `0.6`. The function reads them from one place that `renderStimulus()` also uses, so the two cannot drift.
- **The edge cases are real.** All of these were found in review of #412, and each one made the helper report a wrong answer:
  - PNG names have to be built as `stimulusPngPath()` builds them. A non-whole seed such as `1.5` broke `%d`, and a label containing `_ori.png` broke substitution.
  - Every base label has to be checked. The first base draws the same parameters whether or not `use_same_parameters` is right.
  - RGB PNGs need one channel.
- **`label` and `seed` describe the archive, not the candidate.** They default to the file's own values, which are right when the file is the original. To test a candidate regenerated with other settings, including another seed, pass the archive's `label` and `seed`: the PNG names carry them. Taking them from a candidate with the wrong seed would look for names the archive does not have. Base labels come from the file, since analysis code needs the file's own.

## Behaviour

- **Comparison.** `ori - inv` equals the noise divided by 0.6 wherever neither image is clipped at 0 or 1, because both hold the same base. A pixel agrees when `abs(ori - inv - noise / 0.6) <= 1 / 255`.
- **Why 1/255.** Measured over 20 trials at 128 px (run and output below), the largest residual on unclipped pixels was 0.5/255. That held for PNGs written by this version and by 1.0.1, at `nscales = 1` and `4`. Pixels that 1.0.1 wrote dark when they overflowed differ by about 256/255; they are 0.13% of pixels at `nscales = 1`. The vignette's 2/255 was looser than needed.
- **Parameters.** Read through `loadStimulusParams()`, which converts legacy `s` files, then each base through `selectStimulusParams()`, which truncates pre-0.3.0 4096-wide matrices to 4092, as `generateCI()` does. `loadStimulusParams()` alone does not truncate, and `generateNoiseImage()` rejects a 4096-element row. Noise is rendered with `generateNoiseImage()`.
- **Missing PNGs.** A row with `missing = TRUE` and `share = NA`, plus one warning naming the first few missing files (`firstFew()`). No PNG found at all stops with an error naming the pattern looked for, and saying that `label` and `seed` must be the archive's.
- **No verdict.** It returns shares, not `TRUE`/`FALSE`. The measured gap between right and wrong settings (at least 98.6% against a median of 2% to 18%) is documented in `@return` and the recipe.
- `png_dir` is required, with no default (AGENTS.md → write paths; it is read-only, but a default folder would be a guess).

## Not changed

- `renderStimulus()` keeps its arithmetic. Only its constants move to a shared internal place. The golden master and the release gate confirm stimuli are unchanged.
- No `.Rdata` field, no numeric output, no existing argument.

## Tests (`tests/testthat/test-check-stimulus-pngs.R`)

- Right settings: every share is 1 on a synthetic set.
- Wrong `nscales`: every share is below 0.2. A candidate with the wrong seed, checked with the archive's `seed`: every share below 0.2. The same candidate without `seed`: the error that names `label` and `seed`.
- A pre-0.3.0-shaped file: a current `nscales = 5` file whose matrices are widened to 4096 columns, as `test-ci-inputs.R:171` builds them (no pre-0.3.0 fixture exists), checked against PNGs rendered from the original: shares of 1, through the 4096 → 4092 truncation.
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

## The run behind the tolerance

This becomes `analyses/stimulus-png-residuals.Rmd`, knitted to `.md` beside it, as the other analyses are. 1.0.1 is installed from its tag into a library of its own, as the release gate does. Per archive: a synthetic face at 128 px, `seed = 42`, 20 trials, written with `save_rdata = FALSE`. The noise is regenerated with this version over a grey base. Residuals are `abs(ori - inv - noise / 0.6)` over pixels where neither image is 0 or 1. The script that produced the numbers:

```r
grey <- tempfile(fileext = ".png"); png::writePNG(matrix(0.5, n, n), grey)
generateStimuli2IFC(list(face = grey), n_trials = 20, img_size = n, seed = 42, nscales = ns,
                    ncores = 1, stimulus_path = path, save_as_png = FALSE,
                    maximize_baseimage_contrast = FALSE)
res <- unlist(lapply(1:20, function(i) {
  ori <- png::readPNG(sprintf("%s/rcic_face_42_%05d_ori.png", archive, i))
  inv <- png::readPNG(sprintf("%s/rcic_face_42_%05d_inv.png", archive, i))
  r <- abs(ori - inv - generateNoiseImage(e$stimuli_params$face[i, ], e$p) / 0.6)
  r[ori > 0 & ori < 1 & inv > 0 & inv < 1]
}))
```

Output:

```
current nscales=1: max residual 0.500/255, share <= 1/255 1.00000, share <= 2/255 1.00000, n=327143
v1.0.1  nscales=1: max residual 256.498/255, share <= 1/255 0.99871, share <= 2/255 0.99871, n=327567
current nscales=4: max residual 0.500/255, share <= 1/255 1.00000, share <= 2/255 1.00000, n=327680
v1.0.1  nscales=4: max residual 0.500/255, share <= 1/255 1.00000, share <= 2/255 1.00000, n=327680
```

The 0.13% of 1.0.1's pixels outside 1/255 are off by about 256/255: pixels above white that 1.0.1 wrote as dark, which this development version writes white (`NEWS.md`, "Pixels above white"). The analysis file will state the configurations measured and no others; `@return` cites it rather than restating the numbers.
