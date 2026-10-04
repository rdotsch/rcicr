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
- **Why 1/255.** Measured over 20 trials at 128 px (run and output below), on PNGs written by this version and by 1.0.1 at `nscales = 1` and `4`: every unclipped pixel's residual was at most 0.5/255, the 8-bit rounding, except the pixels 1.0.1 wrote dark on overflow. Those 424 pixels (0.13%, all at `nscales = 1` from 1.0.1) are off by 255.5/255 to 256.5/255. No residual fell between the two. The vignette's 2/255 was looser than needed.
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

This becomes `analyses/stimulus-png-residuals.Rmd`, knitted to `.md` beside it, as the other analyses are. These are the scripts that produced the numbers, as run.

1.0.1 is installed from its tag into a library of its own:

```sh
git archive v1.0.1 | tar -x -C v101
R CMD INSTALL --no-test-load -l lib101 v101
```

`write_archive.R` writes 20 stimuli of a synthetic face, without the `.Rdata` file. It is run once under each version:

```r
library(rcicr); out <- commandArgs(TRUE)[1]; ns <- as.integer(commandArgs(TRUE)[2]); n <- 128
x <- matrix(seq(-1, 1, length.out = n), n, n, byrow = TRUE); y <- -matrix(seq(-1, 1, length.out = n), n, n)
face <- exp(-(x^2 / 0.45 + y^2 / 0.75)); face <- (face - min(face)) / diff(range(face))
dir.create(out); b <- file.path(out, "base.png"); png::writePNG(face, b)
invisible(capture.output(generateStimuli2IFC(list(face = b), n_trials = 20, img_size = n, stimulus_path = out, seed = 42, nscales = ns, ncores = 1, save_rdata = FALSE)))
```

`residuals.R` regenerates the noise with this version over a grey base and compares:

```r
library(rcicr); args <- commandArgs(TRUE); archive <- args[1]; ns <- as.integer(args[2]); n <- 128
q <- function(e) { capture.output(v <- suppressWarnings(suppressMessages(e))); v }
grey <- tempfile(fileext = ".png"); png::writePNG(matrix(0.5, n, n), grey)
path <- tempfile(); q(generateStimuli2IFC(list(face = grey), n_trials = 20, img_size = n, seed = 42, nscales = ns, ncores = 1, stimulus_path = path, save_as_png = FALSE, maximize_baseimage_contrast = FALSE))
e <- new.env(); load(list.files(path, "Rdata$", full.names = TRUE), envir = e)
res <- unlist(lapply(1:20, function(i) {
  f <- file.path(archive, sprintf("rcic_face_42_%05d_ori.png", i)); ori <- png::readPNG(f); inv <- png::readPNG(sub("_ori", "_inv", f))
  if (length(dim(ori)) == 3) { ori <- ori[, , 1]; inv <- inv[, , 1] }
  r <- abs(ori - inv - generateNoiseImage(e$stimuli_params$face[i, ], e$p) / 0.6); u <- ori > 0 & ori < 1 & inv > 0 & inv < 1; r[u] }))
cat(sprintf("max residual %.3f/255, share <= 1/255 %.5f, share <= 2/255 %.5f, n=%d\n", max(res) * 255, mean(res <= 1/255 + 1e-12), mean(res <= 2/255 + 1e-12), length(res)))
```

Driver:

```sh
for ns in 1 4; do
  Rscript write_archive.R cur$ns $ns                     # this version
  R_LIBS=lib101 Rscript write_archive.R arch$ns $ns      # 1.0.1
  Rscript residuals.R cur$ns $ns; Rscript residuals.R arch$ns $ns
done
```

Output:

```
current nscales=1: max residual 0.500/255, share <= 1/255 1.00000, share <= 2/255 1.00000, n=327143
v1.0.1  nscales=1: max residual 256.498/255, share <= 1/255 0.99871, share <= 2/255 0.99871, n=327567
current nscales=4: max residual 0.500/255, share <= 1/255 1.00000, share <= 2/255 1.00000, n=327680
v1.0.1  nscales=4: max residual 0.500/255, share <= 1/255 1.00000, share <= 2/255 1.00000, n=327680
```

For the 1.0.1, `nscales = 1` archive, the residuals above 1/255, from the same script with the last line replaced by a filter on `res > 1/255`:

```
beyond 1/255: 424 pixels, min 255.503/255, max 256.498/255
```

Those are pixels above white that 1.0.1 wrote as dark, and that this development version writes white (`NEWS.md`, "Pixels above white"). The analysis file will state the configurations measured and no others; `@return` cites it rather than restating the numbers.
