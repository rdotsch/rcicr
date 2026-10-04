How well stimulus PNGs identify the settings that made them
================

- [Two installed versions](#two-installed-versions)
- [The archives](#the-archives)
- [The residuals behind 1/255](#the-residuals-behind-1255)
- [Shares with right and wrong
  settings](#shares-with-right-and-wrong-settings)
- [What this shows, and what it does
  not](#what-this-shows-and-what-it-does-not)

`checkStimulusPNGs2IFC()` compares a stimulus `.Rdata` file with
stimulus PNGs. An original and an inverted stimulus are rendered over
the same base image, so wherever neither is clipped, `ori - inv` is the
trial’s noise divided by 0.6. A pixel agrees when its residual is within
1/255.

This document measures two things on synthetic stimulus sets written by
this checkout and by rcicr 1.0.1:

1.  the residual between `ori - inv` and the noise a regenerated file
    records, which sets the 1/255 tolerance;
2.  the share of agreeing pixels per trial, with the right settings and
    with a wrong `nscales`, seed, noise type, Gabor `sigma`, base order
    or `use_same_parameters`.

It measures only the configurations below: 128 pixels, 20 trials, seed
42 and a synthetic base face. It does not establish shares for other
image sizes or base images.

## Two installed versions

Each version is installed from git into a library of its own. Archives
are written by running the version’s own `generateStimuli2IFC()` in a
separate R process. Candidates are regenerated, and scored, by the
checkout loaded in this session.

``` r
repo <- normalizePath("..")
lib_current <- tempfile("lib_current")
lib_v101 <- tempfile("lib_v101")
v101_source <- tempfile("v101")
invisible(lapply(c(lib_current, lib_v101, v101_source), dir.create))
system2("sh", c("-c", shQuote(sprintf("git -C %s archive v1.0.1 | tar -x -C %s",
                                      shQuote(repo), shQuote(v101_source)))))
r_cmd <- file.path(R.home("bin"), "R")
install <- function(lib, source) {
  status <- system2(r_cmd, c("CMD", "INSTALL", "--no-test-load", "-l", shQuote(lib), shQuote(source)),
                    stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0)
}
install(lib_current, repo)
install(lib_v101, v101_source)

# Parallel workers run library(rcicr), and find the checkout through R_LIBS.
Sys.setenv(R_LIBS = lib_current)
.libPaths(c(lib_current, .libPaths()))
library(rcicr)
c(current = format(packageVersion("rcicr")),
  v101 = format(packageVersion("rcicr", lib.loc = lib_v101)))
     current         v101 
"1.5.0.9000"      "1.0.1" 
```

## The archives

Each version writes the same five stimulus sets of a synthetic face,
with PNGs and no `.Rdata` file.

``` r
n <- 128
n_trials <- 20
archive_seed <- 42
x <- matrix(seq(-1, 1, length.out = n), n, n, byrow = TRUE)
y <- -matrix(seq(-1, 1, length.out = n), n, n)
faces <- list(face = exp(-(x^2 / 0.45 + y^2 / 0.75)), other = exp(-(x^2 / 0.75 + y^2 / 0.45)))
faces <- lapply(faces, function(f) (f - min(f)) / diff(range(f)))

archive_settings <- list(
  nscales1 = list(bases = "face", nscales = 1),
  nscales4 = list(bases = "face", nscales = 4),
  sinusoid = list(bases = "face"),
  gabor = list(bases = "face", noise_type = "gabor"),
  two_bases = list(bases = c("face", "other"), use_same_parameters = FALSE)
)

write_archive <- function(lib, settings) {
  dir <- tempfile("archive")
  dir.create(dir)
  base_files <- vapply(settings$bases, function(b) {
    file <- file.path(dir, paste0("base_", b, ".png"))
    png::writePNG(faces[[b]], file)
    file
  }, character(1))
  args <- c(list(base_face_files = as.list(base_files), n_trials = n_trials, img_size = n,
                 stimulus_path = dir, seed = archive_seed, ncores = 1, save_rdata = FALSE),
            settings[setdiff(names(settings), "bases")])
  call_file <- tempfile(fileext = ".rds")
  saveRDS(args, call_file)
  script <- tempfile(fileext = ".R")
  writeLines(c("args <- readRDS(commandArgs(TRUE)[1])",
               "invisible(capture.output(do.call(rcicr::generateStimuli2IFC, args)))"), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(script), shQuote(call_file)),
                    env = paste0("R_LIBS=", shQuote(lib)), stdout = FALSE, stderr = FALSE)
  stopifnot(status == 0)
  dir
}

archives <- list(current = lapply(archive_settings, write_archive, lib = lib_current),
                 v101 = lapply(archive_settings, write_archive, lib = lib_v101))
count_pngs <- function(dir) length(list.files(dir, "_(ori|inv)\\.png$"))
sapply(archives, vapply, count_pngs, integer(1))
          current v101
nscales1       40   40
nscales4       40   40
sinusoid       40   40
gabor          40   40
two_bases      80   80
```

Candidates are what a researcher without the `.Rdata` file regenerates:
the same call over a grey stand-in for every base image, without PNGs
and without contrast maximisation.

``` r
grey <- tempfile(fileext = ".png")
png::writePNG(matrix(0.5, n, n), grey)

candidate <- function(bases = "face", seed = archive_seed, ...) {
  dir <- tempfile("candidate")
  base_files <- stats::setNames(rep(list(grey), length(bases)), bases)
  generate <- function() {
    generateStimuli2IFC(base_files, n_trials = n_trials, img_size = n, stimulus_path = dir,
                        seed = seed, ncores = 1, save_as_png = FALSE,
                        maximize_baseimage_contrast = FALSE, ...)
  }
  invisible(capture.output(suppressMessages(generate())))
  list.files(dir, "\\.Rdata$", full.names = TRUE)
}
```

## The residuals behind 1/255

For every pixel neither image clips, the residual is
`abs(ori - inv - noise / 0.6)`, with the noise from a candidate
regenerated with the right settings.

``` r
residuals_of <- function(archive, rdata, base = "face") {
  stored <- new.env()
  load(rdata, envir = stored)
  unlist(lapply(seq_len(n_trials), function(trial) {
    ori <- png::readPNG(sprintf("%s/rcic_%s_%d_%05d_ori.png", archive, base, archive_seed, trial))
    inv <- png::readPNG(sprintf("%s/rcic_%s_%d_%05d_inv.png", archive, base, archive_seed, trial))
    noise <- generateNoiseImage(stored$stimuli_params[[base]][trial, ], stored$p)
    unclipped <- ori > 0 & ori < 1 & inv > 0 & inv < 1
    abs(ori - inv - noise / 0.6)[unclipped]
  }))
}

residual_rows <- do.call(rbind, lapply(c("current", "v101"), function(version) {
  do.call(rbind, lapply(c(1, 4), function(ns) {
    r <- residuals_of(archives[[version]][[paste0("nscales", ns)]], candidate(nscales = ns)) * 255
    data.frame(version = version, nscales = ns, pixels = length(r),
               max_residual = round(max(r), 3),
               beyond_1 = sum(r > 1 + 1e-9),
               between_0.5_and_255 = sum(r > 0.5 + 1e-9 & r < 255),
               min_beyond_1 = if (any(r > 1 + 1e-9)) round(min(r[r > 1 + 1e-9]), 3) else NA)
  }))
}))
knitr::kable(residual_rows)
```

| version | nscales | pixels | max_residual | beyond_1 | between_0.5_and_255 | min_beyond_1 |
|:--------|--------:|-------:|-------------:|---------:|--------------------:|-------------:|
| current |       1 | 327143 |        0.500 |        0 |                   0 |           NA |
| current |       4 | 327680 |        0.500 |        0 |                   0 |           NA |
| v101    |       1 | 327567 |      256.498 |      424 |                   0 |      255.503 |
| v101    |       4 | 327680 |        0.500 |        0 |                   0 |           NA |

Residuals are in units of 1/255. Where a residual is at most 0.5/255,
`ori - inv` and the regenerated noise differ only by the PNG writer’s
8-bit rounding. Pixels beyond 1/255 are the ones rcicr 1.0.1 wrote dark
when they lay above white, which later versions write white (`NEWS.md`,
“Pixels above white”); no residual falls between 0.5/255 and those. A
tolerance of 1/255 therefore accepts all rounding, with headroom, and
rejects only the overflow pixels.

## Shares with right and wrong settings

`checkStimulusPNGs2IFC()` is given the archive’s label and seed, so that
the PNGs are found whatever the candidate holds. Each row summarises the
per-trial shares of one base label.

``` r
checks <- list(
  list(archive = "nscales1", setting = "right settings", args = list(nscales = 1)),
  list(archive = "sinusoid", setting = "right settings"),
  list(archive = "sinusoid", setting = "nscales = 4", args = list(nscales = 4)),
  list(archive = "sinusoid", setting = "seed = 43", args = list(seed = 43)),
  list(archive = "sinusoid", setting = "noise_type = \"gabor\"", args = list(noise_type = "gabor")),
  list(archive = "gabor", setting = "right settings", args = list(noise_type = "gabor")),
  list(archive = "gabor", setting = "sigma = 24", args = list(noise_type = "gabor", sigma = 24)),
  list(archive = "gabor", setting = "sigma = 20", args = list(noise_type = "gabor", sigma = 20)),
  list(archive = "gabor", setting = "noise_type = \"sinusoid\""),
  list(archive = "two_bases", setting = "right settings", bases = c("face", "other"),
       args = list(use_same_parameters = FALSE)),
  list(archive = "two_bases", setting = "bases in the other order", bases = c("other", "face"),
       args = list(use_same_parameters = FALSE)),
  list(archive = "two_bases", setting = "use_same_parameters = TRUE", bases = c("face", "other"))
)

score <- function(check) {
  bases <- if (is.null(check$bases)) "face" else check$bases
  rdata <- do.call(candidate, c(list(bases = bases), check$args))
  do.call(rbind, lapply(c("current", "v101"), function(version) {
    result <- checkStimulusPNGs2IFC(rdata, archives[[version]][[check$archive]], label = "rcic",
                                seed = archive_seed)
    shares <- split(result$share, factor(result$base, unique(result$base)))
    data.frame(archive = check$archive, candidate = check$setting, version = version,
               base = names(shares),
               min = round(vapply(shares, min, numeric(1)), 4),
               median = round(vapply(shares, stats::median, numeric(1)), 4),
               max = round(vapply(shares, max, numeric(1)), 4),
               trials_below_1 = vapply(shares, function(share) sum(share < 1), integer(1)),
               row.names = NULL)
  }))
}
```

``` r
share_rows <- do.call(rbind, lapply(checks, score))
knitr::kable(share_rows)
```

| archive   | candidate                  | version | base  |    min | median |    max | trials_below_1 |
|:----------|:---------------------------|:--------|:------|-------:|-------:|-------:|---------------:|
| nscales1  | right settings             | current | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| nscales1  | right settings             | v101    | face  | 0.9860 | 1.0000 | 1.0000 |              2 |
| sinusoid  | right settings             | current | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| sinusoid  | right settings             | v101    | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| sinusoid  | nscales = 4                | current | face  | 0.0209 | 0.0237 | 0.0735 |             20 |
| sinusoid  | nscales = 4                | v101    | face  | 0.0209 | 0.0237 | 0.0735 |             20 |
| sinusoid  | seed = 43                  | current | face  | 0.0215 | 0.0258 | 0.0295 |             20 |
| sinusoid  | seed = 43                  | v101    | face  | 0.0215 | 0.0258 | 0.0295 |             20 |
| sinusoid  | noise_type = “gabor”       | current | face  | 0.0246 | 0.0294 | 0.0334 |             20 |
| sinusoid  | noise_type = “gabor”       | v101    | face  | 0.0246 | 0.0294 | 0.0334 |             20 |
| gabor     | right settings             | current | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| gabor     | right settings             | v101    | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| gabor     | sigma = 24                 | current | face  | 0.9500 | 0.9834 | 0.9963 |             20 |
| gabor     | sigma = 24                 | v101    | face  | 0.9500 | 0.9834 | 0.9963 |             20 |
| gabor     | sigma = 20                 | current | face  | 0.3507 | 0.4223 | 0.4990 |             20 |
| gabor     | sigma = 20                 | v101    | face  | 0.3507 | 0.4223 | 0.4990 |             20 |
| gabor     | noise_type = “sinusoid”    | current | face  | 0.0233 | 0.0298 | 0.0328 |             20 |
| gabor     | noise_type = “sinusoid”    | v101    | face  | 0.0233 | 0.0298 | 0.0328 |             20 |
| two_bases | right settings             | current | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| two_bases | right settings             | current | other | 1.0000 | 1.0000 | 1.0000 |              0 |
| two_bases | right settings             | v101    | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| two_bases | right settings             | v101    | other | 1.0000 | 1.0000 | 1.0000 |              0 |
| two_bases | bases in the other order   | current | other | 0.0224 | 0.0250 | 0.0272 |             20 |
| two_bases | bases in the other order   | current | face  | 0.0219 | 0.0253 | 0.0293 |             20 |
| two_bases | bases in the other order   | v101    | other | 0.0224 | 0.0250 | 0.0272 |             20 |
| two_bases | bases in the other order   | v101    | face  | 0.0219 | 0.0253 | 0.0293 |             20 |
| two_bases | use_same_parameters = TRUE | current | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| two_bases | use_same_parameters = TRUE | current | other | 0.0224 | 0.0250 | 0.0272 |             20 |
| two_bases | use_same_parameters = TRUE | v101    | face  | 1.0000 | 1.0000 | 1.0000 |              0 |
| two_bases | use_same_parameters = TRUE | v101    | other | 0.0224 | 0.0250 | 0.0272 |             20 |

## What this shows, and what it does not

In these configurations:

- **The right settings** agree on every compared pixel of every trial,
  with one exception: the 1.0.1 `nscales = 1` archive, whose overflow
  pixels lower some trials’ shares below 1.
- **A wrong discrete setting** (`nscales`, seed, noise type, base order,
  `use_same_parameters`) agrees on 2.1% to 7.4% of pixels in every trial
  of each base label it changes.
- **A wrong Gabor `sigma` close to the right one** agrees on most
  pixels: `sigma = 24` against an archive made with 25 agrees on 95% or
  more in every trial, though on no trial completely. `sigma` is
  continuous, so a share near 1 does not by itself rule out a nearby
  value.
- **No fixed threshold separates right from wrong across archives.** The
  right settings scored as low as 0.986 on the 1.0.1 overflow archive,
  below the 0.996 that `sigma = 24` reached on another. On the same
  archive, the right candidate scored higher than every wrong one on
  every trial of each base label the wrong setting changes.
- PNGs from both versions give the same shares in every row except the
  `nscales = 1` one.

Two rows show why every base label is checked. With
`use_same_parameters` wrong, the first base draws the same parameters
either way and matches; only the second base shows the mismatch. With
the bases in the other order, both labels are scored against the other
base’s noise.

None of these ranges is a bound for other image sizes, base images or
settings, and `checkStimulusPNGs2IFC()` returns shares, not a verdict.
