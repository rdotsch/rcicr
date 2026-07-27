## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 4,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(rcicr)

## ----generate-base-face-------------------------------------------------------
set.seed(42)
base_face_path <- tempfile(fileext = ".png")
png::writePNG(matrix(runif(64 * 64), 64, 64), base_face_path)

## ----generate-stimuli, results = "hide"---------------------------------------
stimulus_path <- tempdir()

generateStimuli2IFC(
  base_face_files = list(face = base_face_path),
  n_trials        = 20,
  img_size        = 64,
  stimulus_path   = stimulus_path,
  seed            = 1,
  ncores          = 1,
  save_as_png     = FALSE # set to TRUE to also write stimulus PNGs to stimulus_path
)

rdata_file <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)[1]

## ----simulate-responses-------------------------------------------------------
responses <- sample(c(1, -1), 20, replace = TRUE)

## ----compute-ci---------------------------------------------------------------
ci <- generateCI(
  stimuli     = 1:20,
  responses   = responses,
  baseimage   = "face",
  rdata       = rdata_file,
  save_as_png = FALSE
)

names(ci)

## ----plot-ci------------------------------------------------------------------
image(ci$combined, col = gray.colors(256), axes = FALSE, asp = 1)
