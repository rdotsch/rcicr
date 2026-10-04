# What the parallel loops send their workers (#405): the basis's patch indices
# as integer, each participant's rows with its own task, and no more workers
# than tasks. None of it may move a number.

test_that("renderingBasis() holds patch indices as integer and renders identically", {
  p <- generateNoisePattern(32, nscales = 2)
  expect_type(p$patchIdx, "double")
  basis <- rcicr:::renderingBasis(p)
  expect_type(basis$patchIdx, "integer")
  expect_identical(dim(basis$patchIdx), dim(p$patchIdx))
  set.seed(3)
  params <- runif(max(p$patchIdx), -1, 1)
  expect_identical(generateNoiseImage(params, basis), generateNoiseImage(params, p))

  legacy <- list(sinusoids = p$patches, sinIdx = p$patchIdx)
  expect_type(rcicr:::renderingBasis(legacy)$sinIdx, "integer")

  # A double index selects as if truncated, so a fractional one renders the same.
  fractional <- p
  fractional$patchIdx[1] <- 1.5
  expect_identical(generateNoiseImage(params, rcicr:::renderingBasis(fractional)),
                   generateNoiseImage(params, fractional))
})

test_that("the stimulus file keeps the basis as it was built", {
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, img_size = 32, n_trials = 4)
  stored <- new.env()
  load(rdata, envir = stored)
  expect_type(stored$p$patchIdx, "double")
})

test_that("startBackend() starts no more workers than there are tasks", {
  expect_null(rcicr:::startBackend(4L, 1L))
  skip_on_cran()
  cl <- rcicr:::startBackend(4L, 2L)
  on.exit(rcicr:::stopClusterSafely(cl), add = TRUE)
  expect_length(cl, 2L)
})

test_that("participants' own rows give identical CIs and PNGs on one core and two", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, img_size = 32, n_trials = 9)
  responses <- c(1, -1, 1, -1, 1, -1, 1, -1, 1)
  # Participant 3 has a single trial, the parameter-vector path.
  pids <- c(1, 1, 1, 1, 2, 2, 2, 2, 3)
  run <- function(n_cores) {
    out <- file.path(dir, paste0("cores", n_cores))
    dir.create(out)
    ci <- generateCI(1:9, responses, "base", rdata, save_as_png = FALSE, participants = pids,
                     save_individual_cis = TRUE, targetpath = out, n_cores = n_cores)
    pngs <- list.files(file.path(out, "individual_cis"), full.names = TRUE)
    list(ci = ci, pngs = basename(pngs), pixels = lapply(pngs, png::readPNG))
  }
  serial <- run(1)
  parallel_run <- run(2)
  expect_identical(serial$ci, parallel_run$ci)
  expect_length(serial$pngs, 3L)
  expect_identical(serial$pngs, parallel_run$pngs)
  expect_identical(serial$pixels, parallel_run$pixels)
})

test_that("the t-test z-map is identical on one core and two", {
  skip_on_cran()
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, img_size = 32, n_trials = 8)
  responses <- rep(c(1, -1), 4)
  run <- function(n_cores) {
    generateCI(1:8, responses, "base", rdata, save_as_png = FALSE, zmap = TRUE,
               zmapmethod = "t.test", zmapdecoration = FALSE,
               zmaptargetpath = file.path(dir, paste0("z", n_cores)), n_cores = n_cores)$zmap
  }
  serial <- run(1)
  expect_false(all(is.na(serial)))
  expect_identical(serial, run(2))
})
