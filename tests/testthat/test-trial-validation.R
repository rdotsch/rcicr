test_that("participant IDs must align with every trial before PNG output", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, n_trials = 20)
  for (n in c(2, 3, 21)) {
    for (container in list(identity, function(x) data.frame(pid = x))) {
      target <- file.path(dir, paste0("invalid-", n))
      expect_error(generateCI(1:20, rep(c(1, -1), 10), "base", rdata,
                     participants = container(rep(c("a", "b"), length.out = n)),
                     save_as_png = TRUE, save_individual_cis = TRUE,
                     targetpath = target, n_cores = 1
                   ), "participants must have one ID per trial")
      expect_false(dir.exists(target))
    }
  }
})

test_that("CI and cumulative APIs reject malformed IDs before aggregation or indexing", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, n_trials = 6)
  bad_ids <- list(c(1.2, 2.8), c(0, 1, 2), c(-1, -2), c(1, NA),
                  c(1, Inf), c(1, -Inf), c(1, NaN), c(1, 7), numeric(),
                  factor(c("1", "2"), levels = c("2", "1")), c("1", "2"), c(TRUE, FALSE))
  for (ids in bad_ids) {
    for (container in list(identity, function(x) data.frame(id = x))) {
      input <- container(ids)
      responses <- rep(1, length(ids))
      expect_error(generateCI(input, responses, "base", rdata,
                              save_as_png = FALSE), "stimuli")
      expect_error(computeCumulativeCICorrelation(input, responses,
                                                  "base", rdata), "stimuli")
    }
  }
  expect_error(generateCI2IFC(c(0, 1, 2), c(1, -1, 1), "base", rdata,
                              save_as_png = FALSE), "stimuli")
  expect_error(computeCumulativeCICorrelation(1:3, c(1, -1), "base", rdata),
               "same length")
})

test_that("parameter selection rejects IDs even when called directly", {
  params <- list(base = matrix(seq_len(36), nrow = 3))
  for (ids in list(c(0, 1, 2), c(1.2, 2.8), c(1, NA), c(1, 4),
                   factor(c("1", "2")), c("1", "2"), c(TRUE, FALSE), numeric())) {
    expect_error(rcicr:::selectStimulusParams(params, "base", ids), "stimuli")
  }
})

test_that("numeric tibbles survive coercion but mixed factor and logical lists do not", {
  good <- rcicr:::coerceTrialVectors(tibble::tibble(id = c(3, 1, 3)),
                                     tibble::tibble(response = c(1, -1, 1)), tibble::tibble(pid = c("b", "a", "b")))
  expect_identical(good$stimuli, c(3, 1, 3))
  expect_identical(good$responses, c(1, -1, 1))
  expect_identical(good$participants, c("b", "a", "b"))
  for (ids in list(list(1, factor("2")), list(1, TRUE),
                   tibble::tibble(id = factor(c("1", "2"))))) {
    expect_error(rcicr:::coerceTrialVectors(ids, c(1, -1), NA), "stimuli")
  }
})

test_that("direct CI noise calls reject response recycling", {
  p <- generateNoisePattern(img_size = 32, nscales = 1)
  params <- matrix(seq_len(2 * max(p$patchIdx)) / 100, nrow = 2)
  for (responses in list(1, c(1, -1, 1), numeric())) {
    expect_error(generateCINoise(params, responses, p), "one parameter row per response")
  }
  expect_error(generateCINoise(params[1, ], c(1, -1), p), "single-trial vector")
  expect_error(generateCINoise(params[FALSE, ], numeric(), p), "one parameter row")
  expect_identical(generateCINoise(params[1, ], -1, p),
                   generateCINoise(params[1, , drop = FALSE], -1, p))
  expected <- generateNoiseImage((params[1, ] - params[2, ]) / 2, p)
  expect_identical(generateCINoise(params, c(1, -1), p), expected)
})

test_that("valid repeated and nonconsecutive trials keep their established weights", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, n_trials = 6)
  e <- new.env()
  load(rdata, envir = e)
  ids <- c(6, 2, 6, 1)
  responses <- c(1, -1, -1, 1)
  # CI weights unique stimuli; the cumulative API deliberately weights trials.
  ci_expected <- generateNoiseImage(
                                    (e$stimuli_params$base[1, ] - e$stimuli_params$base[2, ]) / 3, e$p)
  run <- function(participants = NA) {
    generateCI(data.frame(id = ids), data.frame(response = responses), "base", rdata,
               participants = participants, save_as_png = FALSE, scaling = "none")
  }
  result <- run()
  expect_equal(result$ci, ci_expected, tolerance = 1e-14)
  for (sentinel in list(NA, rep(NA, 4), rep(NA, 2))) {
    expect_identical(run(sentinel), result)
  }
  final <- generateNoiseImage(colMeans(e$stimuli_params$base[ids, ] * responses), e$p)
  expected <- vapply(seq_along(ids), function(i) {
    selected <- seq_len(i)
    partial <- generateNoiseImage(colMeans(
                                           e$stimuli_params$base[ids[selected], , drop = FALSE] * responses[selected]), e$p)
    cor(as.vector(partial), as.vector(final))
  }, numeric(1))
  actual <- computeCumulativeCICorrelation(data.frame(id = ids),
                                           data.frame(response = responses), "base", rdata)
  expect_equal(actual, expected, tolerance = 1e-14)
  expect_false(isTRUE(all.equal(final, ci_expected)))
  single <- generateCI(6, -1, "base", rdata, save_as_png = FALSE, scaling = "none")
  expect_identical(single$ci, generateNoiseImage(-e$stimuli_params$base[6, ], e$p))
  expect_equal(computeCumulativeCICorrelation(6, -1, "base", rdata), 1)
})

test_that("unequal participant groups retain individual PNGs and group values", {
  skip_on_cran()
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, n_trials = 20)
  responses <- rep(c(1, -1), 10)
  ids <- c(rep("b", 4), rep("a", 16))
  e <- new.env()
  load(rdata, envir = e)
  expected <- lapply(c("a", "b"), function(id) {
    rows <- which(ids == id)
    generateNoiseImage(colMeans(e$stimuli_params$base[rows, ] * responses[rows]), e$p)
  })
  group <- (expected[[1]] + expected[[2]]) / 2
  wrong <- generateCI(1:20, responses, "base", rdata,
                      participants = rep(c("a", "b"), 10), save_as_png = FALSE,
                      scaling = "none", n_cores = 1)$ci
  expect_false(isTRUE(all.equal(group, wrong)))
  for (cores in c(1, 2)) {
    target <- file.path(dir, paste0("valid-", cores))
    out <- generateCI(1:20, responses, "base", rdata,
                      participants = data.frame(pid = ids), save_as_png = FALSE,
                      save_individual_cis = TRUE, targetpath = target, scaling = "none",
                      individual_scaling = "constant", individual_scaling_constant = 100,
                      n_cores = cores)
    expect_equal(out$ci, group, tolerance = 1e-14)
    for (i in seq_along(expected)) {
      reference <- file.path(dir, paste0("expected-", i, ".png"))
      png::writePNG(((expected[[i]] + 100) / 200 + e$base_faces$base) / 2, reference)
      actual <- file.path(target, "individual_cis", paste0("ci_", c("a", "b")[i], ".png"))
      expect_identical(png::readPNG(actual), png::readPNG(reference))
    }
    expect_error(generateCI(1:20, responses, "base", rdata,
                            participants = c("a", "b"), save_as_png = FALSE, n_cores = cores),
                 "participants must have one ID per trial")
  }
})

test_that("valid trial selection still reads released stimulus files", {
  for (version in c("1.0.1", "1.0.1-gabor", "1.1.0")) {
    path <- test_path("fixtures", paste0("legacy-rdata-", version, ".Rdata"))
    e <- new.env()
    load(path, envir = e)
    expected <- generateNoiseImage(
                                   (e$stimuli_params$base[1, ] - e$stimuli_params$base[2, ]) / 3, e$p)
    actual <- generateCI(c(6, 2, 6, 1), c(1, -1, -1, 1), "base", path,
                         save_as_png = FALSE, scaling = "none")
    expect_equal(actual$ci, expected, tolerance = 1e-14, info = version)
    single <- generateCI(2, -1, "base", path, save_as_png = FALSE, scaling = "none")
    expect_identical(single$ci, generateNoiseImage(-e$stimuli_params$base[2, ], e$p))
  }
})
