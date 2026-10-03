# Every reference path makes the same checks (#377): iter before anything is
# simulated, and a stored reference before it is scored against.

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

test_that("an invalid iter stops the default path before anything is simulated", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  before <- tools::md5sum(rdata)
  local_mocked_bindings(referenceNorms = function(...) stop("simulated"))
  for (iter in list(2.5, 0, -1, NA, c(10, 20), "10")) {
    expect_error(generateReferenceDistribution2IFC(rdata, iter = iter, ncores = 1),
                 "iter must be a positive integer", info = format(iter))
  }
  expect_identical(tools::md5sum(rdata), before)
  # A valid iter reaches the simulation, so the check is not what stopped the others.
  expect_error(quietly(generateReferenceDistribution2IFC(rdata, iter = 3, ncores = 1)), "simulated")
})

test_that("a stored default reference that is not finite, or has a MAD of 0, stops", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  ci <- quietly(generateCI(1:8, c(1, -1, -1, 1, 1, 1, -1, 1), "base", rdata, save_as_png = FALSE))
  plant <- function(norms) {
    e <- new.env()
    load(rdata, envir = e)
    e$reference_norms <- norms
    e$reference_norms_seed <- 5
    save(list = ls(e), file = rdata, envir = e)
  }
  plant(c(1, NA, 2))
  expect_error(quietly(computeInfoVal2IFC(ci, rdata)),
               "Invalid cached reference in the stimulus file")
  plant(c(1, 2, 2, 2, 3))
  expect_error(quietly(computeInfoVal2IFC(ci, rdata)), "has a MAD of 0")
  plant(c(1, 2, 3, 4, 5))
  norm_ci <- norm(matrix(ci$ci), "f")
  expect_identical(quietly(computeInfoVal2IFC(ci, rdata)), (norm_ci - 3) / mad(1:5))
})

test_that("a per-base reference with a MAD of 0 stops instead of giving Inf", {
  independent <- make_independent_fixture(withr::local_tempdir())
  ci <- quietly(generateCI(1:12, rep(c(1, -1), 6), "second", independent, save_as_png = FALSE))
  e <- new.env()
  load(independent, envir = e)
  e$reference_norms_by_base <- list(second = list(norms = c(4, 4, 4, 5), response_seed = 5))
  save(list = ls(e), file = independent, envir = e)
  expect_error(quietly(computeInfoVal2IFC(ci, independent, baseimage = "second")),
               "The reference for baseimage second has a MAD of 0")
})
