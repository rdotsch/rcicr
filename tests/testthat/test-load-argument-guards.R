# Every function that reads a stimulus set does load(rdata), and load() assigns
# straight into the calling function's frame -- so an object in the file silently
# overwrites an argument of the same name.
#
# This has bitten twice, both times from the .Rdata side: a field was added to
# the file and captured an argument that had been there all along (the noise
# `sigma` took over generateCI()'s z-map blur `sigma`; fixed in #146). No field
# written today collides with the arguments below, so these tests pass against
# the unguarded code too if the collision is not planted -- which is exactly why
# each one plants it explicitly.

# Planting is done with mutate_rdata() from helper-fixtures.R, which the
# "did not contain X" tests in test-error-paths.R also use, in its removing
# direction.

test_that("computeInfoVal2IFC scores the caller's CI, not one found in the .Rdata", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  seed_reference_norms(rdata)

  mine  <- list(ci = matrix(0.10, 32, 32))
  decoy <- list(ci = matrix(0.90, 32, 32))

  # target_ci is read at the very end, to compute the CI norm -- and it is read
  # after a *second* load() on the cache-writing path, so both loads need the
  # restore. A file carrying this name would return a plausible number for the
  # wrong image rather than an error, which is the worst shape this can take.
  mutate_rdata(rdata, target_ci = decoy)

  norms <- local({
    e <- new.env()
    load(rdata, envir = e)
    e$reference_norms
  })
  expected_mine  <- (norm(matrix(mine$ci), "f")  - median(norms)) / mad(norms)
  expected_decoy <- (norm(matrix(decoy$ci), "f") - median(norms)) / mad(norms)

  val <- suppressWarnings(computeInfoVal2IFC(target_ci = mine, rdata = rdata))

  expect_equal(val, expected_mine)

  # Discriminating half: the two CIs must actually score differently, or the
  # assertion above would hold whichever image were used.
  expect_false(isTRUE(all.equal(expected_mine, expected_decoy)))
  expect_false(isTRUE(all.equal(val, expected_decoy)))
})

test_that("computeInfoVal2IFC keeps its own rdata path across the load", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)
  seed_reference_norms(rdata)

  # An older generateReferenceDistribution2IFC() saved its own `rdata` argument
  # into the file, so this name really does occur in the wild.
  mutate_rdata(rdata, rdata = file.path(dir, "does-not-exist.Rdata"))

  expect_no_error(
    suppressWarnings(computeInfoVal2IFC(target_ci = list(ci = matrix(0.1, 32, 32)),
                                        rdata = rdata))
  )
})

test_that("computeCumulativeCICorrelation keeps its arguments across the load", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  responses <- rep(c(1, -1), 3)

  reference <- suppressWarnings(computeCumulativeCICorrelation(
    stimuli = 1:6, responses = responses, baseimage = "base", rdata = rdata
  ))

  # Plant collisions for three arguments at once: the step size, the label used
  # to look up the parameters, and the stimulus indices themselves.
  mutate_rdata(rdata, step = 99L, baseimage = "not-a-real-label",
               stimuli = integer(0))

  planted <- suppressWarnings(computeCumulativeCICorrelation(
    stimuli = 1:6, responses = responses, baseimage = "base", rdata = rdata
  ))

  expect_equal(planted, reference)

  # Without the guard, `baseimage` alone would take the lookup to a label that
  # does not exist and error out, so this cannot pass vacuously.
  expect_gt(length(reference), 0)
})
