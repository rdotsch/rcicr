test_that("computeInfoVal2IFC uses a pre-seeded reference distribution and never calls yesno()", {
  # ref_lookup in R/computeInfoVal2IFC.R has all its data rows commented out, so
  # the yesno() interactive-prompt branch is currently unreachable. Pre-seeding
  # reference_norms takes the function down the "already have a reference
  # distribution" path; mocking yesno() here is a defensive regression guard in
  # case ref_lookup gets populated again in the future without re-checking this.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  seed_reference_norms(rdata_path, n = 50, seed = 1)

  testthat::local_mocked_bindings(
    yesno = function(...) stop("yesno() should not be called"),
    .package = "yesno"
  )

  target_ci <- list(ci = matrix(withr::with_seed(2, rnorm(32 * 32)), 32, 32))
  iv <- computeInfoVal2IFC(target_ci = target_ci, rdata = rdata_path)

  expect_type(iv, "double")
  expect_length(iv, 1)

  e <- new.env()
  load(rdata_path, envir = e)
  expected <- (norm(matrix(target_ci$ci), "f") - median(e$reference_norms)) / mad(e$reference_norms)
  expect_equal(iv, expected)
})

test_that("response_seed regenerates the null even when the .Rdata already has one, and never caches it", {
  # The trap this guards. computeInfoVal2IFC() only calls the generator when
  # reference_norms is absent or force_gen_ref_dist is TRUE -- and the first
  # call writes reference_norms into the file. So a response_seed that was
  # merely passed through would be accepted, documented, and silently ignored on
  # every call after the first. That is the same shape as force_gen_ref_dist
  # being ignored in this exact branch, and as plotZmap()'s documented-but-
  # unapplied `mask`. Asserting the seeded result *differs* from the cached one
  # is what makes this test non-vacuous.
  tmp <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  seed_reference_norms(rdata_path, n = 50, seed = 1)

  cached <- new.env()
  load(rdata_path, envir = cached)

  target_ci <- list(ci = matrix(withr::with_seed(2, rnorm(32 * 32)), 32, 32))

  iv_cached <- computeInfoVal2IFC(target_ci = target_ci, rdata = rdata_path)
  iv_seeded <- suppressWarnings(computeInfoVal2IFC(
    target_ci = target_ci, rdata = rdata_path, iter = 8, response_seed = 99
  ))

  expect_false(identical(iv_cached, iv_seeded))

  # The one-off check must not have become this stimulus set's reference
  # distribution: the file's norms are untouched, and a later default call
  # returns exactly what it did before.
  after <- new.env()
  load(rdata_path, envir = after)
  expect_identical(after$reference_norms, cached$reference_norms)
  expect_equal(computeInfoVal2IFC(target_ci = target_ci, rdata = rdata_path), iv_cached)
})
