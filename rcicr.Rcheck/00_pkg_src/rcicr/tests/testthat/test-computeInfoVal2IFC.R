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
