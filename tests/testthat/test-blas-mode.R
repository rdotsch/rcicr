test_that("the check log says whether norms and InfoVal were compared exactly", {
  if (!reference_blas()) {
    skip(paste("non-reference BLAS or LAPACK: exact comparisons of norms and InfoVal",
               "relaxed to 1e-12 (", extSoftVersion()[["BLAS"]], ")"))
  }
  expect_true(reference_blas())
})
