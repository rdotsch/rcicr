test_that("matlab is not imported wholesale, so sum() keeps base semantics", {
  ns <- asNamespace("rcicr")

  expect_false("matlab" %in% names(getNamespaceImports(ns)))

  # The hazard itself: matlab::sum() returns column sums for a matrix rather
  # than a single total, so an @import matlab anywhere in R/ silently changes
  # what an unqualified sum() means throughout the package.
  expect_identical(get("sum", envir = ns), base::sum)
})
