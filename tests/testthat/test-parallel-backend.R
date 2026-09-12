# parallel renders a worker's connection arguments into a command line as text,
# so the parent's formatting options decide what the child parses. Under
# OutDec = "," with scipen forcing scientific notation the port was written
# 1,1e+04 and read as NA, and every cluster start failed -- after the full
# connection timeout, four minutes at the default (issue #316).
#
# The assertions are the fix rather than the failure: reproducing the stall
# costs those four minutes on every run. A cluster that starts under those
# options is the property that was broken, and the caller's options surviving
# the call is what makes neutralising them safe.
#
# CI runners have OutDec = ".", so nothing else in the suite can hold this.

test_that("startBackend starts a cluster under any OutDec and scipen", {
  withr::local_options(OutDec = ",", scipen = -10)

  cl <- rcicr:::startBackend(2L)
  on.exit(rcicr:::stopClusterSafely(cl), add = TRUE)

  expect_false(is.null(cl))
  expect_length(cl, 2L)

  # A cluster that exists but cannot carry work would satisfy the above.
  expect_identical(unlist(parallel::clusterEvalQ(cl, 1L + 1L)), c(2L, 2L))
})

test_that("startBackend leaves the caller's formatting options alone", {
  withr::local_options(OutDec = ",", scipen = -10)

  cl <- rcicr:::startBackend(2L)
  rcicr:::stopClusterSafely(cl)

  expect_identical(getOption("OutDec"), ",")
  expect_identical(getOption("scipen"), -10)
})
