# The t.test z-map's noise-image loop must return the same numbers on two cores
# as on one.
#
# This is the only path that runs that loop: with `participants` given the z-map
# short-circuits to the per-participant CIs instead of building a noise image
# per trial. Every existing exercise of it forces one core --
# test-regression-baseline.R:105-119 and tools/compare-harness.R:295-308 both
# pass n_cores = 1 -- and test-parallel-progress.R runs two but asserts only
# that a progress bar appeared.

test_that("the t.test z-map is identical on one core and two", {
  skip_on_cran()
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, img_size = 32, n_trials = 12)
  responses <- rep(c(1, -1), 6)

  zmap_at <- function(n_cores) {
    generateCI(1:12, responses, "base", rdata, save_as_png = FALSE,
               zmap = TRUE, zmapmethod = "t.test", zmapdecoration = FALSE,
               zmaptargetpath = file.path(dir, paste0("z", n_cores)),
               threshold = 0, n_cores = n_cores)$zmap
  }

  serial <- zmap_at(1)
  parallel_zmap <- zmap_at(2)

  # Exact, not tolerance-based: neither loop in this package draws random
  # numbers, so there is no worker RNG stream to diverge from the sequential
  # one. A tolerance becoming necessary would mean that stopped being true.
  expect_identical(serial, parallel_zmap)

  # Not vacuous: a constant or all-NA z-map would compare equal trivially.
  expect_equal(dim(serial), c(32L, 32L))
  expect_gt(length(unique(as.vector(serial))), 1)
  expect_false(anyNA(serial))
})

test_that("computeZmapTTest reuses the participant CIs rather than the trials", {
  skip_if_not_installed("withr")

  # The cross-branch dependency, now an argument. Passing a pid_cis stack must
  # produce the z-map that stack implies, and must not consult params/responses
  # at all -- so deliberately wrong ones are passed here. Stage 0's end-to-end
  # pin (test-generateCI-paths.R) covers the same switch from the outside.
  set.seed(1)
  stack <- array(rnorm(32 * 32 * 4, mean = 0.02, sd = 0.01), c(32, 32, 4))
  ci <- apply(stack, c(1, 2), mean)

  from_stack <- rcicr:::computeZmapTTest(
    ci = ci, params = NULL, responses = NULL, p = NULL, pid_cis = stack,
    img_size = 32, n_cores = 1
  )

  pmap <- apply(stack, 1:2, function(x) unlist(t.test(x)['p.value']))
  expect_equal(from_stack, sign(ci) * abs(qnorm(pmap / 2)), ignore_attr = TRUE)
})
