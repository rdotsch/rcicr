# The participant CI loop must return the same numbers on two cores as on one.
#
# test-parallel-equivalence.R looks like it covers this and does not: its title
# says "identical stimuli and CIs" and a comment says it checks the CI "since
# generateCI() has its own parallel loop", but its two generateCI() calls set
# neither n_cores nor participants, so they take the single-shot
# generateCINoise() path, which has no foreach loop in it at all. The ncores
# variation there reaches generateStimuli2IFC() only.

test_that("the participant loop gives identical CIs on one core and two", {
  skip_on_cran()
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, img_size = 32, n_trials = 12)
  responses <- rep(c(1, -1), 6)
  pids <- rep(c(1, 2, 3), each = 4)

  run <- function(n_cores) {
    generateCI(1:12, responses, "base", rdata, save_as_png = FALSE,
               participants = pids, n_cores = n_cores)
  }

  serial <- run(1)
  parallel_run <- run(2)

  # Exact, not tolerance-based: neither loop in this package draws random
  # numbers -- stimulus parameters are drawn under set.seed() before the loop
  # starts -- so there is no worker RNG stream to diverge. If this ever needs a
  # tolerance, that assumption has stopped holding and is the thing to look at.
  expect_identical(serial$ci, parallel_run$ci)
})

test_that("the per-participant stack survives the parallel path in order", {
  skip_on_cran()
  skip_if_not_installed("withr")

  # Separate from the group CI on purpose. The group CI is a *mean across
  # participants*, so it is invariant to their order: a .combine that returned
  # the participants shuffled would leave it bit-identical and only show up
  # here. This is also what the t.test z-map consumes, via the
  # `noiseimages <- pid.cis` short-circuit, so an ordering change would reach a
  # researcher's z-map while every group-level check stayed green.
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir, img_size = 32, n_trials = 12)
  responses <- rep(c(1, -1), 6)
  pids <- rep(c(1, 2, 3), each = 4)

  stack <- function(n_cores) {
    e <- new.env()
    load(rdata, envir = e)
    out <- rcicr:::computeParticipantCIs(
      params = e$stimuli_params$base[1:12, ], responses = responses,
      participants = pids, p = e$p, base = e$base_faces$base,
      baseimage = "base", img_size = 32, mask = NA, n_cores = n_cores,
      save_individual_cis = FALSE, targetpath = NULL,
      individual_scaling = "independent", individual_scaling_constant = 0.1,
      antiCI = FALSE
    )
    out$pid_cis
  }

  serial <- stack(1)
  parallel_stack <- stack(2)

  expect_identical(dim(serial), c(32L, 32L, 3L))
  expect_identical(serial, parallel_stack)

  # Not vacuous: the participants must actually differ from one another, or an
  # ordering bug would be invisible to the comparison above.
  expect_false(isTRUE(all.equal(serial[, , 1], serial[, , 2])))
})
