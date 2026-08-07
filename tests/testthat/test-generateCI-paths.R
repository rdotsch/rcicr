# Coverage for the two branches of generateCI() that a plain call never reaches:
# the `participants` two-step averaging, and the two z-map methods.
#
# Both were previously untested. Between them they account for most of the
# conditional logic in generateCI(), including a second parallel loop.

# --------------------------------------------------------------------------
# participants: per-participant CIs, then average
# --------------------------------------------------------------------------

test_that("participants averages per-participant CIs rather than pooling trials", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 12, nscales = 1, seed = 1)

  responses <- rep(c(1, -1), 6)
  pid <- rep(c("a", "b", "c"), each = 4)

  group <- generateCI(
    stimuli = 1:12, responses = responses, baseimage = "base", rdata = rdata,
    participants = pid, save_as_png = FALSE, n_cores = 1
  )

  # Compute each participant's CI independently, then average by hand. The whole
  # point of the `participants` argument is that this is *not* the same thing as
  # pooling every trial, whenever participants contribute unequal trial counts.
  per_pid <- lapply(split(seq_along(pid), pid), function(rows) {
    generateCI(
      stimuli = rows, responses = responses[rows], baseimage = "base",
      rdata = rdata, save_as_png = FALSE, n_cores = 1
    )$ci
  })
  by_hand <- Reduce(`+`, per_pid) / length(per_pid)

  expect_equal(group$ci, by_hand)
})

test_that("participants and pooled trials diverge when the design is unbalanced", {
  # Guards the reason the argument exists. With a balanced design the two routes
  # coincide and this test would be vacuous, so the design here is deliberately
  # lopsided: one participant contributes 8 trials, the other 2.
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 10, nscales = 1, seed = 1)

  responses <- c(rep(1, 8), rep(-1, 2))
  pid <- c(rep("a", 8), rep("b", 2))

  by_participant <- generateCI(
    stimuli = 1:10, responses = responses, baseimage = "base", rdata = rdata,
    participants = pid, save_as_png = FALSE, n_cores = 1
  )$ci

  pooled <- generateCI(
    stimuli = 1:10, responses = responses, baseimage = "base", rdata = rdata,
    save_as_png = FALSE, n_cores = 1
  )$ci

  # Averaging participant CIs weights each person equally; pooling weights each
  # trial equally. The heavy contributor dominates the pooled version.
  expect_false(isTRUE(all.equal(by_participant, pooled)))
})

test_that("save_individual_cis writes one PNG per participant", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 12, nscales = 1, seed = 1)

  target <- file.path(tmp, "cis")
  generateCI(
    stimuli = 1:12, responses = rep(c(1, -1), 6), baseimage = "base",
    rdata = rdata, participants = rep(c("a", "b", "c"), each = 4),
    save_individual_cis = TRUE, save_as_png = FALSE, targetpath = target,
    n_cores = 1
  )

  written <- list.files(file.path(target, "individual_cis"), pattern = "\\.png$")
  expect_setequal(written, c("ci_a.png", "ci_b.png", "ci_c.png"))
})

# --------------------------------------------------------------------------
# z-maps
# --------------------------------------------------------------------------

test_that("zmapmethod = 'quick' returns a thresholded z-map", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 12, nscales = 1, seed = 1)

  res <- generateCI(
    stimuli = 1:12, responses = rep(c(1, -1), 6), baseimage = "base",
    rdata = rdata, save_as_png = FALSE, zmap = TRUE, zmapmethod = "quick",
    threshold = 1, zmapdecoration = FALSE, n_cores = 1,
    zmaptargetpath = file.path(tmp, "zmaps")
  )

  expect_true("zmap" %in% names(res))
  expect_equal(dim(res$zmap), c(32, 32))

  # Everything below the threshold is blanked, everything kept is above it.
  kept <- res$zmap[!is.na(res$zmap)]
  expect_true(length(kept) > 0)
  expect_true(all(abs(kept) >= 1))
})

test_that("zmapmethod = 't.test' returns a z-map of per-pixel test statistics", {
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 12, nscales = 1, seed = 1)

  res <- generateCI(
    stimuli = 1:12, responses = rep(c(1, -1), 6), baseimage = "base",
    rdata = rdata, save_as_png = FALSE, zmap = TRUE, zmapmethod = "t.test",
    zmapdecoration = FALSE, n_cores = 1,
    zmaptargetpath = file.path(tmp, "zmaps")
  )

  expect_equal(dim(res$zmap), c(32, 32))
  expect_false(anyNA(res$zmap))

  # This method is not thresholded, and its sign is taken from the CI itself.
  expect_equal(sign(res$zmap), sign(res$ci))
})

test_that("a z-map is written where zmaptargetpath asks for it", {
  # Regression test. zmaptargetpath was documented and accepted but never passed
  # through to plotZmap(), so every z-map landed in ./zmaps relative to the
  # working directory no matter what the caller asked for. Assert the intended
  # behaviour, never the old one.
  tmp <- withr::local_tempdir()
  withr::local_dir(tmp)
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 12, nscales = 1, seed = 1)

  target <- file.path(tmp, "my_zmaps")
  generateCI(
    stimuli = 1:12, responses = rep(c(1, -1), 6), baseimage = "base",
    rdata = rdata, save_as_png = FALSE, zmap = TRUE, threshold = 1,
    zmaptargetpath = target, zmapdecoration = FALSE, n_cores = 1
  )

  expect_true(file.exists(file.path(target, "base.png")))
  expect_false(dir.exists(file.path(tmp, "zmaps")))
})
