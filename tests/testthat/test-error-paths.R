# Backlog item 37. The suite was strong on the numeric paths and thin on the
# failure side: 9 assertions against 33 stop()/warning() calls. An unexercised
# guard is indistinguishable from one that works, which is how items 6, 23 and
# 28 each stayed broken for years -- and item 38 was a guard nothing ran that
# pasted 8,190 characters of pixel values into its own message.
#
# Two conventions here, both deliberate:
#
#  * Match a short load-bearing fragment of each message, never the whole
#    sentence. Item 14 will rewrite this wording, and tests that pin it turn
#    that into a slog.
#  * Assert nothing platform-dependent. These tests ship in the tarball and run
#    on CRAN's check farm. Our own message text is stable; locale-formatted
#    numbers, paths and anything read back from a graphics device are not.

# --------------------------------------------------------------------------
# generateCI(): the argument the user is most likely to get wrong
# --------------------------------------------------------------------------

test_that("generateCI rejects stimuli and responses of different lengths", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  err <- expect_error(
    generateCI(stimuli = 1:6, responses = c(1, -1, 1), baseimage = "base",
               rdata = rdata, save_as_png = FALSE),
    "same length"
  )

  # Reporting both lengths is what makes this message actionable -- knowing
  # they differ is not enough to find which vector is wrong.
  msg <- conditionMessage(err)
  expect_match(msg, "6")
  expect_match(msg, "3")

  # And it fires before the file is read, so a mismatch is reported even when
  # the rdata path is nonsense.
  expect_error(
    generateCI(stimuli = 1:6, responses = c(1, -1, 1), baseimage = "base",
               rdata = file.path(dir, "no-such-file.Rdata"), save_as_png = FALSE),
    "same length"
  )
})

# --------------------------------------------------------------------------
# generateCI(): the four ".Rdata did not contain X" guards
#
# These are the returning-user path the CRAN reinstatement is for: a stimulus
# set generated years ago, by a version that saved a different set of names.
# Reached by removing the name from a good fixture rather than by hand-building
# a broken one, so the file stays realistic in every other respect.
# --------------------------------------------------------------------------

test_that("generateCI reports which variable an .Rdata file is missing", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()

  expect_missing <- function(removed, pattern) {
    d <- withr::local_tempdir()
    rdata <- make_fixture_rdata(d)
    mutate_rdata(rdata, .remove = removed)
    expect_error(
      generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3), baseimage = "base",
                 rdata = rdata, save_as_png = FALSE),
      pattern
    )
  }

  # The noise basis. The guard is `!exists('s') & !exists('p')`, so both names
  # have to go -- the fixture carries p, and s only exists in pre-0.3.3 files.
  expect_missing(c("s", "p"), "did not contain s or p")

  expect_missing("base_faces", "did not contain base_faces")
  expect_missing("stimuli_params", "did not contain stimuli_params")

  # img_size is not a formal of generateCI(), so it can only come from the file
  # -- unlike sigma, which is both, and which is what item 32 turned on.
  expect_missing("img_size", "did not contain img_size")
})

test_that("generateCI names the unknown base image and lists the ones it has", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  err <- expect_error(
    generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
               baseimage = "no-such-label", rdata = rdata, save_as_png = FALSE),
    "no-such-label"
  )

  # Listing what the file does contain is the difference between a dead end and
  # a one-line fix, since the label is a key the user chose at generation time.
  expect_match(conditionMessage(err), "base")
})

# --------------------------------------------------------------------------
# computeCumulativeCICorrelation(): the same guards, none of them covered
# --------------------------------------------------------------------------

test_that("computeCumulativeCICorrelation reports which variable is missing", {
  skip_if_not_installed("withr")

  expect_missing <- function(removed, pattern) {
    d <- withr::local_tempdir()
    rdata <- make_fixture_rdata(d)
    mutate_rdata(rdata, .remove = removed)
    expect_error(
      computeCumulativeCICorrelation(stimuli = 1:6, responses = rep(c(1, -1), 3),
                                     baseimage = "base", rdata = rdata),
      pattern
    )
  }

  expect_missing(c("s", "p"), "did not contain s or p")
  expect_missing("base_faces", "did not contain base_faces")
  expect_missing("stimuli_params", "did not contain stimuli_params")
})

test_that("computeCumulativeCICorrelation names the unknown base image", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  expect_error(
    computeCumulativeCICorrelation(stimuli = 1:6, responses = rep(c(1, -1), 3),
                                   baseimage = "no-such-label", rdata = rdata),
    "no-such-label"
  )
})

# --------------------------------------------------------------------------
# applyMask(): all four import failures
#
# plotZmap() has one of these covered; generateCI() has none -- and generateCI()
# is the copy that has always been live. Called directly, as test-plotZmap.R
# already does, because reaching them through generateCI() would test the
# wiring rather than the guard.
# --------------------------------------------------------------------------

test_that("applyMask rejects a mask that is neither a string nor a matrix", {
  expect_error(rcicr:::applyMask(matrix(1, 8, 8), mask = TRUE, img_size = 8),
               "neither a string nor a matrix")
  expect_error(rcicr:::applyMask(matrix(1, 8, 8), mask = list(), img_size = 8),
               "neither a string nor a matrix")

  # A vector is a double but has no dim, so it takes the same branch.
  expect_error(rcicr:::applyMask(matrix(1, 8, 8), mask = rep(0, 64), img_size = 8),
               "neither a string nor a matrix")
})

test_that("applyMask rejects a mask of the wrong size and says both sizes", {
  err <- expect_error(
    rcicr:::applyMask(matrix(1, 8, 8), mask = matrix(0, 4, 4), img_size = 8),
    "same dimensions"
  )

  msg <- conditionMessage(err)
  expect_match(msg, "8")
  expect_match(msg, "4")
})

test_that("applyMask rejects a mask that is not binary", {
  mask <- matrix(0, 8, 8)
  mask[1, 1] <- 0.5

  expect_error(rcicr:::applyMask(matrix(1, 8, 8), mask = mask, img_size = 8),
               "other than 0 or 1")
})

test_that("applyMask rejects a PNG whose colour channels genuinely differ", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  path <- file.path(dir, "rgb.png")

  # Three channels that are not identical, so the greyscale conversion cannot
  # apply. The neighbouring case -- greyscale stored as RGB -- must still work,
  # and is asserted below, because this stop() once ran unconditionally and
  # rejected convertible images too.
  rgb <- array(0, dim = c(8, 8, 3))
  rgb[, , 1] <- 1
  png::writePNG(rgb, path)

  expect_error(rcicr:::applyMask(matrix(1, 8, 8), mask = path, img_size = 8),
               "not a greyscale image")

  grey_as_rgb <- array(0, dim = c(8, 8, 3))
  png::writePNG(grey_as_rgb, file.path(dir, "grey.png"))
  expect_no_error(
    rcicr:::applyMask(matrix(1, 8, 8), mask = file.path(dir, "grey.png"),
                      img_size = 8)
  )
})

# --------------------------------------------------------------------------
# generateStimuli2IFC(): base image import
# --------------------------------------------------------------------------

test_that("generateStimuli2IFC rejects a base image that is not a PNG or JPEG", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  path <- file.path(dir, "base.txt")
  writeLines("not an image", path)

  expect_error(
    generateStimuli2IFC(base_face_files = list(base = path), n_trials = 2,
                        img_size = 32, stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE),
    "PNG or JPEG"
  )
})

test_that("generateStimuli2IFC rejects a non-square base image and reports its size", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  path <- file.path(dir, "oblong.png")
  png::writePNG(matrix(0.5, 16, 32), path)

  err <- expect_error(
    generateStimuli2IFC(base_face_files = list(base = path), n_trials = 2,
                        img_size = 32, stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE),
    "not square"
  )

  # A different branch from the img_size mismatch covered in test-fixed-bugs.R:
  # this one fires whatever img_size was asked for.
  msg <- conditionMessage(err)
  expect_match(msg, "16")
  expect_match(msg, "32")
})

test_that("generateStimuli2IFC rejects base_face_files that is not a list", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  png <- make_square_png(file.path(dir, "base.png"), size = 32)

  # This one calls stop() with no arguments after writing its explanation to
  # stderr(), so the condition carries an *empty* message and no regexp can
  # match it. Asserting only that it errors is the most this path allows;
  # giving it a real message belongs to item 14.
  expect_error(
    suppressWarnings(
      generateStimuli2IFC(base_face_files = png, n_trials = 2, img_size = 32,
                          stimulus_path = dir, ncores = 1,
                          save_as_png = FALSE, save_rdata = FALSE)
    )
  )
})

# --------------------------------------------------------------------------
# computeCumulativeCICorrelation(targetci = ...): never exercised
#
# The existing test covers only the no-targetci route, where the last
# correlation is 1 by construction -- true whatever the function does with the
# argument. Supplying a target CI is what the function is actually for.
# --------------------------------------------------------------------------

test_that("computeCumulativeCICorrelation correlates against a supplied target CI", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  stimuli <- 1:6
  responses <- rep(c(1, -1), 3)

  self <- suppressWarnings(computeCumulativeCICorrelation(
    stimuli = stimuli, responses = responses, baseimage = "base", rdata = rdata
  ))

  # The self-referential run ends at 1: the cumulative CI over every trial is
  # the final CI. That is the assertion the old test rested on.
  expect_equal(self[length(self)], 1)

  # A different target must not. Build one from the same stimulus set with the
  # responses reversed, so it is a real classification image rather than noise.
  target_ci <- suppressWarnings(generateCI(
    stimuli = stimuli, responses = rev(responses), baseimage = "base",
    rdata = rdata, save_as_png = FALSE
  ))

  against_target <- suppressWarnings(computeCumulativeCICorrelation(
    stimuli = stimuli, responses = responses, baseimage = "base", rdata = rdata,
    targetci = target_ci
  ))

  expect_length(against_target, length(self))
  expect_false(isTRUE(all.equal(against_target[length(against_target)], 1)))

  # And the value is the correlation with the supplied CI, not merely "not 1".
  full_ci <- suppressWarnings(generateCI(
    stimuli = stimuli, responses = responses, baseimage = "base",
    rdata = rdata, save_as_png = FALSE
  ))
  expect_equal(
    against_target[length(against_target)],
    cor(as.vector(full_ci$ci), as.vector(target_ci$ci))
  )
})
