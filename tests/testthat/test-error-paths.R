# The suite was strong on the numeric paths and thin on the failure side:
# 9 assertions against 33 stop()/warning() calls. An unexercised guard is
# indistinguishable from one that works, which is how the `mask` argument, the
# plotZmap() mask and the single-trial 4096-parameter truncation each stayed
# broken for years -- and how a guard nothing ran came to paste 8,190
# characters of pixel values into its own error message.
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
  # -- unlike sigma, which is both, and which is why load() could overwrite the
  # argument with the file's field (the z-map sigma bug, fixed in #146).
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

test_that("applyMask accepts integer and logical matrix masks", {
  ci <- matrix(seq_len(64), 8, 8)

  integer_mask <- matrix(1L, 8, 8)
  integer_mask[1, 1] <- 0L
  integer_masked <- rcicr:::applyMask(ci, mask = integer_mask, img_size = 8)
  expect_true(is.na(integer_masked[1, 1]))
  expect_equal(integer_masked[!is.na(integer_masked)], ci[!is.na(integer_masked)])

  logical_mask <- matrix(TRUE, 8, 8)
  logical_mask[2, 2] <- FALSE
  logical_masked <- rcicr:::applyMask(ci, mask = logical_mask, img_size = 8)
  expect_true(is.na(logical_masked[2, 2]))
  expect_equal(logical_masked[!is.na(logical_masked)], ci[!is.na(logical_masked)])
})

# The region applyMask() actually masked, as a logical matrix. The accept tests
# below assert this rather than only expect_no_error(): a channel collapse that
# picked the wrong plane would still not error, and the release gate cannot
# catch it either -- its mask fixture is single-channel, so it never enters the
# collapse branch at all (tools/compare-release-output.R, make_mask()).
masked_region <- function(mask, n = 8) {
  is.na(rcicr:::applyMask(matrix(1, n, n), mask = mask, img_size = n))
}

# Asymmetric on both axes, so a transpose or a flip fails too.
mask_pattern <- function(n = 8) {
  p <- matrix(1, n, n)
  p[1:3, 1:5] <- 0
  p
}

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

  pattern <- mask_pattern()
  grey_as_rgb <- array(0, dim = c(8, 8, 3))
  for (i in 1:3) grey_as_rgb[, , i] <- pattern
  png::writePNG(grey_as_rgb, file.path(dir, "grey.png"))
  expect_equal(masked_region(file.path(dir, "grey.png")), pattern == 0)
})

test_that("applyMask ignores the alpha channel of a 2-channel PNG", {
  # png::readPNG() decodes an 8-bit greyscale-plus-alpha PNG to a 2-channel
  # array. This used to crash applyMask() unconditionally (subscript out of
  # bounds indexing channel 3), even though plotZmap()'s own inline mask code
  # already accepted the identical-planes case. Alpha is never colour
  # information and must be ignored regardless of whether it is uniform.
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()

  pattern <- mask_pattern()

  identical_alpha <- array(0, dim = c(8, 8, 2))
  identical_alpha[, , 1] <- pattern
  identical_alpha[, , 2] <- 1
  path_identical <- file.path(dir, "ga_identical.png")
  png::writePNG(identical_alpha, path_identical)
  expect_equal(masked_region(path_identical), pattern == 0)

  # Alpha is the inverse of the greyscale plane, so collapsing to the wrong
  # plane inverts the masked region rather than erroring -- silently, which is
  # the failure this assertion exists to catch.
  inverted_alpha <- array(0, dim = c(8, 8, 2))
  inverted_alpha[, , 1] <- pattern
  inverted_alpha[, , 2] <- 1 - pattern
  path_inverted <- file.path(dir, "ga_inverted.png")
  png::writePNG(inverted_alpha, path_inverted)
  expect_equal(masked_region(path_inverted), pattern == 0)

  # A non-binary alpha additionally proves alpha is dropped before the
  # 0/1 check, which it would otherwise fail.
  fractional_alpha <- array(0, dim = c(8, 8, 2))
  fractional_alpha[, , 1] <- pattern
  fractional_alpha[, , 2] <- 0.5
  path_fractional <- file.path(dir, "ga_fractional.png")
  png::writePNG(fractional_alpha, path_fractional)
  expect_equal(masked_region(path_fractional), pattern == 0)
})

test_that("applyMask ignores a differing alpha channel of a 4-channel PNG", {
  # An RGBA mask with a genuine RGB pattern and a constant (fully opaque)
  # alpha plane already succeeded through applyMask() before this change --
  # confirmed on unmodified main. This pins that it still does: alpha is
  # dropped before the colour-channel comparison, whether or not it matches
  # the other channels.
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()

  pattern <- mask_pattern()

  rgba <- array(0, dim = c(8, 8, 4))
  for (i in 1:3) rgba[, , i] <- pattern
  rgba[, , 4] <- 1
  path <- file.path(dir, "rgba.png")
  png::writePNG(rgba, path)
  expect_equal(masked_region(path), pattern == 0)

  # As above: an inverted alpha makes a wrong-plane collapse show up as an
  # inverted region instead of an error.
  rgba_inverted <- array(0, dim = c(8, 8, 4))
  for (i in 1:3) rgba_inverted[, , i] <- pattern
  rgba_inverted[, , 4] <- 1 - pattern
  path_inverted <- file.path(dir, "rgba_inverted.png")
  png::writePNG(rgba_inverted, path_inverted)
  expect_equal(masked_region(path_inverted), pattern == 0)
})

test_that("applyMask still rejects a 4-channel PNG whose RGB planes differ", {
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()

  rgba <- array(0, dim = c(8, 8, 4))
  rgba[, , 1] <- 1
  rgba[, , 4] <- 1
  path <- file.path(dir, "rgba_rgb_differs.png")
  png::writePNG(rgba, path)

  expect_error(
    rcicr:::applyMask(matrix(1, 8, 8), mask = path, img_size = 8),
    "not a greyscale image"
  )
})

test_that("applyMask rejects a multi-channel PNG mask with a singleton dimension", {
  # A 1-by-8 RGB mask against a 2-by-4 target: dropping the channel dimension
  # with a bare `[, , 1]` also drops the singleton row dimension, and the
  # resulting vector's NULL dim() makes the size check vacuously TRUE, so the
  # mismatched mask used to be applied by linear indexing instead of rejected.
  skip_if_not_installed("withr")
  dir <- withr::local_tempdir()

  rgb <- array(0, dim = c(1, 8, 3))
  path <- file.path(dir, "singleton.png")
  png::writePNG(rgb, path)

  expect_error(
    rcicr:::applyMask(matrix(1, 2, 4), mask = path, img_size = c(2, 4),
                      context = "z-map"),
    "same dimensions"
  )
})

test_that("applyMask accepts a rectangular img_size and reports both dims", {
  # plotZmap() is not restricted to square zmaps and passes
  # img_size = c(nrow(zmap), ncol(zmap)); applyMask()'s dim() == img_size
  # comparison already handles this via recycling, but the size-mismatch
  # message used to hardcode img_size for both dimensions.
  masked <- rcicr:::applyMask(matrix(1, 4, 6), mask = matrix(1, 4, 6),
                              img_size = c(4, 6))
  expect_equal(dim(masked), c(4, 6))

  err <- expect_error(
    rcicr:::applyMask(matrix(1, 4, 6), mask = matrix(0, 4, 4),
                      img_size = c(4, 6)),
    "same dimensions"
  )
  msg <- conditionMessage(err)
  expect_match(msg, "4 x 6")
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

  # This used to write its explanation to stderr() and then call stop() with no
  # arguments, so the condition carried an *empty* message: a caller could not
  # report why it failed, and this test could assert nothing beyond "it errors".
  # Match on the message, not just the failure -- an unpatterned expect_error()
  # passes when the function breaks for an entirely unrelated reason (#180).
  err <- expect_error(
    generateStimuli2IFC(base_face_files = png, n_trials = 2, img_size = 32,
                        stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE),
    "base_face_files must be a named list"
  )
  expect_match(conditionMessage(err), "character")
})

test_that("generateStimuli2IFC rejects an unnamed or empty base_face_files", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  png <- make_square_png(file.path(dir, "base.png"), size = 32)

  gen <- function(files) {
    generateStimuli2IFC(base_face_files = files, n_trials = 2, img_size = 32,
                        stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE)
  }

  # All three used to get past the type check and fail inside the parallel
  # workers with "attempt to select less than one element in get1index", naming
  # neither the argument nor the file -- issue #124's shape.
  expect_error(gen(list()), "base_face_files is empty")
  expect_error(gen(list(png)), "must be named")
  expect_error(gen(list(base = png, png)), "must be named")

  # A duplicate name silently dropped every entry after the first.
  expect_error(gen(list(base = png, base = png)), "duplicate names")
})

test_that("generateStimuli2IFC names the base image it cannot use", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  ok <- make_square_png(file.path(dir, "base.png"), size = 32)

  gen <- function(files) {
    generateStimuli2IFC(base_face_files = files, n_trials = 2, img_size = 32,
                        stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE)
  }

  # Every message identifies the offending entry by its name, so a run over
  # several base images says which one is wrong.
  expect_error(gen(list(good = ok, oblong = 42)),
               'Base image "oblong" must be a single file name')

  missing <- file.path(dir, "absent.png")
  err <- expect_error(gen(list(good = ok, gone = missing)),
                      'Base image "gone" does not exist')
  expect_match(conditionMessage(err), "absent\\.png")

  # Present, correctly named, and not actually a PNG: the reader's own
  # complaint is kept, with the file it came from attached.
  corrupt <- file.path(dir, "corrupt.png")
  writeLines("not a PNG", corrupt)
  err <- expect_error(gen(list(bad = corrupt)), 'Base image "bad"')
  expect_match(conditionMessage(err), "could not be read")
  expect_match(conditionMessage(err), "not in PNG format")
})

test_that("generateStimuli2IFC rejects a base image with no contrast", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()
  flat <- file.path(dir, "flat.png")
  png::writePNG(matrix(0.5, 32, 32), flat)

  # (img - min) / (max - min) is 0/0 on a uniform image. The all-NaN base face
  # used to be saved into the .Rdata and inherited by every CI computed from
  # the set, with no error and no warning (#176).
  err <- expect_error(
    generateStimuli2IFC(base_face_files = list(flat = flat), n_trials = 2,
                        img_size = 32, stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE),
    "no contrast"
  )
  expect_match(conditionMessage(err), "maximize_baseimage_contrast = FALSE")

  # The way out named in the message has to actually work: with the rescale
  # off, a flat base image is usable and its stimuli are finite.
  stims <- expect_no_error(
    generateStimuli2IFC(base_face_files = list(flat = flat), n_trials = 2,
                        img_size = 32, stimulus_path = dir, ncores = 1,
                        maximize_baseimage_contrast = FALSE,
                        return_as_dataframe = TRUE,
                        save_as_png = FALSE, save_rdata = FALSE)
  )
  expect_true(all(is.finite(as.matrix(stims))))
})

test_that("generateStimuli2IFC saves a base image the CI pipeline can use", {
  skip_if_not_installed("withr")

  # The NaN in #176 was invisible at generation time and surfaced only in the
  # saved .Rdata, so assert on what was actually written rather than on the
  # error alone.
  dir <- withr::local_tempdir()
  rdata <- make_fixture_rdata(dir)

  env <- new.env()
  load(rdata, envir = env)
  expect_true(all(is.finite(env$base_faces[["base"]])))
})

test_that("generateStimuli2IFC picks the reader from the extension, not the path", {
  skip_if_not_installed("withr")

  dir <- withr::local_tempdir()

  # A JPEG under a directory called "png". The old test was
  # grepl('png|PNG', filename) against the whole path, so this file was handed
  # to png::readPNG() and died with "file is not in PNG format".
  sub <- file.path(dir, "png")
  dir.create(sub)
  path <- file.path(sub, "face.jpg")
  # A ramp, not a flat field: a uniform base image is rejected for having no
  # contrast to maximize (#176), which would fail this test for the wrong reason.
  jpeg::writeJPEG(matrix(rep(seq(0, 1, length.out = 32), each = 32), 32, 32), path)

  expect_no_error(
    generateStimuli2IFC(base_face_files = list(base = path), n_trials = 2,
                        img_size = 32, stimulus_path = dir, ncores = 1,
                        save_as_png = FALSE, save_rdata = FALSE)
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
