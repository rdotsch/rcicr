test_that("plotZmap writes a PNG file", {
  tmp <- withr::local_tempdir()
  zmap <- matrix(5, 8, 8)

  plotZmap(
    zmap = zmap, sigma = 3, threshold = 3, decoration = FALSE,
    targetpath = tmp, filename = "zmap", size = 200
  )

  expect_true(file.exists(file.path(tmp, "zmap.png")))
})

test_that("decoration = FALSE works with a background image", {
  # Regression test. This branch used `if (bgimage != '')`, and bgimage is
  # normally an image matrix, so the condition had length img_size^2 -- an error
  # on R >= 4.2 rather than a silent first-element match. Every generateCI() call
  # passes a background image, so the whole undecorated path was dead. Same root
  # cause as the `mask` bug (BACKLOG.md item 6).
  tmp <- withr::local_tempdir()

  plotZmap(
    zmap = matrix(5, 8, 8), bgimage = matrix(0.5, 8, 8), sigma = 3,
    threshold = 3, decoration = FALSE, targetpath = tmp,
    filename = "withbg", size = 200
  )

  expect_true(file.exists(file.path(tmp, "withbg.png")))
})

test_that("decoration = FALSE works on a device too small for default margins", {
  # Regression test. plot.new() was called before par(mar = c(0,0,0,0)), and
  # plot.new() rejects a device that cannot fit the current margins. Since
  # generateCI() sizes the device to img_size, small stimulus sets could not
  # produce a z-map at all ("figure margins too large").
  tmp <- withr::local_tempdir()

  plotZmap(
    zmap = matrix(5, 8, 8), sigma = 3, threshold = 3, decoration = FALSE,
    targetpath = tmp, filename = "small", size = 32
  )

  expect_true(file.exists(file.path(tmp, "small.png")))
})

# The three tests above assert only that a file appeared, which was the right
# bar for two hard-error regressions but says nothing about what was drawn: a
# plotZmap() that wrote a blank image would pass all of them. These read the
# PNG back and check the thresholding actually took effect. Rendering onto a
# uniform background makes "nothing was painted" exactly equal to "the output
# is the background", which is a comparison rather than an eyeball.
render_zmap <- function(dir, zmap, name, threshold = 3) {
  plotZmap(
    zmap = zmap, bgimage = matrix(0.5, 8, 8), sigma = 3, threshold = threshold,
    decoration = FALSE, targetpath = dir, filename = name, size = 64
  )
  png::readPNG(file.path(dir, paste0(name, ".png")))
}

test_that("sub-threshold z-scores are not painted and supra-threshold ones are", {
  tmp <- withr::local_tempdir()

  below <- render_zmap(tmp, matrix(1, 8, 8), "below")
  above <- render_zmap(tmp, matrix(5, 8, 8), "above")

  # Every cell is under threshold, so the z-map contributes nothing and the
  # result is the bare background: a single grey value across the image.
  expect_length(unique(as.vector(below)), 1)

  # ...and when the z-map does clear the threshold, the output is no longer the
  # background. Without this the assertion above is satisfied by a function
  # that never draws anything at all.
  expect_false(isTRUE(all.equal(above, below)))
})

test_that("the threshold is applied per pixel, not to the map as a whole", {
  # Half the z-map above threshold, half below. The painted half must match a
  # fully-painted render and the unpainted half must match a bare background,
  # which pins *where* the thresholding lands rather than just that some of it
  # happened.
  tmp <- withr::local_tempdir()

  split <- matrix(1, 8, 8)
  split[, 1:4] <- 5

  mixed <- render_zmap(tmp, split, "mixed")
  above <- render_zmap(tmp, matrix(5, 8, 8), "all_above")
  below <- render_zmap(tmp, matrix(1, 8, 8), "all_below")

  expect_equal(mixed[, 1:32, ], above[, 1:32, ])
  expect_equal(mixed[, 33:64, ], below[, 33:64, ])
})

test_that("z-scores beyond the negative threshold are painted too", {
  # The threshold is applied to abs(zmap), so a strongly negative region is
  # signal and must be drawn. A `zmap < threshold` implementation would drop
  # it, and nothing else in the suite would notice.
  tmp <- withr::local_tempdir()

  negative <- render_zmap(tmp, matrix(-5, 8, 8), "negative")
  below <- render_zmap(tmp, matrix(1, 8, 8), "neutral")

  expect_false(isTRUE(all.equal(negative, below)))
})

test_that("raising the threshold removes regions from the z-map", {
  # Monotonicity: the same z-map drawn at a stricter threshold must lose the
  # region between the two thresholds and fall back to background there.
  tmp <- withr::local_tempdir()
  zmap <- matrix(4, 8, 8)

  lenient <- render_zmap(tmp, zmap, "lenient", threshold = 3)
  strict <- render_zmap(tmp, zmap, "strict", threshold = 10)
  background <- render_zmap(tmp, matrix(0, 8, 8), "background", threshold = 3)

  expect_false(isTRUE(all.equal(lenient, background)))
  expect_equal(strict, background)
})

test_that("mismatched mask and zmap dimensions error", {
  tmp <- withr::local_tempdir()

  expect_error(
    plotZmap(
      zmap = matrix(0, 8, 8), mask = matrix(0, 4, 4), sigma = 3,
      targetpath = tmp, size = 200
    ),
    "not the same size"
  )
})
