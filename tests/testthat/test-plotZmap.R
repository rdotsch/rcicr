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
