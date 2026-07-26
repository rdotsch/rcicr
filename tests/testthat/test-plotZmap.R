test_that("plotZmap writes a PNG file", {
  tmp <- withr::local_tempdir()
  zmap <- matrix(5, 8, 8)

  plotZmap(
    zmap = zmap, sigma = 3, threshold = 3, decoration = FALSE,
    targetpath = tmp, filename = "zmap", size = 200
  )

  expect_true(file.exists(file.path(tmp, "zmap.png")))
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
