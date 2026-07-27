test_that("autoscale picks the constant that covers the widest range across all CIs and scales accordingly", {
  cis <- list(
    a = list(ci = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2), base = matrix(0.5, 2, 2)),
    b = list(ci = matrix(c(-0.05, 0.3, -0.1, 0.2), 2, 2), base = matrix(0.5, 2, 2))
  )

  result <- autoscale(cis, save_as_pngs = FALSE)

  # range across all cis: min = -0.2, max = 0.3 -> constant = max(abs(min), max) = 0.3
  constant <- 0.3

  expect_equal(result$a$scaled, (cis$a$ci + constant) / (2 * constant))
  expect_equal(result$b$scaled, (cis$b$ci + constant) / (2 * constant))
})

test_that("autoscale leaves $combined untouched and only rewrites $scaled", {
  # This is intended behaviour, not an oversight, and the test exists to keep it
  # that way. autoscale() rescales the noise; it does not re-derive the overlay.
  # Whatever combination the caller made before this call survives it, so an
  # existing analysis script that plots $combined keeps plotting the same image.
  # $scaled is the field that carries the autoscaled result.
  #
  # If this test fails because someone made autoscale() rewrite $combined, that
  # is a silent change to what published scripts render -- fix the code, not the
  # assertion.
  original_combined <- (matrix(c(-0.2, 0.1, 0, 0.05), 2, 2) + matrix(0.5, 2, 2)) / 2
  cis <- list(
    a = list(
      ci       = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2),
      base     = matrix(0.5, 2, 2),
      scaled   = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2), # as left by 'none' scaling
      combined = original_combined
    )
  )

  result <- autoscale(cis, save_as_pngs = FALSE)

  expect_identical(result$a$combined, original_combined)
  expect_false(isTRUE(all.equal(result$a$scaled, cis$a$scaled)))
})

test_that("autoscale writes a PNG built from the autoscaled noise, not from $combined", {
  # The file on disk is the autoscaled overlay even though the in-memory
  # $combined is not: (scaled + base) / 2 is recomputed at save time. Anyone
  # wanting that image in R should build it the same way.
  tmp <- withr::local_tempdir()
  cis <- list(
    a = list(
      ci       = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2),
      base     = matrix(0.5, 2, 2),
      combined = matrix(0, 2, 2) # deliberately wrong, must not be used
    )
  )

  result <- autoscale(cis, save_as_pngs = TRUE, targetpath = tmp)

  written <- png::readPNG(file.path(tmp, "a_autoscaled.png"))
  expect_equal(written, (result$a$scaled + result$a$base) / 2, tolerance = 1 / 255)
})
