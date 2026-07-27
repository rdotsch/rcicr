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

test_that("autoscale refreshes $combined so it agrees with the new $scaled", {
  # $combined used to be left at whatever the caller passed in, so after
  # autoscale() the two fields described different images. batchGenerateCI()
  # scales with 'none' before handing over, so its $combined was an overlay of
  # the *unscaled* noise -- and it disagreed with the PNG the same call wrote.
  cis <- list(
    a = list(
      ci       = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2),
      base     = matrix(0.5, 2, 2),
      scaled   = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2), # 'none' scaling
      combined = (matrix(c(-0.2, 0.1, 0, 0.05), 2, 2) + matrix(0.5, 2, 2)) / 2
    )
  )

  result <- autoscale(cis, save_as_pngs = FALSE)

  expect_equal(result$a$combined, (result$a$scaled + result$a$base) / 2)
  expect_false(isTRUE(all.equal(result$a$combined, cis$a$combined)))
})

test_that("autoscale still works when there is no base image to combine with", {
  # The documented example passes only $ci and $base; a caller may reasonably
  # pass neither $combined nor $base, and that must not error.
  cis <- list(a = list(ci = matrix(c(-0.2, 0.1, 0, 0.05), 2, 2)))

  result <- autoscale(cis, save_as_pngs = FALSE)

  expect_equal(result$a$scaled, (cis$a$ci + 0.2) / 0.4)
  expect_null(result$a$combined)
})
