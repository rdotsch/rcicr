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
