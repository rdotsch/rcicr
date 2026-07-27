test_that("deg2rad converts degrees to radians", {
  expect_equal(deg2rad(0), 0)
  expect_equal(deg2rad(180), pi)
  expect_equal(deg2rad(90), pi / 2)
  expect_equal(deg2rad(360), 2 * pi)
  expect_equal(deg2rad(-90), -pi / 2)
})

test_that("deg2rad is vectorized", {
  expect_equal(deg2rad(c(0, 180, 360)), c(0, pi, 2 * pi))
})
