test_that("generateSinusoid returns a matrix of the requested size", {
  s <- generateSinusoid(16, 2, 90, 0, 1.0)
  expect_equal(dim(s), c(16, 16))
})

test_that("generateSinusoid values are bounded by contrast", {
  s1 <- generateSinusoid(16, 2, 45, 0, 1.0)
  expect_true(all(abs(s1) <= 1 + 1e-8))

  s2 <- generateSinusoid(16, 2, 45, 0, 0.5)
  expect_true(all(abs(s2) <= 0.5 + 1e-8))
})

test_that("generateSinusoid with contrast 0 is all zero", {
  expect_equal(generateSinusoid(8, 2, 45, 0, 0), matrix(0, 8, 8))
})

test_that("angle 0 and angle 90 are transposes of each other", {
  # At angle=0, cos(angle)=1/sin(angle)=0, so only `sinepatch` (not transposed) is used.
  # At angle=90, cos(angle)=0/sin(angle)=1, so only t(sinepatch) is used.
  s0 <- generateSinusoid(16, 2, 0, 0, 1)
  s90 <- generateSinusoid(16, 2, 90, 0, 1)
  expect_equal(s0, t(s90), tolerance = 1e-10)
})

test_that("a phase shift of pi negates the pattern", {
  s <- generateSinusoid(16, 2, 30, 0, 1)
  s_shifted <- generateSinusoid(16, 2, 30, pi, 1)
  expect_equal(s_shifted, -s, tolerance = 1e-10)
})
