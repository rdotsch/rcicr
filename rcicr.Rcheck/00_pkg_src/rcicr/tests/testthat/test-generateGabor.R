test_that("generateGabor returns a matrix of the requested size", {
  g <- generateGabor(16, 2, 45, 0, 25, 1)
  expect_equal(dim(g), c(16, 16))
})

test_that("generateGabor equals sinusoid times gaussian mask", {
  img_size <- 16
  cycles <- 2
  angle <- 45
  phase <- 0
  sigma <- 25
  contrast <- 1

  x0 <- scales::rescale(1:img_size, to = c(-.5, .5))
  gauss <- matlab::meshgrid(x0, x0)
  gauss_mask <- exp(-(((gauss$x^2) + (gauss$y^2)) / (2 * (sigma / img_size)^2)))

  expected <- generateSinusoid(img_size, cycles, angle, phase, contrast) * gauss_mask
  expect_equal(generateGabor(img_size, cycles, angle, phase, sigma, contrast), expected)
})

test_that("generateGabor with contrast 0 is all zero", {
  expect_equal(generateGabor(8, 2, 0, 0, 10, 0), matrix(0, 8, 8))
})
