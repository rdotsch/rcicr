# The noise basis without matlab or scales (#208) is the one they built.
# Reference values were captured from the matlab-based code (rcicr 1.5.0.9000
# at bdafa78) into fixtures/noise-pattern-matlab-reference.rds. The gate's
# smallest patch is 4 pixels, so these pin what it cannot: a 1-pixel patch,
# where img_size == 2^(nscales - 1).

ref <- readRDS(test_path("fixtures", "noise-pattern-matlab-reference.rds"))

test_that("a 1-pixel sinusoid and Gabor patch keep their matlab-era values", {
  sinusoid <- mapply(function(a, ph) generateSinusoid(1, 2, a, ph, 1), ref$grid$angle, ref$grid$phase)
  gabor <- mapply(function(a, ph) generateGabor(1, 1.5, a, ph, 25, 1), ref$grid$angle, ref$grid$phase)
  expect_equal(sinusoid, ref$sinusoid_size1, tolerance = 1e-12)
  expect_equal(gabor, ref$gabor_size1, tolerance = 1e-12)
  # A single pixel holds `cycles`, not 0; 1.25 cycles tells the two apart.
  expect_equal(generateSinusoid(1, 1.25, 0, 0, 1)[1, 1], sin(2 * pi * 1.25))
})

test_that("a 16px, 5-scale noise pattern, finest patch 1 pixel, is unchanged", {
  for (case in list(list(ref$pattern_sinusoid, generateNoisePattern(16, nscales = 5)),
                    list(ref$pattern_gabor, generateNoisePattern(16, nscales = 5, noise_type = "gabor",
                                                                 sigma = 25)))) {
    expect_identical(case[[2]]$patchIdx, case[[1]]$patchIdx)
    expect_equal(case[[2]]$patches, case[[1]]$patches, tolerance = 1e-12)
  }
})
