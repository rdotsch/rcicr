# Every export in the latent module is gated, and nothing else is.

latent_exports <- c('latentGeneratorPCA')

test_that("the latent exports refuse to run with the option unset", {
  withr::local_options(rcicr.experimental = NULL)

  for (fn in latent_exports) {
    err <- expect_error(do.call(fn, list()))
    expect_match(conditionMessage(err), 'experimental', info = fn)
    expect_match(conditionMessage(err), 'options(rcicr.experimental = TRUE)',
                 fixed = TRUE, info = fn)
  }
})

test_that("the gate opens only for TRUE", {
  for (value in list(NULL, FALSE, 'TRUE', 1, NA)) {
    withr::local_options(rcicr.experimental = value)
    expect_error(requireExperimental('someFunction'), 'experimental')
  }

  withr::local_options(rcicr.experimental = TRUE)
  expect_true(requireExperimental('someFunction'))
})

test_that("the existing exports are not gated", {
  withr::local_options(rcicr.experimental = NULL)

  # generateNoisePattern() is representative: cheap, pure, and untouched by the
  # latent module. A gate leaking onto it would stop every stored script.
  expect_silent(p <- generateNoisePattern(img_size = 8, nscales = 1))
  expect_named(p, c('patches', 'patchIdx', 'noise_type', 'generator_version'))
})
