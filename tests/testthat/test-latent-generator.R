# The generator contract, checked against deliberately broken backends. A
# generator that fails the contract must fail on the call rather than partway
# through a render loop.

make_generator <- function(..., probe = TRUE) {
  defaults <- list(
    kind = 'stub',
    latent_dim = 2,
    img_size = 4,
    space = 'z',
    latent_mean = c(0, 0),
    latent_sd = c(1, 1),
    render = function(latents) array(0.5, dim = c(nrow(latents), 4, 4)),
    fingerprint = 'stub:1'
  )
  do.call(rcicrGenerator, utils::modifyList(defaults, list(...)))
}

test_that("a well formed generator validates and renders", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- make_generator()

  expect_s3_class(generator, 'rcicr_generator')
  expect_s3_class(generator, 'rcicr_generator_stub')
  expect_identical(validateGenerator(generator), generator)

  rendered <- renderLatent(generator, matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE))
  expect_equal(dim(rendered), c(2, 4, 4))
})

test_that("renderLatent accepts a bare vector as one latent", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- make_generator()
  expect_equal(dim(renderLatent(generator, c(0, 0))), c(1, 4, 4))
})

test_that("a non-generator is rejected by class", {
  err <- expect_error(validateGenerator(list(latent_dim = 2)))
  expect_match(conditionMessage(err), 'latentGenerator')
})

test_that("a missing field is named", {
  generator <- make_generator()
  generator$fingerprint <- NULL

  expect_error(validateGenerator(generator), 'fingerprint')
})

test_that("latent_mean and latent_sd must match latent_dim", {
  expect_error(validateGenerator(make_generator(latent_mean = c(0, 0, 0))),
               'latent_mean has length 3')
  expect_error(validateGenerator(make_generator(latent_sd = 1)),
               'latent_sd has length 1')
})

test_that("a negative or non-finite scale is rejected", {
  expect_error(validateGenerator(make_generator(latent_sd = c(1, -1))), 'negative')
  expect_error(validateGenerator(make_generator(latent_mean = c(0, NA))), 'non-finite')
})

test_that("a render() returning the wrong shape is rejected on the call", {
  wrong_size <- make_generator(render = function(latents) array(0.5, dim = c(nrow(latents), 8, 8)))
  expect_error(validateGenerator(wrong_size), '1 by 4 by 4')

  not_an_array <- make_generator(render = function(latents) 0.5)
  expect_error(validateGenerator(not_an_array), 'render\\(\\) must return an array')
})

test_that("a render() leaving [0, 1] is rejected, naming the range", {
  out_of_range <- make_generator(render = function(latents) array(1.5, dim = c(nrow(latents), 4, 4)))
  err <- expect_error(validateGenerator(out_of_range))
  expect_match(conditionMessage(err), 'outside')
  expect_match(conditionMessage(err), '1.5')

  non_finite <- make_generator(render = function(latents) array(NaN, dim = c(nrow(latents), 4, 4)))
  expect_error(validateGenerator(non_finite), 'non-finite')
})

test_that("latents of the wrong width are rejected before reaching render()", {
  withr::local_options(rcicr.experimental = TRUE)
  generator <- make_generator(render = function(latents) stop("render should not be reached"))
  expect_error(renderLatent(generator, matrix(0, nrow = 1, ncol = 5), validate = FALSE),
               'latent_dim is 2')
})

test_that("fingerprints separate different content and repeat for the same", {
  expect_identical(fingerprintNumeric(1:10), fingerprintNumeric(1:10))
  expect_false(identical(fingerprintNumeric(1:10), fingerprintNumeric(2:11)))

  # Ordering has to move the fingerprint, or a permuted component matrix would
  # pass as the generator that made the stimuli.
  expect_false(identical(fingerprintNumeric(c(1, 2, 3)), fingerprintNumeric(c(3, 2, 1))))
})

test_that("validating a generator does not need the experimental option", {
  # validateGenerator() renders the centre of the space to check the contract.
  # That is the package checking itself, so it must work whether or not the
  # caller has opened the gate, or building a generator would fail inside its
  # own validation.
  withr::local_options(rcicr.experimental = NULL)

  expect_silent(validateGenerator(make_generator()))
  expect_error(renderLatent(make_generator(), c(0, 0)), 'experimental')
})
