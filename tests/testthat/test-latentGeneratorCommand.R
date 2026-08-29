# The external-program backend. Exercised end to end against a real subprocess,
# because the protocol between rcicr and the program is the whole of what this
# function is, and a mocked call would test none of it.

write_renderer <- function(dir, body) {
  script <- file.path(dir, 'renderer.R')
  writeLines(c(
    'args <- commandArgs(trailingOnly = TRUE)',
    'latents <- as.matrix(utils::read.csv(args[1], header = FALSE))',
    'outdir <- args[2]',
    body
  ), script)
  script
}

rscript <- function() file.path(R.home('bin'), 'Rscript')

# A renderer whose images depend on the latents, so a wrong row order or a
# dropped row shows up rather than passing as identical grey squares.
faithful_renderer <- c(
  'for (i in seq_len(nrow(latents))) {',
  '  value <- (latents[i, 1] + 1) / 2',
  '  png::writePNG(matrix(max(0, min(1, value)), 8, 8),',
  '                file.path(outdir, sprintf("%05d.png", i)))',
  '}'
)

test_that("a generator built on an external program satisfies the contract", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, faithful_renderer),
    latent_dim = 3, img_size = 8, latent_sd = 1
  )

  expect_s3_class(generator, 'rcicr_generator')
  expect_s3_class(generator, 'rcicr_generator_command')
  expect_identical(generator$latent_dim, 3L)
  expect_silent(validateGenerator(generator))
})

test_that("every latent is rendered, in the order given", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, faithful_renderer),
    latent_dim = 2, img_size = 8, latent_sd = 1
  )

  latents <- matrix(c(-1, 0, 0, 0, 1, 0), ncol = 2, byrow = TRUE)
  rendered <- renderLatent(generator, latents)

  expect_equal(dim(rendered), c(3, 8, 8))
  expect_equal(c(rendered[1, 1, 1], rendered[2, 1, 1], rendered[3, 1, 1]),
               c(0, 0.5, 1), tolerance = 1 / 255)
})

test_that("a colour image comes back as greyscale", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  colour <- c(
    'for (i in seq_len(nrow(latents))) {',
    '  rgb <- array(0, dim = c(8, 8, 3))',
    '  rgb[, , 1] <- 1; rgb[, , 2] <- 0.5; rgb[, , 3] <- 0',
    '  png::writePNG(rgb, file.path(outdir, sprintf("%05d.png", i)))',
    '}'
  )

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, colour),
    latent_dim = 2, img_size = 8, latent_sd = 1
  )

  rendered <- renderLatent(generator, c(0, 0))
  expect_equal(dim(rendered), c(1, 8, 8))
  expect_equal(rendered[1, 1, 1], 0.5, tolerance = 1 / 255)
})

test_that("a program that fails reports its own output", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  failing <- c('cat("the model file was not found\n"); quit(status = 3)')
  script <- write_renderer(dir, failing)

  # validateGenerator() renders once, so the failure surfaces at construction.
  err <- expect_error(latentGeneratorCommand(
    command = rscript(), args = script, latent_dim = 2, img_size = 8, latent_sd = 1
  ))
  expect_match(conditionMessage(err), 'exit status 3')
  expect_match(conditionMessage(err), 'the model file was not found')
})

test_that("a program that writes too few images says which one is missing", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  short <- c('png::writePNG(matrix(0.5, 8, 8), file.path(outdir, "00001.png"))')

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, short),
    latent_dim = 2, img_size = 8, latent_sd = 1
  )

  err <- expect_error(renderLatent(generator, matrix(0, nrow = 3, ncol = 2)))
  expect_match(conditionMessage(err), 'no image for latent 2')
})

test_that("a program writing the wrong size is refused rather than resized", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  wrong <- c(
    'for (i in seq_len(nrow(latents))) {',
    '  png::writePNG(matrix(0.5, 16, 16), file.path(outdir, sprintf("%05d.png", i)))',
    '}'
  )

  err <- expect_error(latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, wrong),
    latent_dim = 2, img_size = 8, latent_sd = 1
  ))
  expect_match(conditionMessage(err), '16 by 16')
  expect_match(conditionMessage(err), 'does not resize')
})

test_that("latent_sd given as one number is spread across the dimensions", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, faithful_renderer),
    latent_dim = 4, img_size = 8, latent_sd = 2.5
  )

  expect_equal(generator$latent_sd, rep(2.5, 4))
  expect_equal(generator$latent_mean, rep(0, 4))
})

test_that("the fingerprint follows the weights file when one is given", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  script <- write_renderer(dir, faithful_renderer)

  weights <- file.path(dir, 'model.bin')
  writeLines('version one', weights)

  build <- function() {
    latentGeneratorCommand(command = rscript(), args = script, latent_dim = 2,
                           img_size = 8, latent_sd = 1, weights = weights)
  }

  first <- build()
  expect_identical(first$fingerprint, build()$fingerprint)

  # Swapping the model behind the same script has to change the fingerprint, or
  # a classification image could be rendered through weights the participants
  # never saw.
  writeLines('version two', weights)
  expect_false(identical(first$fingerprint, build()$fingerprint))

  expect_error(latentGeneratorCommand(command = rscript(), args = script,
                                      latent_dim = 2, img_size = 8, latent_sd = 1,
                                      weights = file.path(dir, 'absent.bin')),
               'weights file does not exist')
})

test_that("a caller's own fingerprint is used as given", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, faithful_renderer),
    latent_dim = 2, img_size = 8, latent_sd = 1, fingerprint = 'ffhq-1024-mine'
  )

  expect_identical(generator$fingerprint, 'ffhq-1024-mine')
})

test_that("the whole pipeline runs on an external generator", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  out <- withr::local_tempdir()

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, faithful_renderer),
    latent_dim = 3, img_size = 8, latent_sd = 1
  )

  stimuli <- generateStimuliLatent2IFC(generator, n_trials = 6, stimulus_path = out,
                                       seed = 1, save_as_png = FALSE)

  # This generator cannot travel inside the file, so the analysis half has to be
  # handed it back. That refusal is the point of the fingerprint.
  expect_error(generateLatentCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
                                rdata = stimuli$rdata, save_as_png = FALSE),
               'too ')

  ci <- generateLatentCI(stimuli = 1:6, responses = rep(c(1, -1), 3),
                         rdata = stimuli$rdata, generator = generator,
                         save_as_png = FALSE)
  expect_equal(dim(ci$ci_image), c(8, 8))
  expect_length(ci$direction, 3)
})

test_that("the shipped StyleGAN helper is installed and documents its protocol", {
  script <- system.file('python', 'rcicr_stylegan.py', package = 'rcicr')

  expect_true(nzchar(script))
  expect_true(file.exists(script))

  # The script and latentGeneratorCommand() have to agree on the protocol, and
  # nothing else checks that they do.
  source_text <- paste(readLines(script, warn = FALSE), collapse = '\n')
  expect_match(source_text, '%05d.png', fixed = TRUE)
  expect_match(source_text, 'latents.csv', fixed = TRUE)
})

test_that("a grey-plus-alpha image keeps its grey and drops its alpha", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  # Two channels is grey plus alpha, not two colours. Averaging them folds
  # opacity into the pixel value, so a fully opaque black pixel comes back as
  # 0.5 rather than 0 and every stimulus this renderer produces is wrong.
  grey_alpha <- c(
    'for (i in seq_len(nrow(latents))) {',
    '  ga <- array(0, dim = c(8, 8, 2))',
    '  ga[, , 1] <- max(0, min(1, (latents[i, 1] + 1) / 2))',
    '  ga[, , 2] <- 1',
    '  png::writePNG(ga, file.path(outdir, sprintf("%05d.png", i)))',
    '}'
  )

  generator <- latentGeneratorCommand(
    command = rscript(), args = write_renderer(dir, grey_alpha),
    latent_dim = 2, img_size = 8, latent_sd = 1
  )

  latents <- matrix(c(-1, 0, 1, 0), ncol = 2, byrow = TRUE)
  rendered <- renderLatent(generator, latents)

  expect_equal(dim(rendered), c(2, 8, 8))
  expect_equal(c(rendered[1, 1, 1], rendered[2, 1, 1]), c(0, 1), tolerance = 1 / 255)

  # An opaque grey image is still fully opaque, so alpha carries no information
  # and averaging it in would pull both values towards 0.5.
  expect_false(isTRUE(all.equal(rendered[1, 1, 1], 0.5, tolerance = 1 / 255)))
})
