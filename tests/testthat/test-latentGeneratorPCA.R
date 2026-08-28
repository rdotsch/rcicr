# The eigenface backend. It is the reason the rest of the latent module can be
# tested at all, so its own reconstruction properties are checked rather than
# assumed.

test_that("a generator is built and satisfies the contract", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)

  expect_s3_class(generator, 'rcicr_generator')
  expect_identical(generator$kind, 'pca')
  expect_identical(generator$latent_dim, 3L)
  expect_identical(generator$img_size, 16L)
  expect_silent(validateGenerator(generator))
})

test_that("img_size is taken from the images when not given", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir, size = 8), n_components = 2)
  expect_identical(generator$img_size, 8L)
})

test_that("the origin of the space renders the mean face", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir)

  generator <- latentGeneratorPCA(files, n_components = 3, img_size = 16)

  faces <- lapply(files, function(f) png::readPNG(f))
  faces <- lapply(faces, function(img) (img - min(img)) / (max(img) - min(img)))
  mean_face <- Reduce(`+`, faces) / length(faces)

  rendered <- renderLatent(generator, generator$latent_mean)[1, , ]
  expect_equal(rendered, mean_face, tolerance = 1e-8)
})

test_that("keeping every available component reconstructs the training faces", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir, n = 5)

  # n - 1 components span the centred data exactly, so projecting a training
  # face onto them and rendering it back has to return that face. If this fails
  # the components or the vectorisation order are wrong.
  generator <- latentGeneratorPCA(files, n_components = 4, img_size = 16)

  original <- png::readPNG(files[3])
  original <- (original - min(original)) / (max(original) - min(original))

  scores <- as.vector(as.vector(original) - generator$state$mean_face) %*% generator$state$components
  expect_equal(renderLatent(generator, scores)[1, , ], original, tolerance = 1e-8)
})

test_that("latent_sd reports the spread of the training faces along each component", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir, n = 8)

  generator <- latentGeneratorPCA(files, n_components = 4, img_size = 16)

  faces <- t(vapply(files, function(f) {
    img <- png::readPNG(f)
    as.vector((img - min(img)) / (max(img) - min(img)))
  }, numeric(16 * 16)))
  scores <- sweep(faces, 2, generator$state$mean_face, '-') %*% generator$state$components

  expect_equal(generator$latent_sd, apply(scores, 2, stats::sd), tolerance = 1e-8)
})

test_that("components are capped at the rank of the centred data", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir, n = 4), n_components = 50, img_size = 16)
  expect_identical(generator$latent_dim, 3L)
})

test_that("rendered pixels stay inside [0, 1] however far the latent goes", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  generator <- latentGeneratorPCA(make_face_set(dir), n_components = 3, img_size = 16)
  far <- matrix(c(50, -50, 50), nrow = 1)

  rendered <- renderLatent(generator, far)
  expect_gte(min(rendered), 0)
  expect_lte(max(rendered), 1)
})

test_that("the fingerprint depends on the images and not on the run", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir)

  first <- latentGeneratorPCA(files, n_components = 3, img_size = 16)
  again <- latentGeneratorPCA(files, n_components = 3, img_size = 16)
  expect_identical(first$fingerprint, again$fingerprint)

  other_dir <- withr::local_tempdir()
  other <- latentGeneratorPCA(make_face_set(other_dir, n = 7), n_components = 3, img_size = 16)
  expect_false(identical(first$fingerprint, other$fingerprint))
})

test_that("a named list of files is accepted, as generateStimuli2IFC takes", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir, n = 4)

  as_list <- latentGeneratorPCA(as.list(stats::setNames(files, paste0('f', 1:4))),
                                n_components = 2, img_size = 16)
  as_vector <- latentGeneratorPCA(files, n_components = 2, img_size = 16)
  expect_identical(as_list$fingerprint, as_vector$fingerprint)
})

test_that("bad input is rejected with a message naming the problem", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir)

  expect_error(latentGeneratorPCA(files[1], img_size = 16), 'at least 2 images')
  expect_error(latentGeneratorPCA(c(files[1], file.path(dir, 'nope.png')), img_size = 16),
               'do not exist')
  expect_error(latentGeneratorPCA(c(files[1], file.path(dir, 'face.bmp')), img_size = 16),
               'do not exist')
  expect_error(latentGeneratorPCA(1:4, img_size = 16), 'character vector')

  # The size check applies to the first image as well as the rest.
  odd_dir <- withr::local_tempdir()
  odd <- make_face_set(odd_dir, n = 3, size = 8)
  expect_error(latentGeneratorPCA(odd, img_size = 16), 'does not resize')
})

test_that("identical images leave no face space to build, on every platform", {
  # Regression: this was decided from the singular values, whose threshold is
  # relative to the largest of them and so becomes zero when every one is zero.
  # macOS ARM64's LAPACK returns values around 1e-18 where Linux returns exact
  # zeros, so the call errored on one platform and returned a generator with a
  # meaningless component on the other.
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  files <- file.path(dir, sprintf('same%02d.png', 1:3))
  for (f in files) {
    make_square_png(f, size = 8, seed = 1)
  }

  expect_error(latentGeneratorPCA(files, n_components = 2, img_size = 8),
               'no variation between them')
})
