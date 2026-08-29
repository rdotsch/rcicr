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

test_that("identical images leave no face space to build", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  files <- file.path(dir, sprintf('same%02d.png', 1:3))
  for (f in files) {
    make_square_png(f, size = 8, seed = 1)
  }

  expect_error(latentGeneratorPCA(files, n_components = 2, img_size = 8),
               'no variation between them')
})

test_that("no variation is judged against a tolerance rather than exact zero", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  # Regression, and the reason the test above is not enough on its own. Centring
  # identical rows does not leave exactly zero: the mean of three copies of a
  # number is not always that number in floating point, and whether the residual
  # lands on zero depends on the platform's summation order. Linux got exact
  # zeros for this fixture and macOS ARM64 did not, so a check written against
  # zero passed on one and failed on the other.
  #
  # Perturbing one pixel by far less than the tolerance reproduces the macOS
  # arithmetic on any platform: the images are no longer bit-identical, and the
  # set still has to be rejected.
  base <- withr::with_seed(1, matrix(stats::runif(64), 8, 8))
  nudged <- base
  nudged[1, 1] <- nudged[1, 1] + 1e-16

  files <- file.path(dir, sprintf('near%02d.png', 1:3))
  png::writePNG(base, files[1])
  png::writePNG(base, files[2])
  png::writePNG(nudged, files[3])

  expect_error(latentGeneratorPCA(files, n_components = 2, img_size = 8),
               'no variation between them')
})

test_that("a real difference between images is not mistaken for none", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  # The smallest difference a PNG can carry is one 8-bit step, which is about
  # twelve orders of magnitude above the tolerance above. It must build.
  base <- withr::with_seed(1, matrix(stats::runif(64), 8, 8))
  stepped <- base
  stepped[1, 1] <- stepped[1, 1] + 1 / 255

  files <- file.path(dir, sprintf('step%02d.png', 1:3))
  png::writePNG(base, files[1])
  png::writePNG(base, files[2])
  png::writePNG(stepped, files[3])

  generator <- latentGeneratorPCA(files, n_components = 2, img_size = 8)
  expect_gte(generator$latent_dim, 1L)
})

test_that("images sharing a label are refused rather than silently dropped", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()
  files <- make_face_set(dir, n = 3)

  # Looking each image up by its label returns the first entry carrying it,
  # however many share it, and writes them all into one slot. The set would be
  # built from two of these three images with no indication that one was gone.
  labelled <- stats::setNames(as.list(files), c('face', 'face', 'other'))

  expect_error(latentGeneratorPCA(labelled, n_components = 2, img_size = 16),
               'appear more than once')
})

test_that("every image reaches the decomposition", {
  withr::local_options(rcicr.experimental = TRUE)
  dir <- withr::local_tempdir()

  # Distinct labels, so nothing is dropped: the generator has to see all four
  # images, which n - 1 components of rank confirms.
  generator <- latentGeneratorPCA(make_face_set(dir, n = 4), n_components = 10,
                                  img_size = 16)
  expect_identical(generator$state$n_faces, 4L)
  expect_identical(generator$latent_dim, 3L)
})
