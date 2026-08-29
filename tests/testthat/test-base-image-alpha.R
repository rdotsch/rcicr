# The greyscale conversion a base image goes through on its way into the
# .Rdata. It decides the pixels of every stimulus a participant sees, so the
# assertions here read the base face back out of the file generateStimuli2IFC()
# wrote rather than testing the conversion in isolation.

write_rgba <- function(file, face, alpha) {
  img <- array(0, dim = c(nrow(face), ncol(face), 4))
  for (k in 1:3) {
    img[, , k] <- face
  }
  img[, , 4] <- alpha
  png::writePNG(img, file)
  invisible(file)
}

write_grey_alpha <- function(file, face, alpha) {
  img <- array(0, dim = c(nrow(face), ncol(face), 2))
  img[, , 1] <- face
  img[, , 2] <- alpha
  png::writePNG(img, file)
  invisible(file)
}

# The base face as it was stored, which is what the stimuli were composed from.
stored_base_face <- function(file, out, maximize) {
  generateStimuli2IFC(list(face = file), n_trials = 2, img_size = 16,
    stimulus_path = out, nscales = 2, ncores = 1,
    maximize_baseimage_contrast = maximize, save_as_png = FALSE
  )

  rdata <- list.files(out, pattern = '[.]Rdata$', full.names = TRUE)
  loaded <- new.env()
  load(rdata[1], envir = loaded)

  loaded$base_faces[['face']]
}

rescale01 <- function(x) (x - min(x)) / diff(range(x))

test_that("alpha is dropped rather than averaged into the base face", {
  dir <- withr::local_tempdir()
  face <- withr::with_seed(1, matrix(stats::runif(256), 16, 16))

  rgba <- write_rgba(file.path(dir, 'rgba.png'), face, 1)
  grey_alpha <- write_grey_alpha(file.path(dir, 'ga.png'), face, 1)

  # Averaging all four channels of an opaque RGBA image gives (r + g + b + 1) / 4,
  # so a black pixel was stored as 0.25 and a grey-plus-alpha one as 0.5. Only
  # the colour channels carry the face.
  expect_equal(stored_base_face(rgba, withr::local_tempdir(), FALSE),
               png::readPNG(rgba)[, , 1])
  expect_equal(stored_base_face(grey_alpha, withr::local_tempdir(), FALSE),
               png::readPNG(grey_alpha)[, , 1])

  # And not what the old conversion produced, so this cannot pass vacuously.
  expect_false(isTRUE(all.equal(stored_base_face(rgba, withr::local_tempdir(), FALSE),
                                apply(png::readPNG(rgba), c(1, 2), mean))))
})

test_that("RGB and greyscale base faces are read exactly as before", {
  dir <- withr::local_tempdir()
  face <- withr::with_seed(2, matrix(stats::runif(256), 16, 16))

  # Channels that really are colours are still averaged, and a two-dimensional
  # image is untouched. This is what keeps the change invisible to every base
  # face without an alpha channel.
  rgb <- file.path(dir, 'rgb.png')
  img <- array(0, dim = c(16, 16, 3))
  for (k in 1:3) {
    img[, , k] <- face
  }
  png::writePNG(img, rgb)

  plain <- file.path(dir, 'grey.png')
  png::writePNG(face, plain)

  expect_equal(stored_base_face(rgb, withr::local_tempdir(), FALSE),
               apply(png::readPNG(rgb), c(1, 2), mean))
  expect_equal(stored_base_face(plain, withr::local_tempdir(), FALSE),
               png::readPNG(plain))
})

test_that("a constant alpha changes nothing once contrast is maximized", {
  dir <- withr::local_tempdir()
  face <- withr::with_seed(3, matrix(stats::runif(256), 16, 16))
  rgba <- write_rgba(file.path(dir, 'opaque.png'), face, 1)

  # Adding a constant and then rescaling to [0, 1] is an affine transform of the
  # same image, so a fully opaque base face comes out identical either way. That
  # is the whole of the compatibility argument, and it is asserted rather than
  # reasoned about: the old conversion's answer is required to match the new one.
  stored <- png::readPNG(rgba)
  got <- stored_base_face(rgba, withr::local_tempdir(), TRUE)

  expect_equal(got, rescale01(stored[, , 1]))
  expect_equal(got, rescale01(apply(stored, c(1, 2), mean)))
})

test_that("a varying alpha no longer bends the base face", {
  dir <- withr::local_tempdir()
  ax <- seq(-1, 1, length.out = 16)
  face <- outer(rev(ax), ax, function(y, x) exp(-((x^2 + y^2) / 0.5)))
  cutout <- outer(rev(ax), ax, function(y, x) as.numeric(x^2 + y^2 < 0.6))

  # A face cut out on a transparent background. Alpha is not constant here, so
  # the affine argument above does not apply and the distortion survived the
  # rescaling. This is the case that reaches a real study.
  file <- write_rgba(file.path(dir, 'cutout.png'), face, cutout)
  stored <- png::readPNG(file)

  got <- stored_base_face(file, withr::local_tempdir(), TRUE)
  expect_equal(got, rescale01(stored[, , 1]))

  folded <- rescale01(apply(stored, c(1, 2), mean))
  expect_false(isTRUE(all.equal(got, folded)))
  expect_gt(max(abs(got - folded)), 0.1)
})
