# The base-image reader shared by both pipelines. Its greyscale conversion is
# the part worth pinning: it decides the pixels of every stimulus a participant
# sees, and it is reached identically from generateStimuli2IFC() and from
# latentGeneratorPCA().

square_rgba <- function(file, face, alpha) {
  img <- array(0, dim = c(nrow(face), ncol(face), 4))
  for (k in 1:3) {
    img[, , k] <- face
  }
  img[, , 4] <- alpha
  png::writePNG(img, file)
}

square_grey_alpha <- function(file, face, alpha) {
  img <- array(0, dim = c(nrow(face), ncol(face), 2))
  img[, , 1] <- face
  img[, , 2] <- alpha
  png::writePNG(img, file)
}

test_that("alpha is dropped rather than averaged into the greyscale value", {
  dir <- withr::local_tempdir()
  face <- withr::with_seed(1, matrix(stats::runif(64), 8, 8))

  rgba <- file.path(dir, 'rgba.png')
  grey_alpha <- file.path(dir, 'ga.png')
  square_rgba(rgba, face, 1)
  square_grey_alpha(grey_alpha, face, 1)

  # Averaging all four channels of an opaque RGBA image gives (r+g+b+1)/4, so a
  # black pixel read back as 0.25 and a grey-plus-alpha one as 0.5. Only the
  # colour channels carry the face.
  expect_equal(readBaseImage(rgba, 'rgba', NULL, FALSE, 'test'),
               png::readPNG(rgba)[, , 1])
  expect_equal(readBaseImage(grey_alpha, 'ga', NULL, FALSE, 'test'),
               png::readPNG(grey_alpha)[, , 1])

  # Not the value the old conversion produced: (r + g + b + 1) / 4 for RGBA and
  # (grey + 1) / 2 for grey plus alpha, both of which pull the image towards
  # white by a constant.
  expect_false(isTRUE(all.equal(readBaseImage(rgba, 'rgba', NULL, FALSE, 'test'),
                                apply(png::readPNG(rgba), c(1, 2), mean))))
})

test_that("RGB and greyscale images are read exactly as before", {
  dir <- withr::local_tempdir()
  face <- withr::with_seed(2, matrix(stats::runif(64), 8, 8))

  # The channels that are colours are still averaged, and a two-dimensional
  # image is untouched. This is what keeps the change invisible to every base
  # face that has no alpha, which is the reproducibility claim in NEWS.md.
  rgb <- file.path(dir, 'rgb.png')
  img <- array(0, dim = c(8, 8, 3))
  img[, , 1] <- face
  img[, , 2] <- face
  img[, , 3] <- face
  png::writePNG(img, rgb)

  plain <- file.path(dir, 'grey.png')
  png::writePNG(face, plain)

  expect_equal(readBaseImage(rgb, 'rgb', NULL, FALSE, 'test'),
               apply(png::readPNG(rgb), c(1, 2), mean))
  expect_equal(readBaseImage(plain, 'grey', NULL, FALSE, 'test'), png::readPNG(plain))
})

test_that("a constant alpha changes nothing once contrast is maximized", {
  dir <- withr::local_tempdir()
  face <- withr::with_seed(3, matrix(stats::runif(64), 8, 8))

  # Adding a constant and then rescaling to [0, 1] is an affine transform of the
  # same image, so a fully opaque base face comes out identical either way. This
  # is why the default path is unaffected, and it is worth asserting rather than
  # reasoning about, because it is the whole of the compatibility argument.
  rgba <- file.path(dir, 'opaque.png')
  square_rgba(rgba, face, 1)

  rescale <- function(x) (x - min(x)) / diff(range(x))
  stored <- png::readPNG(rgba)

  rescaled <- readBaseImage(rgba, 'opaque', NULL, TRUE, 'test')
  expect_equal(rescaled, rescale(stored[, , 1]))
  expect_equal(rescaled, rescale(apply(stored, c(1, 2), mean)))
})

test_that("a varying alpha no longer bends the base face", {
  dir <- withr::local_tempdir()
  ax <- seq(-1, 1, length.out = 16)
  face <- outer(rev(ax), ax, function(y, x) exp(-((x^2 + y^2) / 0.5)))
  cutout <- outer(rev(ax), ax, function(y, x) as.numeric(x^2 + y^2 < 0.6))

  # A face cut out on a transparent background: alpha is not constant, so the
  # affine argument above does not apply and the distortion survives the
  # rescaling. This is the case that reaches a real study.
  file <- file.path(dir, 'cutout.png')
  square_rgba(file, face, cutout)

  # Compared against the colour channel the PNG actually stores, so 8-bit
  # quantization is on both sides and the only thing under test is whether alpha
  # was folded in.
  stored <- png::readPNG(file)
  rescale <- function(x) (x - min(x)) / diff(range(x))

  got <- readBaseImage(file, 'cutout', NULL, TRUE, 'test')
  expect_equal(got, rescale(stored[, , 1]))

  # And the old behaviour is genuinely different here, unlike the opaque case
  # above: averaging alpha in survives the rescaling when alpha is not constant.
  folded <- rescale(apply(stored, c(1, 2), mean))
  expect_false(isTRUE(all.equal(got, folded)))
  expect_gt(max(abs(got - folded)), 0.1)
})
