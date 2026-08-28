#' Build an eigenface generator from a set of face images
#'
#' Constructs a generative face model in base R, for use as the renderer in the
#' experimental latent-space reverse correlation pipeline.
#'
#' The model is an ordinary principal component analysis of the images: the
#' latent space is spanned by the leading eigenfaces, the origin is the mean
#' face, and rendering a latent is the mean face plus that latent's weighted sum
#' of components, clamped to \code{[0, 1]}. It is linear, so it is a far weaker
#' face model than a generative adversarial network, and the faces it renders
#' look like blurred averages rather than photographs.
#'
#' It exists so the method does not depend on a GPU. Every function in the
#' latent module can be run, tested and checked against this generator with no
#' Python, no network and no additional package, which is what makes the rest of
#' the module verifiable. For hyper-realistic stimuli, supply a StyleGAN through
#' one of the other backends instead; nothing else about the pipeline changes.
#'
#' The number of components is capped at one fewer than the number of images,
#' because that is the rank of the centred data: 40 images can support at most
#' 39 components however many are asked for.
#'
#' @export
#' @param base_face_files Character vector of image file paths, or a named list
#'   of them as \code{\link{generateStimuli2IFC}} takes. Accepts JPEG and PNG
#'   images, recognised by a \code{.png}, \code{.jpg} or \code{.jpeg} extension.
#'   Each image must be square, and all of them the same size. Two images is the
#'   minimum; a useful face space needs tens of them, aligned on the eyes, as
#'   any eigenface model does.
#' @param n_components Number of components to keep, and so the dimensionality
#'   of the latent space. Capped at \code{length(base_face_files) - 1}.
#' @param img_size Number specifying the pixel size the images are expected to
#'   be. Defaults to NULL, which takes the size from the first image and then
#'   requires every other image to match it. rcicr does not resize images.
#' @param maximize_baseimage_contrast Boolean specifying whether each image's
#'   pixel values are rescaled to maximize its contrast before the analysis,
#'   matching \code{\link{generateStimuli2IFC}}. An image with no contrast at
#'   all is rejected with an error rather than becoming an all-\code{NaN} row.
#' @return An \code{rcicr_generator} object: a list whose \code{render} element
#'   turns a matrix of latents, one row per stimulus, into an array of
#'   \code{nrow} by \code{img_size} by \code{img_size} greyscale pixel values.
#'   Its \code{latent_sd} element gives the standard deviation of the training
#'   faces along each component, so a perturbation size is expressed in units of
#'   how much real faces vary.
#' @examples
#' # This function is part of the experimental latent-space module.
#' options(rcicr.experimental = TRUE)
#'
#' # synthetic square greyscale images stand in for a set of aligned face photos
#' faces <- replicate(6, tempfile(fileext = ".png"))
#' for (i in seq_along(faces)) {
#'   png::writePNG(matrix(runif(16 * 16), 16, 16), faces[i])
#' }
#'
#' generator <- latentGeneratorPCA(faces, n_components = 3, img_size = 16)
#' generator$latent_dim
latentGeneratorPCA <- function(base_face_files, n_components = 50, img_size = NULL, maximize_baseimage_contrast = TRUE) {

  requireExperimental('latentGeneratorPCA')

  files <- normalizeImagePaths(base_face_files)

  if (length(files) < 2) {
    msg <- paste0(
      'latentGeneratorPCA() needs at least 2 images to have ',
      'anything to vary: a single face has no variation around itself. It ',
      'was given ', length(files), '.'
    )
    stop(msg, call. = FALSE)
  }

  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    msg <- paste0(
      'These image file(s) do not exist: ', paste(missing_files, collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  unreadable <- files[is.na(vapply(files, baseImageFormat, character(1)))]
  if (length(unreadable) > 0) {
    msg <- paste0(
      'These file(s) are neither PNG nor JPEG, by their ', 'extension: ',
      paste(unreadable, collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  # The size check applies to the first image too, so a mismatch anywhere in the
  # set is reported the same way rather than the first image silently setting a
  # size the rest fail against.
  if (is.null(img_size)) {
    img_size <- nrow(readBaseImage(files[1], names(files)[1], NULL, FALSE, 'latentGeneratorPCA'))
  }

  faces <- readBaseImages(
    stats::setNames(as.list(files), names(files)),
    img_size, maximize_baseimage_contrast, 'latentGeneratorPCA'
  )

  # One row per face, one column per pixel. as.vector() reads a matrix in
  # column-major order, which is the order render() reverses below.
  pixels_by_face <- t(vapply(faces, as.vector, numeric(img_size * img_size)))

  mean_face <- colMeans(pixels_by_face)
  centred <- sweep(pixels_by_face, 2, mean_face, '-')

  # Against a tolerance, never against zero. Centring identical rows does not
  # leave exactly zero: the mean of three copies of a number is not always that
  # number in floating point, so x - mean(c(x, x, x)) comes out around 1e-17 for
  # some pixel values and exactly 0 for others, and which happens depends on the
  # summation order the platform uses. Testing for zero here passed on Linux and
  # failed on macOS ARM64, and so did reading the same thing off the singular
  # values, whose threshold is relative to the largest of them and so collapses
  # to zero in exactly this case.
  #
  # The scale is the pixel values themselves, which is the one quantity in this
  # calculation that a degenerate set does not send to zero.
  no_variation <- .Machine$double.eps * max(abs(pixels_by_face)) * max(dim(centred))

  if (max(abs(centred)) <= no_variation) {
    msg <- paste0(
      'These images have no variation between them: every one ',
      'is identical, so there is no face space to build.'
    )
    stop(msg, call. = FALSE)
  }

  decomposition <- svd(centred)

  # The centred data has rank at most n - 1, so components beyond that are noise
  # with a singular value of zero and a latent_sd of zero, which would make
  # perturbing along them a no-op. The tolerance is the conventional one for a
  # numerical rank: the largest singular value, scaled by the larger dimension.
  tolerance <- .Machine$double.eps * max(dim(centred)) * decomposition$d[1]
  keep <- min(n_components, nrow(pixels_by_face) - 1, sum(decomposition$d > tolerance))

  components <- decomposition$v[, seq_len(keep), drop = FALSE]

  # Latents are in raw score units, so render() is a plain matrix product, and
  # latent_sd reports how far the training faces spread along each component.
  # A caller's sigma is then read in units of real face variation.
  latent_sd <- decomposition$d[seq_len(keep)] / sqrt(nrow(pixels_by_face) - 1)

  generator <- pcaGenerator(mean_face, components, latent_sd, img_size,
                            nrow(pixels_by_face), files)

  validateGenerator(generator)

  return(generator)
}

# Both the named-list form generateStimuli2IFC() takes and a plain character
# vector are accepted, so a caller with a directory listing does not have to
# name every file. Names are only used in error messages here.
normalizeImagePaths <- function(base_face_files) {
  if (is.list(base_face_files)) {
    lengths_ok <- vapply(base_face_files, function(f) is.character(f) && length(f) == 1, logical(1))
    if (!all(lengths_ok)) {
      msg <- paste0(
        'Each element of base_face_files must be a single file ', 'name. Element(s) ',
        paste(which(!lengths_ok), collapse = ', '), ' are not.'
      )
      stop(msg, call. = FALSE)
    }
    base_face_files <- unlist(base_face_files)
  }

  if (!is.character(base_face_files)) {
    msg <- paste0(
      'base_face_files must be a character vector of image file ',
      'names, or a named list of them. It is of class ',
      paste(class(base_face_files), collapse = '/'), '.'
    )
    stop(msg, call. = FALSE)
  }

  if (is.null(names(base_face_files))) {
    names(base_face_files) <- basename(base_face_files)
  }

  return(base_face_files)
}

# The generator itself, separated from the analysis that produces it so a
# stimulus file can rebuild the identical renderer without the original images.
pcaGenerator <- function(mean_face, components, latent_sd, img_size, n_faces, base_face_files) {
  render <- function(latents) {
    pixels <- sweep(latents %*% t(components), 2, mean_face, '+')

    # An eigenface reconstruction leaves [0, 1] readily, and the contract is
    # that a generator hands back displayable pixels.
    pixels[pixels < 0] <- 0
    pixels[pixels > 1] <- 1

    out <- array(0, dim = c(nrow(latents), img_size, img_size))
    for (i in seq_len(nrow(latents))) {
      out[i, , ] <- matrix(pixels[i, ], img_size, img_size)
    }

    return(out)
  }

  keep <- ncol(components)

  return(rcicrGenerator(
    kind = 'pca',
    latent_dim = keep,
    img_size = img_size,
    space = 'pca',
    latent_mean = rep(0, keep),
    latent_sd = latent_sd,
    render = render,
    fingerprint = paste0('pca:', fingerprintNumeric(mean_face, components, latent_sd)),
    state = list(
      mean_face = mean_face,
      components = components,
      n_faces = n_faces,
      base_face_files = base_face_files
    )
  ))
}
