#' Renders latents into images through a generator
#'
#' Turn one latent vector, or a matrix of them, into the greyscale images a
#' generator produces for them.
#'
#' Every other function in the experimental latent-space module calls this to
#' get from latents to pixels. It is exported because a caller holding a
#' generator has the same need: to look at the face at the centre of the space
#' before starting, to render a latent a search arrived at, or to see what a
#' particular direction does.
#'
#' Latents are given one per row. A generator renders a whole matrix in one
#' call, because an external one pays its startup cost per call, so rendering
#' many latents together is much faster than looping.
#'
#' @export
#' @param generator A generator built by \code{\link{latentGeneratorPCA}} or
#'   another of rcicr's latentGenerator functions.
#' @param latents A matrix of latents, one row per image, or a plain vector for
#'   a single one.
#' @param validate Boolean specifying whether the generator is checked against
#'   its contract first (default: TRUE). The functions in this module set it to
#'   FALSE inside a render loop, having checked once already.
#' @return An array of \code{nrow(latents)} by \code{img_size} by
#'   \code{img_size} greyscale pixel values in \code{[0, 1]}. A single latent
#'   still comes back as a three-dimensional array, so \code{[1, , ]} is the
#'   image.
#' @examples
#' # This function is part of the experimental latent-space module.
#' options(rcicr.experimental = TRUE)
#'
#' faces <- replicate(6, tempfile(fileext = ".png"))
#' for (i in seq_along(faces)) {
#'   png::writePNG(matrix(runif(16 * 16), 16, 16), faces[i])
#' }
#' generator <- latentGeneratorPCA(faces, n_components = 3, img_size = 16)
#'
#' # The face at the centre of the space: for an eigenface generator, the
#' # average of the images it was built from.
#' mean_face <- renderLatent(generator, generator$latent_mean)[1, , ]
#' dim(mean_face)
#'
#' # Two standard deviations along the first component, and back the other way.
#' along <- rbind(c(2, 0, 0), c(-2, 0, 0)) * generator$latent_sd[1]
#' dim(renderLatent(generator, along))
#' @name renderLatent
NULL
