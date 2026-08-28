#' Build a generator that renders through an external program
#'
#' Wrap any program that can turn latent vectors into images, so a generative
#' adversarial network running under Python can be used as the renderer for the
#' experimental latent-space pipeline.
#'
#' rcicr never loads the model. It writes the latents to a file, runs the
#' program, and reads back the images the program wrote. Anything that can be
#' invoked from a command line qualifies, whatever language it is written in and
#' whatever it needs installed, and rcicr gains no dependency from it.
#'
#' \subsection{The protocol}{
#'
#' rcicr calls
#' \code{<command> <args> <latents file> <output directory>} and waits.
#'
#' The latents file holds one comma-separated row per stimulus, with no header
#' and no row names, \code{latent_dim} numbers per row.
#'
#' The program must write one PNG per row into the output directory, named
#' \code{00001.png}, \code{00002.png} and so on in the order of the rows, each
#' \code{img_size} by \code{img_size} pixels. Colour images are averaged to
#' greyscale on the way back in.
#'
#' A script under \code{inst/python} in the installed package shows the whole of
#' it for StyleGAN; \code{system.file('python', 'rcicr_stylegan.py', package =
#' 'rcicr')} finds it.
#' }
#'
#' \subsection{Identifying the model}{
#'
#' A stimulus file records a fingerprint so a classification image cannot later
#' be rendered through a different model than the one participants saw. rcicr
#' cannot see inside this program, so the fingerprint has to come from you.
#' Passing \code{weights} is the reliable way: the file's own checksum changes
#' whenever the model does. Without it the fingerprint covers only the command
#' and the dimensions, which will not notice a model swapped underneath the same
#' script.
#' }
#'
#' @export
#' @param command The program to run: an executable name on the path, or a path
#'   to one.
#' @param latent_dim Number of dimensions in one latent vector.
#' @param img_size Pixel size of the images the program writes. They must be
#'   square.
#' @param args Character vector of arguments passed to the program before the
#'   latents file and the output directory.
#' @param space String naming the space the latents live in, recorded for
#'   reference: \code{z}, \code{w} or \code{w+} for a StyleGAN.
#' @param latent_mean The centre of the space, length \code{latent_dim}.
#'   Defaults to NULL, which uses zero. For a StyleGAN's W space, pass the
#'   average latent rather than zero.
#' @param latent_sd The spread of the space along each dimension, length
#'   \code{latent_dim} or a single number applied to all of them. This is what
#'   perturbation sizes are expressed in, so a wrong value makes
#'   \code{latent_sigma} mean something other than it says.
#' @param weights Path to the model's weight file, used to fingerprint it.
#'   Defaults to NULL. See Identifying the model above.
#' @param fingerprint String identifying this generator, if you would rather set
#'   it yourself. Defaults to NULL, which derives one.
#' @param timeout Seconds to wait for the program before giving up. Zero waits
#'   forever.
#' @return An \code{rcicr_generator} object, as \code{\link{latentGeneratorPCA}}
#'   returns.
#' @examples
#' # This function is part of the experimental latent-space module.
#' options(rcicr.experimental = TRUE)
#'
#' # Any program will do. This one ignores the latents and writes grey squares,
#' # which is enough to show the protocol.
#' script <- tempfile(fileext = ".R")
#' writeLines(c(
#'   'args <- commandArgs(trailingOnly = TRUE)',
#'   'latents <- as.matrix(read.csv(args[1], header = FALSE))',
#'   'for (i in seq_len(nrow(latents))) {',
#'   '  png::writePNG(matrix(0.5, 8, 8), file.path(args[2], sprintf("%05d.png", i)))',
#'   '}'
#' ), script)
#'
#' generator <- latentGeneratorCommand(
#'   command = file.path(R.home("bin"), "Rscript"),
#'   args = script, latent_dim = 4, img_size = 8, latent_sd = 1
#' )
#' generator$latent_dim
latentGeneratorCommand <- function(command, latent_dim, img_size, args = character(), space = 'w', latent_mean = NULL, latent_sd = 1, weights = NULL, fingerprint = NULL, timeout = 0) {

  requireExperimental('latentGeneratorCommand')

  latent_dim <- as.integer(latent_dim)
  img_size <- as.integer(img_size)

  if (is.null(latent_mean)) {
    latent_mean <- rep(0, latent_dim)
  }
  if (length(latent_sd) == 1) {
    latent_sd <- rep(latent_sd, latent_dim)
  }

  if (is.null(fingerprint)) {
    fingerprint <- commandFingerprint(command, args, latent_dim, img_size, space, weights)
  }

  render <- function(latents) {
    return(renderThroughCommand(latents, command, args, img_size, timeout))
  }

  generator <- rcicrGenerator(
    kind = 'command',
    latent_dim = latent_dim,
    img_size = img_size,
    space = space,
    latent_mean = latent_mean,
    latent_sd = latent_sd,
    render = render,
    fingerprint = fingerprint,
    state = list(command = command, args = args, weights = weights)
  )

  validateGenerator(generator)

  return(generator)
}

renderThroughCommand <- function(latents, command, args, img_size, timeout) {
  workdir <- tempfile('rcicr_render')
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  outdir <- file.path(workdir, 'images')
  dir.create(outdir)
  latents_file <- file.path(workdir, 'latents.csv')
  utils::write.table(latents, latents_file, sep = ',', row.names = FALSE,
                     col.names = FALSE)

  log_file <- file.path(workdir, 'output.log')
  status <- system2(command, c(args, shQuote(latents_file), shQuote(outdir)),
                    stdout = log_file, stderr = log_file, timeout = timeout)

  if (!identical(status, 0L)) {
    msg <- paste0(
      'The rendering program failed (exit status ', status, '). Its output was:\n',
      paste(utils::tail(readLines(log_file, warn = FALSE), 20), collapse = '\n')
    )
    stop(msg, call. = FALSE)
  }

  return(readRenderedImages(outdir, nrow(latents), img_size, log_file))
}

readRenderedImages <- function(outdir, n_latents, img_size, log_file) {
  out <- array(0, dim = c(n_latents, img_size, img_size))

  for (i in seq_len(n_latents)) {
    file <- file.path(outdir, sprintf('%05d.png', i))

    if (!file.exists(file)) {
      msg <- paste0(
        'The rendering program wrote no image for latent ', i, '. It must write ',
        n_latents, ' PNGs named 00001.png upwards into the output directory it ',
        'is given. Its output was:\n',
        paste(utils::tail(readLines(log_file, warn = FALSE), 20), collapse = '\n')
      )
      stop(msg, call. = FALSE)
    }

    img <- png::readPNG(file)
    if (length(dim(img)) == 3) {
      img <- apply(img[, , seq_len(min(3, dim(img)[3])), drop = FALSE], c(1, 2), mean)
    }

    if (nrow(img) != img_size || ncol(img) != img_size) {
      msg <- paste0(
        'The rendering program wrote a ', nrow(img), ' by ', ncol(img),
        ' image for latent ', i, ', but this generator was built with img_size = ',
        img_size, '. rcicr does not resize.'
      )
      stop(msg, call. = FALSE)
    }

    out[i, , ] <- img
  }

  return(out)
}

# Without the weight file there is nothing here that changes when the model
# does, so the fingerprint catches a different script or different dimensions
# and not a different set of weights behind the same script.
commandFingerprint <- function(command, args, latent_dim, img_size, space, weights) {
  parts <- c(command, args, latent_dim, img_size, space)

  if (!is.null(weights)) {
    if (!file.exists(weights)) {
      stop(paste0('The weights file does not exist: ', weights, '.'), call. = FALSE)
    }
    parts <- c(parts, unname(tools::md5sum(weights)))
  }

  return(paste0('command:', paste(parts, collapse = '|')))
}
