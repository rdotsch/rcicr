#' Plots a Z-map
#'
#' Plots a Z-map given a matrix of z-scores that maps onto a specified base image.
#'
#' This function takes in a matrix of z-scores (as returned by generateCI) and an Rdata file containing a base image. It returns a Z-map image in PNG format.
#' Unlisted additional arguments will be passed to raster::plot. For example, a different color palette can be specified using the \code{col} argument. See raster::plot for details.
#'
#' @section Reproducibility across platforms:
#' The z-scores themselves are ordinary R arithmetic and do not depend on your operating system, and neither do classification images, scaling or informational value. The PNG this function writes is different: it is drawn through a graphics device, and graphics devices differ between platforms both in colour management and in whether they write an alpha channel. The same z-map rendered on Linux and on macOS gives visibly identical figures whose files are not byte-identical -- macOS renders a mid-grey background at roughly 0.573 where the cairo device gives 0.502.
#'
#' So when you are checking that an analysis reproduces, compare the numbers rather than the rendered figures. A z-map PNG that differs pixel-for-pixel on a colleague's machine is not a different result, and regenerating figures on another platform is safe.
#'
#' This applies only to \code{plotZmap()}, which is the only function in the package that opens a graphics device. Every other PNG written by \code{rcicr} -- stimuli, classification images, autoscaled classification images -- is written straight from the pixel array by \code{png::writePNG()} and carries no such dependence.
#'
#' @export
#' @importFrom raster raster plot
#' @importFrom grDevices png dev.off
#' @importFrom graphics rasterImage par plot.new plot.window
#' @param zmap A matrix containing z-scores that map onto a given base image. zmap and baseimage must have the same dimensions.
#' @param bgimage A matrix containing the grayscale image to use as a background. This should be either the base image or the final CI. If not this argument is not given, only the Z-map will be drawn.
#' @param sigma The sigma of the smoothing that was applied to the CI to create the Z-map.
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param mask Optional. A binary matrix with the same dimensions as zmap: cells that are 0 (or FALSE) are masked, cells that are 1 (or TRUE) are kept. Can also be the filename (as a string) of a black and white PNG image, in which case black (0) is masked and white (1) is kept. This is the same convention \code{generateCI()} uses for its own \code{mask} argument, so a mask can be used with both. Note that earlier versions of this documentation stated the opposite for the matrix form; the description here matches what the code does.
#' @param decoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param targetpath String specifying the directory to save the Z-map PNG to. Required: this function exists to write a file, and has no default path. It is created if it does not exist. Use \code{tempdir()} if you only want to try the function out.
#' @param filename Optional string to specify a file name for the Z-map PNG.
#' @param size Integer specifying the width and height of the PNG image (default: 512).
#' @param ... Additional arguments to be passed to raster::plot. Only applied when decoration is TRUE.
#' @return Nothing. It writes a Z-map image.
#' @examples
#' set.seed(1)
#' zmap <- matrix(rnorm(64, sd = 5), 8, 8)
#' plotZmap(zmap, sigma = 3, threshold = 3, decoration = FALSE,
#'          targetpath = tempdir(), size = 200)
plotZmap <- function(zmap, bgimage = '', sigma, threshold = 3, mask = NULL, decoration = TRUE, targetpath, filename = 'zmap', size = 512, ...) {

  # targetpath is required, not defaulted: a default path writes to the user's
  # filespace uninvited, which CRAN policy does not allow.
  if (missing(targetpath)) {
    stop(paste0('No targetpath was given. plotZmap() writes a PNG, so it needs a ',
                'directory to write it to: supply targetpath = <a directory>. ',
                'Use tempdir() if you only want to try the function out.'))
  }

  # Create target directory
  dir.create(targetpath, recursive = TRUE, showWarnings = FALSE)

  # If a mask is specified, import and check it
  if (!(is.null(mask))) {
    # Read in the mask from a PNG image if specified
    if (!is.matrix(mask)) {
      mask <- png::readPNG(mask)
    }
    # Are mask and zmap the same size?
    if (nrow(zmap) == dim(mask)[1] & ncol(zmap) == dim(mask)[2]) {
      # Are all the values either 0/1, or TRUE/FALSE?
      if (all(mask %in% c(0, 1)) | all(mask %in% c(TRUE, FALSE))) {
        # If we have more than 1 layer (i.e. the PNG was not greyscale but RGB or
        # CMYK), are all the layers identical? If so, remove superfluous layers
        if (length(dim(mask)) != 2) {
          iden <- c()
          for (i in 2:dim(mask)[3]) {
            if (identical(mask[,,i-1], mask[,,i])) {
              iden <- c(iden, TRUE)
            } else {
              iden <- c(iden, FALSE)
            }
          }
          if (all(iden)) {
            mask <- mask[,,1]
          } else {
            stop('Error in importing Z-map mask: color channels are not identical.')
          }
        }
      } else {
        stop('Error in importing Z-map mask: pixel values are not limited to black (0) and white (1).')
      }
    } else {
      stop('Error in importing Z-map mask: mask and Z-map are not the same size.')
    }
    # Convert to boolean: 0 (black, FALSE) marks the region to mask, for a
    # matrix and a PNG alike.
    #
    # The single convention is deliberate and was verified against the working
    # sibling rather than taken from the docs: generateCI()'s applyMask() does
    # `mask_matrix == 0` unconditionally (R/generateCI.R), so the same mask
    # passed to generateCI() and plotZmap() must mask the same region. Both
    # functions' roxygen claimed 1/TRUE meant masked for a matrix while black
    # (0) meant masked for a PNG -- two opposite conventions in one sentence,
    # and contradicted by the only implementation that ever ran. The docs were
    # corrected to match the code, not the other way round: applyMask()'s
    # behaviour is what users' existing masks were built against.
    #
    # This conversion used to be two in-place assignments --
    # `mask[mask == 0] <- TRUE` then `mask[mask == 1] <- FALSE` -- which set
    # *every* cell to FALSE whatever the input: the first assignment coerces
    # TRUE to 1, so the second matches the cells it had just set as well as the
    # originally-true ones. A swap without a temporary.
    mask <- mask == 0
  }

  # Apply threshold
  zmap[abs(zmap) < threshold] <- NA

  # Apply the mask. Masked cells drop out of the z-map exactly as sub-threshold
  # cells do, which is what the documentation has promised since 2016 --
  # commit 18e07cb landed the import half as "add mask import ... (todo:
  # applying the mask)" and the todo was never picked up, so the argument was
  # validated and then silently discarded in every released version.
  # The mask was validated and then never applied; see NEWS.md 1.2.0.
  if (!is.null(mask)) {
    zmap[mask] <- NA
  }

  # Plot
  png(filename = paste0(targetpath, '/', filename, '.png'), width = size, height = size)

  # Close the device however this function exits, so a failure while plotting
  # cannot leak it or leave a half-written PNG behind.
  on.exit(dev.off())

  # With decoration
  if (decoration) {
    # Initial (dummy) plot; sets up plot with initial dimensions + scale, title, label
    raster::plot(raster(zmap), axes = FALSE, box = FALSE, main = paste0('Z-map of ', filename),
                 xlab = paste0('sigma = ', sigma, ', threshold = ', threshold),
                 col = viridis::viridis(100), ...)
    # Add bgimage (if specified) and superimpose Z-map on top of it
    if (!(identical(bgimage, ''))) {
      rasterImage(bgimage, 0, 0, 1, 1)
      raster::plot(raster(zmap), add = TRUE, ...)
    }
    # If no bgimage was specified, draw a boundary box around the Z-map
    if (identical(bgimage, '')) {
      box <- matrix(NA, nrow(zmap) + 1, ncol(zmap) + 1)
      box[c(1, nrow(zmap) + 1), ] <- 0
      box[, c(1, ncol(zmap) + 1)] <- 0
      rasterImage(box, 0, 0, 1, 1)
    }
    # Without decoration
  }
  if (!decoration) {
    # Initialize plot without margins. The order matters: plot.new() validates
    # that the current margins fit inside the device, and the default margins
    # (c(5.1, 4.1, 4.1, 2.1) *lines*) do not fit a small one. Setting them after
    # plot.new(), as this did, meant any device below roughly 100 px failed with
    # "figure margins too large" -- and generateCI() sizes the device to
    # img_size, so small stimulus sets could not produce a z-map at all.
    # Rendered output at usual sizes is unchanged: par(mar=) still takes effect
    # before plot.window() either way.
    #
    # par() returns the previous values of exactly the parameters being set, so
    # restoring `oldpar` restores `mar` and nothing else. par(no.readonly=TRUE)
    # is not usable here: it also captures derived parameters such as `pin`,
    # which plot.window() below invalidates, and restoring them errors.
    oldpar <- par(mar = c(0, 0, 0, 0))
    on.exit(par(oldpar), add = TRUE, after = FALSE)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')

    # If specified, add bgimage. Must be identical(), not !=: bgimage is normally
    # an image matrix, so `bgimage != ''` is a condition of length img_size^2.
    # R >= 4.2 makes that an error rather than silently using the first element,
    # so this branch could not run at all with a background image -- which is
    # every call from generateCI(), as it always passes the combined CI. The
    # decoration = TRUE branch above already used identical(); this one was
    # missed. Same root cause as the `mask` bug fixed in 1.1.0.
    if (!identical(bgimage, '')) {
      rasterImage(bgimage, 0, 0, 1, 1)
    }
    # Add Z-map
    raster::plot(raster(zmap), add = TRUE, legend = FALSE)
  }
  # The device is closed by the on.exit() handler registered after png().
}
