#' Plots a Z-map
#'
#' Plots a Z-map given a matrix of z-scores that maps onto a specified base image.
#'
#' This function takes in a matrix of z-scores (as returned by generateCI) and an Rdata file containing a base image. It returns a Z-map image in PNG format.
#' Unlisted additional arguments will be passed to raster::plot. For example, a different color palette can be specified using the \code{col} argument. See raster::plot for details.
#'
#' @export
#' @import dplyr
#' @importFrom raster raster plot
#' @importFrom grDevices png
#' @importFrom graphics rasterImage par plot.new plot.window
#' @param zmap A matrix containing z-scores that map onto a given base image. zmap and baseimage must have the same dimensions.
#' @param bgimage A matrix containing the grayscale image to use as a background. This should be either the base image or the final CI. If not this argument is not given, only the Z-map will be drawn.
#' @param sigma The sigma of the smoothing that was applied to the CI to create the Z-map.
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param mask Optional. A boolean matrix with the same dimensions as zmap. If a cell evaluates to TRUE, the corresponding zmap pixel will be masked. Can also be the filename (as a string) of a black and white PNG image to be converted to a matrix (black = masked).
#' @param decoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param targetpath String specifying path to save the Z-map PNG to.
#' @param filename Optional string to specify a file name for the Z-map PNG.
#' @param size Integer specifying the width and height of the PNG image (default: 512).
#' @param ... Additional arguments to be passed to raster::plot. Only applied when decoration is TRUE.
#' @return Nothing. It writes a Z-map image.
plotZmap <- function(zmap, bgimage = '', sigma, threshold = 3, mask = NULL, decoration = T, targetpath = 'zmaps', filename = 'zmap', size = 512, ...) {

  # Create target directory
  dir.create(targetpath, recursive = T, showWarnings = F)

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
    # Convert to boolean
    mask[mask == 0] <- TRUE
    mask[mask == 1] <- FALSE
  }

  # Apply threshold
  zmap[abs(zmap) < threshold] <- NA

  # Plot
  png(filename = paste0(targetpath, '/', filename, '.png'), width = size, height = size)

  # With decoration
  if (decoration) {
    # Initial (dummy) plot; sets up plot with initial dimensions + scale, title, label
    raster::plot(raster(zmap), axes = F, box = F, main = paste0('Z-map of ', filename),
                 xlab = paste0('sigma = ', sigma, ', threshold = ', threshold),
                 col = viridis::viridis(100), ...)
    # Add bgimage (if specified) and superimpose Z-map on top of it
    if (!(identical(bgimage, ''))) {
      rasterImage(bgimage, 0, 0, 1, 1)
      raster::plot(raster(zmap), add = T, ...)
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
    # Initialize plot without margins
    plot.new()
    par(mar = c(0, 0, 0, 0))
    plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')

    # If specified, add bgimage
    if (bgimage != '') {
      rasterImage(bgimage, 0, 0, 1, 1)
    }
    # Add Z-map
    raster::plot(raster(zmap), add = T, legend = F)
  }
  dev.off()
}