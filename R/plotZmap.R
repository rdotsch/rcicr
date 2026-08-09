# The palette raster::plot() used when none was given, reproduced here so that
# dropping the dependency does not change what a z-map looks like. It came from
# the RasterLayer plot method's `if (missing(col)) col <- rev(terrain.colors(255))`.
zmapDefaultPalette <- function() rev(terrain.colors(255))

# Draws a z-map matrix the way raster::plot(raster(zmap)) drew it: row 1 at the
# top, column 1 at the left, over the extent 0..1 that rasterImage() paints the
# background into. image() lays matrix rows along x and counts y upward, so the
# matrix needs transposing and row-reversing.
#
# The coordinates are cell EDGES (n + 1 of them), not centres. Given centres,
# image() infers the edges by extending half a cell past the first and last, so
# the z-map would cover slightly more than 0..1 and sit shifted against the
# background drawn underneath it.
drawZmapLayer <- function(zmap, col, add = FALSE, ...) {
  z <- t(zmap[nrow(zmap):1, , drop = FALSE])
  # A fully masked or entirely sub-threshold z-map is all NA, where range() of
  # nothing gives c(Inf, -Inf) and a warning from each end. raster::plot() drew
  # that empty map silently; pinning zlim keeps the call well-defined.
  finite <- z[is.finite(z)]
  # ylab is pinned empty because image() otherwise labels the axis with a
  # deparse of whatever expression was handed to y.
  args <- list(x = seq(0, 1, length.out = nrow(z) + 1),
               y = seq(0, 1, length.out = ncol(z) + 1),
               z = z, col = col,
               zlim = if (length(finite)) range(finite) else c(0, 1),
               add = add, useRaster = TRUE, ylab = '')

  # Merged by name, not passed alongside: R stops with `formal argument "zlim"
  # matched by multiple actual arguments` if a caller names something already
  # set above. This covers every default at once.
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(image, args)
}

# The colour bar raster::plot() drew beside a decorated map. Nothing in base R
# draws one, and fields::image.plot would only trade one dependency for another.
drawZmapLegend <- function(zmap, col, zlim = NULL, breaks = NULL) {
  finite <- zmap[is.finite(zmap)]

  # A bar wherever there is a scale to label. raster::plot() drew one for a
  # constant z-map and for an all-NA map given an explicit zlim, and drew none
  # only when neither a caller scale nor a surviving value existed.
  if (is.null(breaks) && is.null(zlim) && !length(finite)) {
    return(invisible(NULL))
  }
  # Labelled with the scale the map was drawn on, so the order matches image()'s
  # own precedence: caller's breaks, then caller's zlim, then the data.
  edges <- if (!is.null(breaks)) {
    breaks
  } else {
    if (is.null(zlim)) zlim <- range(finite)
    # A constant z-map gives a degenerate range, which seq() would collapse to a
    # bar of zero height. image() widens such a range around the value rather
    # than refusing it; the same rule is applied here.
    if (diff(zlim) == 0) {
      zlim <- if (zlim[1] == 0) c(-1, 1) else zlim[1] + c(-0.4, 0.4) * abs(zlim[1])
    }
    seq(zlim[1], zlim[2], length.out = length(col) + 1)
  }

  oldpar <- par(fig = c(0.86, 0.90, 0.30, 0.70), mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(par(oldpar), add = TRUE, after = FALSE)

  # No useRaster here, unlike the map: unequally spaced breaks make these y
  # coordinates an irregular grid, which it refuses outright. A single column of
  # rectangles gains nothing from rasterising anyway.
  image(x = c(0, 1), y = edges,
        z = matrix(seq_along(col), nrow = 1), col = col,
        axes = FALSE, xlab = '', ylab = '')
  axis(4, las = 1, cex.axis = 0.8, tcl = -0.2, mgp = c(3, 0.4, 0))
  invisible(NULL)
}

#' Plots a Z-map
#'
#' Plots a Z-map given a matrix of z-scores that maps onto a specified base image.
#'
#' This function takes in a matrix of z-scores (as returned by generateCI) and an Rdata file containing a base image. It returns a Z-map image in PNG format.
#' Unlisted additional arguments will be passed to \code{graphics::image}. For example, a different color palette can be specified using the \code{col} argument. See \code{graphics::image} for details. Versions up to and including 1.2.3 passed these to the raster package's plot method instead; \code{col} works the same way in both, but an argument specific to that method will no longer be understood.
#'
#' @section Reproducibility across platforms:
#' The z-scores themselves are ordinary R arithmetic and do not depend on your operating system, and neither do classification images, scaling or informational value. The PNG this function writes is different: it is drawn through a graphics device, and graphics devices differ between platforms both in colour management and in whether they write an alpha channel. The same z-map rendered on Linux and on macOS gives visibly identical figures whose files are not byte-identical -- macOS renders a mid-grey background at roughly 0.573 where the cairo device gives 0.502.
#'
#' So when you are checking that an analysis reproduces, compare the numbers rather than the rendered figures. A z-map PNG that differs pixel-for-pixel on a colleague's machine is not a different result, and regenerating figures on another platform is safe.
#'
#' This applies only to \code{plotZmap()}, which is the only function in the package that opens a graphics device. Every other PNG written by \code{rcicr} -- stimuli, classification images, autoscaled classification images -- is written straight from the pixel array by \code{png::writePNG()} and carries no such dependence.
#'
#' @export
#' @importFrom grDevices png dev.off terrain.colors
#' @importFrom graphics rasterImage par plot.new plot.window image axis
#' @param zmap A matrix containing z-scores that map onto a given base image. zmap and baseimage must have the same dimensions.
#' @param bgimage A matrix containing the grayscale image to use as a background. This should be either the base image or the final CI. If not this argument is not given, only the Z-map will be drawn.
#' @param sigma The sigma of the smoothing that was applied to the CI to create the Z-map.
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param mask Optional. A binary matrix with the same dimensions as zmap: cells that are 0 (or FALSE) are masked, cells that are 1 (or TRUE) are kept. Can also be the filename (as a string) of a black and white PNG image, in which case black (0) is masked and white (1) is kept. This is the same convention \code{generateCI()} uses for its own \code{mask} argument, so a mask can be used with both. Note that earlier versions of this documentation stated the opposite for the matrix form; the description here matches what the code does.
#' @param decoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param targetpath String specifying the directory to save the Z-map PNG to. Required: this function exists to write a file, and has no default path. It is created if it does not exist. Use \code{tempdir()} if you only want to try the function out.
#' @param filename Optional string to specify a file name for the Z-map PNG.
#' @param size Integer specifying the width and height of the PNG image (default: 512).
#' @param pointsize Integer specifying the text size of the decoration, in points (default: 12, the graphics device's own default). Margins are measured in lines of text, so this also sets how much of the image the decoration takes up: a decorated z-map needs roughly \code{12.3 * pointsize} pixels, about 160px at the default. Below that \code{plotZmap()} stops and says so. Lowering it is what makes a decorated z-map possible on a small device -- at the cost of a smaller map, since the margins shrink but the labels still need room. Ignored when \code{decoration = FALSE}, which has no margins and works at any size.
#' @param ... Additional arguments to be passed to \code{graphics::image}. Only applied when decoration is TRUE.
#' @return Nothing. It writes a Z-map image.
#' @examples
#' set.seed(1)
#' zmap <- matrix(rnorm(64, sd = 5), 8, 8)
#' plotZmap(zmap, sigma = 3, threshold = 3, decoration = FALSE,
#'          targetpath = tempdir(), size = 200)
# pointsize sits after ... so it can only be supplied by name: placed before it,
# a tenth positional argument would bind here instead of reaching image().
plotZmap <- function(zmap, bgimage = '', sigma, threshold = 3, mask = NULL, decoration = TRUE, targetpath, filename = 'zmap', size = 512, ..., pointsize = 12) {

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
  outfile <- paste0(targetpath, '/', filename, '.png')
  png(filename = outfile, width = size, height = size, pointsize = pointsize)

  # Close the device however this function exits, so a failure while plotting
  # cannot leak it. png() creates the file when it opens, so an unfinished
  # render leaves a stub behind unless it is removed as well -- the fit check
  # below stops between those two points.
  drawn <- FALSE
  on.exit({
    dev.off()
    if (!drawn) unlink(outfile)
  })

  # With decoration
  if (decoration) {
    zmapMargins <- c(5.1, 4.1, 4.1, 6.1)

    # Margins are measured in lines of text, so what the decoration needs scales
    # with pointsize and can exceed the device: base R then stops in plot.new()
    # with `figure margins too large`, naming neither this package nor the way
    # out. Asking the open device for its own line height and size gets the
    # boundary exactly right -- 150px fits at pointsize 12 and 140px does not --
    # and stays correct if the margins above ever change.
    needed <- c(sum(zmapMargins[c(2, 4)]), sum(zmapMargins[c(1, 3)])) * par('csi')
    if (any(par('din') <= needed)) {
      stop(paste0(
        'plotZmap() cannot draw a decorated z-map on a ', size, 'px image at ',
        'pointsize ', pointsize, ': the margins, labels and colour scale need ',
        'at least ', ceiling(max(needed) * size / par('din')[1]) + 1, 'px. ',
        'Use a larger size, a smaller pointsize, or decoration = FALSE, which ',
        'draws the z-map alone and works at any size.'))
    }

    # Widen the right margin so the colour bar has somewhere to go. raster::plot()
    # reserved that space itself; drawn into the default margins the bar lands on
    # top of the map's right-hand edge.
    #
    # Restored immediately via on.exit rather than at the end of the function:
    # CRAN's review of 1.2.1 asked for exactly that pattern inside functions.
    oldpar <- par(mar = zmapMargins)
    on.exit(par(oldpar), add = TRUE, after = FALSE)

    # A col in ... replaces the palettes below rather than arriving alongside
    # them, which errors before anything is drawn. ?plotZmap has always
    # documented col as the way to change the palette, and it errored in every
    # released version.
    dots <- list(...)
    user_col <- dots$col
    dots$col <- NULL
    # add is structural rather than a default: the overlay has to be added to
    # the map below it, and the first layer has to start one.
    dots$add <- NULL
    base_col <- if (is.null(user_col)) viridis::viridis(100) else user_col
    # Absent a caller palette the overlay renders in the default one rather than
    # the viridis set above, as raster::plot() did -- it was passed no col
    # either. Changing that would restyle every existing figure.
    overlay_col <- if (is.null(user_col)) zmapDefaultPalette() else user_col

    # Same merge as in drawZmapLayer(): these are defaults, and a caller who
    # names one of them replaces it instead of colliding with it.
    layer_args <- list(axes = FALSE, asp = 1,
                       main = paste0('Z-map of ', filename),
                       xlab = paste0('sigma = ', sigma, ', threshold = ', threshold))
    layer_args[names(dots)] <- dots

    # Initial (dummy) plot; sets up plot with initial dimensions + scale, title, label
    do.call(drawZmapLayer, c(list(zmap, col = base_col), layer_args))
    # Add bgimage (if specified) and superimpose Z-map on top of it.
    if (!(identical(bgimage, ''))) {
      rasterImage(bgimage, 0, 0, 1, 1)
      overlay_args <- dots
      overlay_args$add <- TRUE
      do.call(drawZmapLayer, c(list(zmap, col = overlay_col), overlay_args))
    }
    # If no bgimage was specified, draw a boundary box around the Z-map
    if (identical(bgimage, '')) {
      box <- matrix(NA, nrow(zmap) + 1, ncol(zmap) + 1)
      box[c(1, nrow(zmap) + 1), ] <- 0
      box[, c(1, ncol(zmap) + 1)] <- 0
      rasterImage(box, 0, 0, 1, 1)
    }
    # The colour bar goes last, and has to: it draws in its own figure region,
    # which leaves the device's user coordinates belonging to the bar rather
    # than to the map. Anything drawn into 0..1 afterwards -- the background,
    # the overlay, the boundary box -- would land in the bar's coordinate space
    # instead. Drawn between the map and the box, the box came out as two thin
    # lines across the middle of the figure.
    drawZmapLegend(zmap, col = if (identical(bgimage, '')) base_col else overlay_col,
                   zlim = dots$zlim, breaks = dots$breaks)
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
    drawZmapLayer(zmap, col = zmapDefaultPalette(), add = TRUE)
  }

  # Reached only when something was drawn, so the handler below keeps the file.
  drawn <- TRUE
  # The device is closed by the on.exit() handler registered after png().
}
