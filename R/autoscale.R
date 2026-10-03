#' Determines optimal scaling constant for a list of CIs
#'
#' @export
#' @import png
#' @param cis List of classification images, each a list of pixel matrices holding at least the noise (\code{$ci}), plus the base image (\code{$base}) if PNGs are to be written.
#' @param save_as_pngs Boolean: combine each autoscaled noise pattern with its base image and save it as a PNG, named after its key in \code{cis}. Every element then needs a name of its own; with \code{FALSE}, names are not needed.
#' @param targetpath Directory to save PNGs to. Required when \code{save_as_pngs = TRUE}; there is no default. The directory is created if it does not exist; to just try the function out, use \code{tempdir()}.
#' @return The input \code{cis} list, with each element's \code{$scaled} matrix replaced by its
#' autoscaled version. The scaling constant is printed to the console, and recorded in each
#' element's \code{scaling} attribute, which keeps the earlier record for \code{$combined};
#' \code{vignette("stored-data")} lists its fields.
#'
#' \strong{Look at \code{$scaled}, not \code{$combined}.} \code{$combined} is returned exactly
#' as it was passed in, on purpose, so that existing scripts that plot it keep producing the same
#' image. Only \code{$scaled} holds the autoscaled result. To show the autoscaled noise over the
#' base image, compute \code{(ci$scaled + ci$base) / 2}: that is exactly what
#' \code{save_as_pngs = TRUE} writes to disk.
#'
#' This catches people out after \code{\link{batchGenerateCI}} or
#' \code{\link{batchGenerateCI2IFC}}. Both scale with \code{'none'} before calling this
#' function, so their \code{$combined} overlays the \emph{unscaled} noise and looks almost blank,
#' while \code{$scaled} is the image you want.
#' @examples
#' cis <- list(
#'   participant1 = list(ci = matrix(runif(64, -0.2, 0.2), 8, 8), base = matrix(0.5, 8, 8)),
#'   participant2 = list(ci = matrix(runif(64, -0.3, 0.3), 8, 8), base = matrix(0.5, 8, 8))
#' )
#' scaled_cis <- autoscale(cis, save_as_pngs = FALSE)
autoscale <- function(cis, save_as_pngs = TRUE, targetpath) {

  # targetpath is required, not defaulted: a default path writes to the user's
  # filespace uninvited, which CRAN policy does not allow.
  if (save_as_pngs && missing(targetpath)) {
    stop(paste0('save_as_pngs is TRUE but no targetpath was given. Supply ',
                'targetpath = <a directory> to say where the PNGs should go, ',
                'or set save_as_pngs = FALSE to autoscale without writing ',
                'them. Use tempdir() if you only want to try the function out.'))
  }

  if (!is.list(cis) || !length(cis)) {
    stop('cis must be a non-empty list of classification images.', call. = FALSE)
  }
  labels <- ciLabels(cis)
  # Names become file names, so only a call that writes needs them unique.
  if (save_as_pngs) requirePngNames(names(cis))

  # Get range of each ci, by position: names may repeat or be absent (#375).
  #
  # na.rm is required: generateCI(mask = ...) leaves masked pixels NA (NEWS 1.2.0).
  ranges <- matrix(0, length(cis), 2)
  for (i in seq_along(cis)) {
    ci_values <- cis[[i]]$ci
    if (all(is.na(ci_values))) {
      stop(paste0('Classification image "', labels[i], '" is entirely NA, so there ',
                  'is no range to scale it against. If it was masked, check that ',
                  'the mask does not cover the whole image.'))
    }
    ranges[i, ] <- range(ci_values, na.rm = TRUE)
  }

  if (abs(min(ranges[, 1])) > max(ranges[, 2])) {
    constant <- abs(min(ranges[, 1]))
  }  else {
    constant <- max(ranges[, 2])
  }

  write(paste0("Using scaling factor constant:", constant), stdout())

  # A zero constant means every CI in the list is exactly zero, so dividing by
  # it would make all of them NaN. They render neutral instead, which is what a
  # zero CI already gets here whenever the list also holds one with signal.
  degenerate <- isTRUE(constant == 0)
  if (degenerate) {
    msg <- paste0('Every classification image in this list is exactly zero, ',
      'so there is no range to scale them against. They are rendered ',
      'as uniform neutral images rather than as NaN. The unscaled ',
      'CIs in $ci are unaffected.'
    )
    warning(msg, call. = FALSE)
  }

  for (i in seq_along(cis)) {
    cis[[i]]$scaled <- if (degenerate) {
      neutralScaling(cis[[i]]$ci, 0.5)
    } else {
      (cis[[i]]$ci + constant) / (2 * constant)
    }
    attr(cis[[i]], 'scaling') <- autoscaledRecord(attr(cis[[i]], 'scaling'), constant)

    # $combined is deliberately left as the caller supplied it: rewriting it
    # would change what existing scripts plot. $scaled is the autoscaled
    # result, and the PNG below is built from it.
    if (save_as_pngs) {
      ci <- (cis[[i]]$scaled + cis[[i]]$base) / 2

      dir.create(targetpath, recursive = TRUE, showWarnings = FALSE)

      png::writePNG(ci, paste0(targetpath, '/', names(cis)[i], '_autoscaled.png'))
    }

  }

  return(cis)
}

requirePngNames <- function(nms) {
  if (is.null(nms) || any(is.na(nms) | !nzchar(nms))) {
    stop('save_as_pngs = TRUE names each PNG after its element of cis, so every element ',
         'needs a name. Name them, or use save_as_pngs = FALSE.', call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop('save_as_pngs = TRUE names each PNG after its element of cis, and these names occur ',
         'more than once: ', paste(unique(nms[duplicated(nms)]), collapse = ', '),
         '. Give every element its own name, or use save_as_pngs = FALSE.', call. = FALSE)
  }
  invisible(NULL)
}

# The top level describes $scaled, which this function rewrites; `combined`
# keeps the record of the scaling $combined was made with, which it does not.
# NULL there means unknown: a CI from an older version, or built by hand.
autoscaledRecord <- function(previous, constant) {
  combined <- if (identical(previous$method, 'autoscale')) {
    previous$combined
  } else if (!is.null(previous)) {
    previous[c('method', 'constant')]
  }
  record <- list(method = 'autoscale', constant = constant, combined = combined)
  if (!is.null(previous$individual)) record$individual <- previous$individual
  record
}
