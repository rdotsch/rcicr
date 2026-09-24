#' Generates classification image
#'
#' Generate a classification image for any reverse correlation task.
#'
#' This function returns the classification image (CI) and, by default, saves it as a PNG. How
#' the CI is scaled for display decides what the image looks like and whether two images can be
#' compared. The default, \code{'independent'}, picks the lowest scaling constant that avoids
#' clipping this particular CI (see \code{'constant'} below for the formula). Each CI therefore
#' gets its own constant, and CIs with different noise ranges cannot be compared by eye.
#'
#' \code{'matched'} scaling matches the range of the CI's pixel intensities to the range of the
#' base image's. This is nonlinear and depends on the ranges of both. It also shifts the zero point
#' of the noise: a pixel that would not have changed the base image before scaling may change it
#' afterwards, and the other way round. Use it as a quick look at how the noise affects the base
#' image, not for reporting.
#'
#' \code{'constant'} scaling does not depend on the base image or the noise range, but the
#' constant is yours to choose, with the \code{scaling_constant} argument. The noise is scaled as
#' \code{scaled <- (ci + constant) / (2 * constant)}. Pixel intensities must lie between 0 and 1;
#' if the scaled noise falls outside that range you get a warning and should pick a higher
#' constant. The higher the constant, the fainter the noise in the resulting image. Use the same
#' constant for every classification image you want to compare.
#'
#' For several classification images, the lowest constant that works for all of them is a good
#' choice. \code{\link{autoscale}} finds it for you.
#'
#' @section Repeated presentations of the same stimulus:
#' When \code{participants} is \code{NA} (the default), repeated presentations of the same
#' stimulus are averaged before the CI is built, so each unique stimulus gets equal weight however
#' often it was shown. If every stimulus was shown equally often, this is the same as weighting
#' each trial equally. If not, it changes what the CI estimates: a stimulus shown three times counts
#' the same as one shown once, not three times as much.
#'
#' So in an unbalanced design, where a participant saw some stimuli more often than others (an
#' adaptive procedure, a crashed session, or by design), the CI is the average response per unique
#' stimulus, not per trial. The two can differ a lot: on an 8-trial set with counts 4/2/1/1 they
#' correlate at 0.77.
#'
#' \code{\link{computeCumulativeCICorrelation}} does \emph{not} average repeats; it weights each
#' trial equally. With unequal counts, the final CI it computes itself therefore differs from the
#' one this function returns. To compare against the CI you will report, pass this function's
#' output as its \code{targetci}.
#'
#' @export
#' @import png
#' @import parallel
#' @import doSNOW
#' @import foreach
#' @importFrom stats aggregate t.test qnorm
#' @importFrom spatstat.geom as.im
#' @importFrom spatstat.explore blur
#' @param stimuli Numeric vector of stimulus numbers, one per response and in the same order. Each must be a positive whole number no larger than the number of trials saved for the selected base image. Numbers may repeat and need not be consecutive. Factors, characters and logicals are rejected: if your data hold stimulus labels, check them against the generated stimulus filenames before converting them to numbers.
#' @param responses Vector of responses in the same order as \code{stimuli}: 1 where the original stimulus was chosen, -1 where the inverted one was.
#' @param baseimage String naming the base image: not its file name, but its key in the \code{base_face_files} list passed to \code{\link{generateStimuli2IFC}}.
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus.
#' @param save_as_png Boolean: also save the CI as a PNG image.
#' @param participants Optional vector with one participant ID per trial, the same length as \code{stimuli} and \code{responses}. When given, the CI is computed in two steps: one CI per participant, then their average. When missing or all \code{NA}, one CI is computed from all trials together. Some, but not all, \code{NA} is an error: give every trial an ID, or remove the trials without one.
#' @param save_individual_cis Boolean: when \code{participants} is given, also save each participant's CI as a PNG image.
#' @param targetpath Directory to save PNGs to. Required when \code{save_as_png = TRUE} or \code{save_individual_cis = TRUE}; there is no default. The directory is created if it does not exist; to just try the function out, use \code{tempdir()}.
#' @param filename Optional file name for the PNG image.
#' @param antiCI Boolean: compute the anti-CI, the classification image with its sign flipped, instead of the CI.
#' @param scaling Scaling method: \code{none}, \code{constant}, \code{matched} or \code{independent} (default). When both individual and group CIs are computed, this applies to the group CI.
#' @param scaling_constant Scaling constant for the noise, used only when \code{scaling = 'constant'}. When both individual and group CIs are computed, this applies to the group CI.
#' @param individual_scaling Scaling method for the individual CIs: \code{none}, \code{constant} or \code{independent} (default).
#' @param individual_scaling_constant Scaling constant for the individual CIs, used only when \code{individual_scaling = 'constant'}.
#' @param mask Optional mask for the CI: a 2D matrix (0 = masked, 1 = kept) or the path to a greyscale PNG image (black = masked, white = kept). Default: \code{NA}, no mask. Documentation up to and including 1.1.0 described the matrix the wrong way round (1 = masked); the code has always masked where the matrix is 0.
#' @param zmap Boolean: also create a z-map (default: \code{FALSE}).
#' @param zmapmethod Method for the z-map: \code{quick} (default) or \code{t.test}.
#' @param zmapdecoration Boolean: draw the z-map with margins, a caption (sigma, threshold) and a scale (default: \code{TRUE}).
#' @param zmappointsize Text size of the z-map decoration, in points (default: 12). Passed to \code{\link{plotZmap}}, which makes the z-map image \code{img_size} pixels wide. The decoration needs roughly \code{12.3 * zmappointsize} pixels on a 72 ppi device and \code{16.4 * zmappointsize} on a 96 ppi one. At the default, a stimulus set smaller than about 160 to 200 pixels is too small for it, and \code{generateCI()} stops with an error naming the minimum for your device. Lower this value to fit the decoration on a small z-map, or set \code{zmapdecoration = FALSE}.
#' @param sigma Amount of smoothing applied when creating the z-map (default: 3).
#' @param threshold Threshold z-score (default: 3). Z-scores below it are not drawn on the z-map.
#' @param zmaptargetpath Directory to save z-map PNGs to. Required when \code{zmap = TRUE}; there is no default. The directory is created if it does not exist; to just try the function out, use \code{tempdir()}.
#' @param n_cores Number of CPU cores used to create the z-map (default: \code{detectCores() - 1}; 2 under \code{R CMD check}, per CRAN policy).
#' @return List of pixel matrices: the raw classification noise (\code{ci}), the scaled noise (\code{scaled}), the base image (\code{base}) and the two combined (\code{combined}).
#' @examples
#' # a synthetic square grayscale image stands in for a real base face photo
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#'
#' stimulus_path <- tempdir()
#' generateStimuli2IFC(
#'   base_face_files = list(face = base_face),
#'   n_trials = 6,
#'   img_size = 32,
#'   stimulus_path = stimulus_path,
#'   seed = 1,
#'   ncores = 1,
#'   nscales = 1,
#'   save_as_png = FALSE
#' )
#' rdata_file <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)[1]
#'
#' responses <- sample(c(1, -1), 6, replace = TRUE)
#' ci <- generateCI(
#'   stimuli = 1:6, responses = responses, baseimage = "face",
#'   rdata = rdata_file, save_as_png = FALSE
#' )
# Main function -----------------------------------------------------------
generateCI <- function(stimuli, responses, baseimage, rdata, participants = NA,
                       save_individual_cis = FALSE, save_as_png = TRUE, filename = '',
                       targetpath, antiCI = FALSE, scaling = 'independent',
                       scaling_constant = 0.1, individual_scaling = 'independent',
                       individual_scaling_constant = 0.1, zmap = FALSE,
                       zmapmethod = 'quick', zmapdecoration = TRUE, sigma = 3,
                       threshold = 3, zmaptargetpath,
                       n_cores = default_ncores(), mask = NA,
                       # Appended, never inserted: a new formal in the middle
                       # rebinds every positional argument after it in scripts
                       # that already exist. Here it would have taken sigma's
                       # value and left the z-map blurred at the default.
                       zmappointsize = 12) {

  # Preprocessing -----------------------------------------------------------

  # targetpath and zmaptargetpath are required, not defaulted: a default path
  # writes to the user's filespace uninvited, which CRAN policy does not allow.
  # Both checks must run before the load(rdata) below, which assigns into this
  # frame and can replace an argument with a value from the file.
  if ((save_as_png || save_individual_cis) && missing(targetpath)) {
    stop(paste0('save_as_png or save_individual_cis is TRUE but no targetpath ',
      'was given. Supply targetpath = <a directory> to say where the ',
      'PNGs should go, or set both to FALSE to compute the ',
      'classification image without writing it. Use tempdir() if you ',
      'only want to try the function out.'
    ))
  }
  if (isTRUE(zmap) && missing(zmaptargetpath)) {
    stop(paste0('zmap is TRUE but no zmaptargetpath was given. Supply ',
      'zmaptargetpath = <a directory> to say where the z-map PNG ',
      'should go. Use tempdir() if you only want to try the function ',
      'out.'
    ))
  }

  # Bind targetpath even when it was not supplied. foreach::getexports() scans
  # the %dopar% body for free variables and get()s each one, including the
  # targetpath inside the save_individual_cis branch below - a branch that
  # cannot run when targetpath is absent. Leaving it unbound aborted every
  # participant-CI call with n_cores > 1, which is the default (#235).
  #
  # Must stay below the missing() checks above, which stop being reliable for an
  # argument once it has been assigned to.
  targetpath <- if (missing(targetpath)) NULL else targetpath

  trials <- coerceTrialVectors(stimuli, responses, participants)
  stimuli <- trials$stimuli
  responses <- trials$responses
  participants <- trials$participants

  # Loaded in a frame of its own, so no field of the .Rdata is ever in scope
  # beside this function's arguments -- see loadStimulusParams() in R/rdata.R for
  # what that prevents.
  loaded <- loadStimulusParams(rdata)
  p <- loaded$p
  base_faces <- loaded$base_faces
  stimuli_params <- loaded$stimuli_params
  img_size <- loaded$img_size

  base <- selectBaseImage(base_faces, baseimage)

  validateStimulusIds(stimuli, nrow(stimuli_params[[baseimage]]))

  if (all(is.na(participants))) {
    aggregated <- aggregateResponses(stimuli, responses)
    stimuli <- aggregated$stimuli
    responses <- aggregated$responses
  }

  params <- selectStimulusParams(stimuli_params, baseimage, stimuli)

  # Generate CI(s) ----------------------------------------------------------

  # Invert parameters if antiCI is to be generated
  if (antiCI == TRUE) {
    params <- -params
  }

  # If "participants" argument is not given, compute one CI based on all data
  if (all(is.na(participants))) {
    ci <- generateCINoise(params, responses, p)
    # Bound so the t.test z-map below can switch on it. NULL is what "no
    # per-participant CIs to reuse, build the stack from trials" looks like.
    pid.cis <- NULL # nolint: object_name_linter.
    # If it is given, create a CI for each participant and a group CI by
    # averaging across participants
  } else {

    participant_cis <- computeParticipantCIs(params, responses, participants, p,
      base, baseimage, img_size, mask, n_cores, save_individual_cis, targetpath,
      individual_scaling, individual_scaling_constant, antiCI
    )
    ci <- participant_cis$ci
    pid.cis <- participant_cis$pid_cis # nolint: object_name_linter.
  }

  # Check if a mask has been set. If so, apply it to the CI
  if (hasMask(mask)) {
    ci <- applyMask(ci, mask, img_size)
  }

  # Apply scaling
  scaled <- applyScaling(base, ci, scaling, scaling_constant)

  # Combine with base image
  combined <- combine(scaled, base)

  # Save CI as PNG
  if (save_as_png) {
    saveToImage(baseimage, combined, targetpath, filename, antiCI)
  }

  # Rename zmap to zmapbool so we can use zmap for the actual zmap
  zmapbool <- zmap
  if (zmapbool) {

    if (zmapmethod == 'quick') {
      zmap <- computeZmapQuick(ci, sigma, threshold, img_size)
    }

    if (zmapmethod == 't.test') {
      zmap <- computeZmapTTest(ci, params, responses, p, pid.cis, img_size,
        n_cores
      )
    }

    # Pass zmap object to plotZmap for plotting. targetpath was previously not
    # forwarded, so the documented zmaptargetpath argument was silently ignored
    # and every z-map went to plotZmap()'s own default ('zmaps', relative to the
    # working directory) no matter what the caller asked for.
    plotZmap(zmap = zmap, bgimage = combined, filename = baseimage,
      sigma = sigma, threshold = threshold, size = img_size,
      decoration = zmapdecoration, pointsize = zmappointsize,
      targetpath = zmaptargetpath
    )
  }

  # Return data
  if (zmapbool) {
    return(list(ci = ci, scaled = scaled, base = base, combined = combined, zmap = zmap))
  } else {
    return(list(ci = ci, scaled = scaled, base = base, combined = combined))
  }
}

# Functions ---------------------------------------------------------------

# Apply masking to a CI
# Has the user actually supplied a mask?
# The `mask` argument defaults to NA, so a plain `!is.na(mask)` test returns a
# whole matrix when a mask *is* supplied - which R >= 4.2 rejects outright with
# "the condition has length > 1". This collapses the test to a single logical.
# The sentinel is specifically an atomic scalar with no dim: matrix(NA, 1, 1)
# and list(NA) are also length 1 and is.na(), and both used to be read as "no
# mask" and discarded in silence rather than reaching the validation that
# rejects them.
# Input: mask (NA, NULL, a string path, or a matrix)
# Output: TRUE if a mask was supplied
hasMask <- function(mask) {
  !is.null(mask) &&
    !(is.atomic(mask) && length(mask) == 1L && is.na(mask) && is.null(dim(mask)))
}

# Input: CI (or z-map), mask (either a string or a matrix), expected target
# size, and a label for the size-mismatch message (what the mask is being
# checked against -- generateCI() and plotZmap() share this helper, and "the
# stimuli" is only accurate for the former)
# Output: masked matrix (input matrix, but masked pixels are NA)
applyMask <- function(ci, mask, img_size = nrow(ci), context = 'stimuli') {
  # If mask argument is a string, treat it as a path to a bitmap and try to read
  # it into a matrix. If it is a matrix, use it. Else, throw an error
  if (typeof(mask) == 'character') {
    mask_matrix <- png::readPNG(mask)

    # Check if the PNG uses a greyscale color palette
    if (length(dim(mask_matrix)) != 2) {
      # A trailing channel is alpha -- not colour information, and never
      # compared -- whenever the total channel count is even: 2 (greyscale +
      # alpha) or 4 (RGBA). Every other channel must agree with channel 1 for
      # the image to be greyscale-as-RGB(A).
      # Thanks https://stackoverflow.com/a/30850654
      n <- dim(mask_matrix)[3]
      n_color <- if (n %in% c(2, 4)) n - 1L else n
      if (n_color > 1 && !all(sapply(2:n_color, function(i) {
        identical(mask_matrix[, , i], mask_matrix[, , 1])
      }))) {
        # Only error if the colour channels genuinely differ. This stop() used
        # to run unconditionally, so even a convertible greyscale-as-RGB PNG
        # failed.
        stop(paste0('This PNG is not encoded with a greyscale color palette and ',
          'could not be converted to this encoding either. In other ',
          'words, this is not a greyscale image.'
        ))
      }
      # `[, , 1]` alone would also drop a singleton *spatial* dimension, leaving
      # a dim-less vector -- and `all(NULL == img_size)` is vacuously TRUE, so a
      # 1-by-8 mask would pass the size check below for a 2-by-4 target and then
      # be applied by linear indexing. plotZmap()'s previous inline code checked
      # the PNG's spatial dimensions before dropping channels and so rejected it.
      spatial <- dim(mask_matrix)[1:2]
      mask_matrix <- matrix(mask_matrix[, , 1],
        nrow = spatial[1], ncol = spatial[2]
      )
    }
  } else if (is.matrix(mask) && length(dim(mask)) == 2) {
    mask_matrix <- mask
  } else {
    stop('The mask argument is neither a string nor a matrix!')
  }

  # Check if mask is of the same size as the target (i.e. img_size). This used
  # to compare against a hardcoded 512, so masks failed for every other
  # stimulus size, and reported img_size - which is not in scope here - in the
  # error message. img_size[1] / img_size[length(img_size)] read correctly
  # whether img_size is a scalar (generateCI()'s calls) or the length-2
  # c(rows, cols) plotZmap() passes for a possibly-rectangular zmap.
  if (!all(dim(mask_matrix) == img_size)) {
    stop(paste0('Mask is not of the same dimensions as the ', context, '! ',
      '(', context, ' dimensions: ', img_size[1], ' x ',
      img_size[length(img_size)],
      '; mask dimensions: ', dim(mask_matrix)[1],
      ' by ', dim(mask_matrix)[2], ').'
    ))
  }

  # Check if the mask is binary
  if (length(mask_matrix) != sum(mask_matrix %in% c(0, 1))) {
    stop('This mask contains values other than 0 or 1!')
  }

  # Convert mask to boolean matrix (black == 0 == masked)
  mask <- mask_matrix == 0

  # Apply the mask to the CI. This replaces all the masked pixels with NA
  ci[mask] <- NA

  # Return the masked CI
  return(ci)
}

# Apply scaling to a CI
# Input: base image, CI, scaling method, constant
# Output: scaled CI
applyScaling <- function(base, ci, scaling, constant) {
  # No scaling
  if (scaling == 'none') {
    scaled <- ci
    # Scaling with a constant scaling factor
  } else if (scaling == 'constant') {
    scaled <- (ci + constant) / (2 * constant)
    if (max(scaled[!is.na(scaled)]) > 1.0 || min(scaled[!is.na(scaled)]) < 0) {
      warning(paste0('Chosen constant value for constant scaling made noise ',
        'of classification image exceed possible intensity range ',
        'of pixels (<0 or >1). Choose a lower value, or clipping ',
        'will occur.'
      ))
    }
    # Scaling using 'matched' method
  } else if (scaling == 'matched') {
    values <- ci[!is.na(ci)]
    if (length(values) > 0 && max(values) == min(values)) {
      # No CI range to map onto the base image's, so the midpoint of that range
      # is the neutral answer. Deliberately not a literal 0.5: this method
      # renders into the base image's own intensity range, and 0.5 can sit
      # outside it -- a base spanning [0, 0.3] would show a no-signal CI
      # brighter than any pixel in the base.
      warnDegenerateScaling(scaling)
      scaled <- neutralScaling(ci, min(base) + (max(base) - min(base)) / 2)
    } else {
      scaled <- min(base) +
        ((max(base) - min(base)) * (ci - min(values)) /
           (max(values) - min(values)))
    }
    # Scaling with maximum scaling factor for the given CI
  } else if (scaling == "independent") {

    # Determine the lowest possible scaling factor constant
    if (abs(range(ci[!is.na(ci)])[1]) > abs(range(ci[!is.na(ci)])[2])) {
      constant <- abs(range(ci[!is.na(ci)])[1])
    } else {
      constant <- abs(range(ci[!is.na(ci)])[2])
    }

    if (isTRUE(constant == 0)) {
      warnDegenerateScaling(scaling)
      scaled <- neutralScaling(ci, 0.5)
    } else {
      scaled <- (ci + constant) / (2 * constant)
    }
    # Print warning when scaling method name is not recognized
  } else {
    warning(paste0('Scaling method \'', scaling, '\' not found. Using none.'))
    scaled <- ci
  }

  # Return the scaled CI
  return(scaled)
}

# Render a classification image that has nothing to scale.
#
# Responses that cancel exactly give an all-zero CI: a result, not a failure.
# The range-based methods would divide by a zero range and return NaN for every
# pixel, so they fill the image with a neutral value instead -- what 'constant'
# scaling already returns for the same input, and what autoscale() already
# returns for a zero CI sitting beside one with signal. Masked pixels keep
# their NA.
neutralScaling <- function(ci, value) {
  scaled <- ci
  scaled[!is.na(ci)] <- value
  return(scaled)
}

warnDegenerateScaling <- function(scaling) {
  msg <- paste0('This classification image has no range to scale: every ',
    'pixel holds the same value, which is what exactly cancelling ',
    'responses produce. Under \'', scaling, '\' scaling it is ',
    'rendered as a uniform neutral image rather than as NaN. The ',
    'unscaled CI in $ci is unaffected.'
  )
  warning(msg, call. = FALSE)
}

# Combine scaled CI with base image
# Input: scaled CI, base image
# Output: CI with base image
combine <- function(scaled, base) {
  return((scaled + base) / 2)
}

# Save a CI to an image file
# Input: base image name, combined CI, target path, filename, CI/antiCI boolean
# Output: nothing (image is saved to file)
saveToImage <- function(baseimage, combined, targetpath, filename, antiCI) {
  # If no filename is specified, default to name of base image
  if (filename == '') {
    filename <- paste0(baseimage)
  }

  # Add ci/antici prefix to filename
  if (antiCI) {
    filename <- paste0('antici_', filename)
  } else {
    filename <- paste0('ci_', filename)
  }

  # Add extension to filename
  filename <- paste0(filename, '.png')

  # Create output directory
  dir.create(targetpath, recursive = TRUE, showWarnings = FALSE)

  # Write CI to image file
  png::writePNG(combined, paste0(targetpath, '/', filename))
}
