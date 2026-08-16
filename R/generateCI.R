#' Generates classification image
#'
#' Generate classification image for any reverse correlation task.
#'
#' This function saves the classification image as PNG to a folder and returns the CI. Your choice of scaling
#' matters. The default, \code{'independent'}, picks the lowest scaling constant that avoids clipping this
#' particular classification image (see \code{'constant'} scaling below for the formula), so it is not
#' comparable across classification images with different noise ranges.
#'
#' \code{'matched'} scaling will match the range of the intensity of the pixels to
#' the range of the base image pixels. This scaling is nonlinear and depends on the range of both base image
#' and noise pattern. It is truly suboptimal, because it shifts the 0 point of the noise (that is, pixels that would
#' not have changed the base image at all before scaling may change the base image after scaling and vice versa). It is
#' however the quick and dirty way to see how the CI noise affects the base image.
#'
#' For more control, use \code{'constant'} scaling, where the scaling is independent of
#' the base image and noise range, but where the choice of constant is arbitrary (provided by the user with
#' the \code{constant} parameter). The noise is then scale as follows: \code{scaled <- (ci + constant) / (2*constant)}.
#' Note that pixels can take intensity values between 0 and 1. If your scaled noise exceeds those values,
#' a warning will be given. You should pick a higher constant (but do so consistently for different classification images
#' that you want to compare). The higher the constant, the less visible the noise will be in the resulting image.
#'
#' When creating multiple classification images a good strategy is to find the lowest constant that works for all
#' classification images. This can be automatized using the \code{autoscale} function.
#'
#' @section Repeated presentations of the same stimulus:
#' When \code{participants} is \code{NA} (the default), repeated presentations of the same
#' stimulus are collapsed before building the CI: each unique stimulus gets equal weight,
#' regardless of how many times it was presented. Where every stimulus was presented the same
#' number of times, this is equivalent to weighting each trial equally and changes nothing. Where
#' repeat counts differ, it changes the estimand: a stimulus presented three times counts the same
#' as one presented once, rather than three times as much.
#'
#' This is worth knowing for unbalanced designs. If a participant saw some stimuli more often than
#' others -- because of an adaptive procedure, a crashed session, or a design choice -- the CI
#' reflects the average response per unique stimulus, not per trial. The difference can be
#' substantial: on an 8-trial set with counts 4/2/1/1 the two weightings correlate at 0.77.
#'
#' \code{\link{computeCumulativeCICorrelation}} does \emph{not} aggregate and weights each trial
#' equally, so its self-computed final CI diverges from the one this function returns under
#' unequal counts. Pass this function's output as \code{targetci} to compare against the CI you
#' will actually report.
#'
#' @export
#' @import png
#' @import parallel
#' @import doSNOW
#' @import foreach
#' @importFrom stats aggregate t.test qnorm
#' @importFrom spatstat.geom as.im
#' @importFrom spatstat.explore blur
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli.
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param save_as_png Optional boolean stating whether to additionally save the CI as PNG image.
#' @param participants Optional vector specifying participant IDs. If specified, will compute the requested CIs in two steps: step 1, compute CI for each participant. Step 2, compute final CI by averaging participant CIs. If unspecified, the function defaults to averaging all data in the stimuli and responses vector.
#' @param save_individual_cis Optional boolean specifying whether individual CIs should be save as PNG images when the \code{participants} parameter is used.
#' @param targetpath String specifying the directory to save PNGs to. Required when \code{save_as_png = TRUE} or \code{save_individual_cis = TRUE}; there is no default path. It is created if it does not exist. Use \code{tempdir()} if you only want to try the function out.
#' @param filename Optional string to specify a file name for the PNG image.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, or \code{independent} (default). This scaling applies to the group-level CIs if both individual-level and group-level CIs are being generated.
#' @param scaling_constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}). This scaling applies to the group-level CIs if both individual-level and group-level CIs are being generated.
#' @param individual_scaling Optional string specifying scaling method for individual CIs: \code{none}, \code{constant}, \code{independent} (default).
#' @param individual_scaling_constant Optional number specifying the value used as constant scaling factor for the noise of individual CIs (only works for \code{individual_scaling='constant'}).
#' @param mask Optional 2D matrix that defines the mask to be applied to the CI (0 = masked, 1 = unmasked). May also be a string specifying the path to a grayscale PNG image (black = masked, white = unmasked). Default: NA. Note the matrix convention was documented the wrong way round (as 1 = masked) up to and including 1.1.0; the code has always masked where the matrix is 0, matching the PNG form, and that is what is described here.
#' @param zmap Boolean specifying whether a z-map should be created (default: FALSE).
#' @param zmapmethod String specifying the method to create the z-map. Can be: \code{quick} (default), \code{t.test}.
#' @param zmapdecoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param zmappointsize Integer specifying the text size of the Z-map decoration, in points (default: 12). Passed to \code{plotZmap()}, which sizes the Z-map image to \code{img_size}. The decoration needs roughly \code{12.3 * zmappointsize} pixels on a 72 ppi device and \code{16.4 * zmappointsize} on a 96 ppi one, so a stimulus set below about 160-200px cannot carry it at the default and \code{generateCI()} stops with an error naming the minimum for the device in use. Lower this to fit the decoration onto a small Z-map, or set \code{zmapdecoration = FALSE}.
#' @param sigma Integer specifying the amount of smoothing to apply when generating the z-maps (default: 3).
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param zmaptargetpath String specifying the directory to save z-map PNGs to. Required when \code{zmap = TRUE}; there is no default path. It is created if it does not exist. Use \code{tempdir()} if you only want to try the function out.
#' @param n_cores Optional integer specifying the number of CPU cores to use to generate the z-map (default: \code{detectCores()-1}; 2 under \code{R CMD check}, per CRAN policy).
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined.
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
  # argument once it has been assigned to, and above captureArgs(), which skips
  # required-and-absent arguments and so would otherwise leave this one exposed
  # to the load() below.
  targetpath <- if (missing(targetpath)) NULL else targetpath

  # Must stay above captureArgs() below, which snapshots what these resolve to.
  trials <- coerceTrialVectors(stimuli, responses, participants)
  stimuli <- trials$stimuli
  responses <- trials$responses
  participants <- trials$participants

  # load() assigns straight into this function's frame, so any object stored in
  # the .Rdata file silently overwrites an argument of the same name (the same
  # hazard handled in generateReferenceDistribution2IFC()). `sigma` is the live
  # case: since 1.1.0 the file carries the *noise* sigma - 25 by default - which
  # replaced this function's z-map blur sigma of 3, so every z-map made from a
  # 1.1.0-generated stimulus set was smoothed with the wrong constant and the
  # sigma the caller passed was ignored. Keep private copies of every argument
  # and restore them after loading, so a field added to the .Rdata later cannot
  # quietly capture another one.
  .args <- captureArgs(environment())

  # Load parameter file (created when generating stimuli)
  load(rdata)

  list2env(.args, envir = environment())

  # Check whether critical variables have been loaded
  if (!exists('s', envir = environment(), inherits = FALSE) &&
        !exists('p', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata did not contain s or p variable.', rdataWriterNote(environment()))
  }

  if (!exists('base_faces', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata did not contain base_faces variable.', rdataWriterNote(environment()))
  }

  if (!exists('stimuli_params', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata did not contain stimuli_params variable.', rdataWriterNote(environment()))
  }

  if (!exists('img_size', envir = environment(), inherits = FALSE)) {
    stop('File specified in rdata did not contain img_size variable.', rdataWriterNote(environment()))
  }

  # Convert s to p (if rdata file originates from pre-0.3.3)
  if (exists('s', envir = environment(), inherits = FALSE)) {
    p <- list(patches = s$sinusoids, patchIdx = s$sinIdx, noise_type = 'sinusoid')
    rm(s)
  }

  base <- selectBaseImage(base_faces, baseimage)

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
    # If it is given, create a CI for each participant and a group CI by
    # averaging across participants
  } else {

    # First generate a CI for each participant, then average across participants
    pids <- as.numeric(factor(participants))
    npids <- length(unique(pids))

    # Initialize progress bar
    pb <- txtProgressBar(min = 1, max = npids, style = 3)

    # Create cluster for parallel processing
    cl <- startBackend(n_cores)
    if (!is.null(cl)) {
      on.exit(stopClusterSafely(cl), add = TRUE)
    }

    # For each weighted stimulus, construct the noise pattern
    pid.cis <- foreach::foreach(obs = 1:npids, # nolint: object_name_linter.
      .combine = 'c',
      .packages = 'rcicr',
      .options.snow = progressOption(pb, cl)
    ) %dopar% {

      # Serial path only; in parallel .options.snow ticks the bar in the parent.
      if (is.null(cl)) setTxtProgressBar(pb, obs)

      # Select only the observations of the current participant
      pid.rows <- pids == obs # nolint: object_name_linter.

      # Construct the noise pattern
      ci <- generateCINoise(params[pid.rows, ], responses[pid.rows], p)

      # Check if individual CIs should be saved. If so, generate and save them
      if (save_individual_cis) {
        if (hasMask(mask)) {
          individual_ci <- applyMask(ci, mask, img_size)
        } else {
          individual_ci <- ci
        }
        scaled <- applyScaling(base, individual_ci, individual_scaling,
          individual_scaling_constant
        )
        combined <- combine(scaled, base)
        saveToImage(baseimage, combined, paste0(targetpath, '/individual_cis'),
          unique(participants)[obs], antiCI
        )
      }

      # Return the CI
      return(ci)
    }
    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }
    cl <- NULL
    dim(pid.cis) <- c(img_size, img_size, npids) # nolint: object_name_linter.

    # Average across participants for final CI and return to original variance
    ci <- apply(pid.cis, c(1, 2), mean) #* sqrt(npids)
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
      # Blur CI
      zmap <- as.matrix(blur(as.im(ci), sigma = sigma))

      # Create z-map
      zmap <- matrix(scale(as.vector(zmap)), img_size, img_size)

      # Apply threshold
      zmap[zmap > -threshold & zmap < threshold] <- NA
    }

    if (zmapmethod == 't.test') {

      # Compute one CI in one single step based on all data
      if (all(is.na(participants))) {
        # Weigh the stimulus parameters of each trial using the given responses
        weightedparameters <- params * responses

        # Get number of observations
        n_observations <- length(responses)

        # Initialize progress bar
        pb <- txtProgressBar(min = 1, max = n_observations, style = 3)

        # Create cluster for parallel processing
        cl <- startBackend(n_cores)
        if (!is.null(cl)) {
          on.exit(stopClusterSafely(cl), add = TRUE)
        }

        # For each weighted stimulus, construct the complementary noise pattern
        noiseimages <- foreach::foreach(obs = 1:n_observations, .combine = 'c',
          .packages = 'rcicr',
          .options.snow = progressOption(pb, cl)
        ) %dopar% {
          noiseimage <- generateNoiseImage(weightedparameters[obs, ], p)
          if (is.null(cl)) setTxtProgressBar(pb, obs)
          return(noiseimage)
        }
        if (!is.null(cl)) {
          parallel::stopCluster(cl)
        }
        cl <- NULL
        dim(noiseimages) <- c(img_size, img_size, n_observations)

      } else {
        noiseimages <- pid.cis
      }

      # Get p value for each pixel
      pmap <- apply(noiseimages, 1:2, function(x) unlist(t.test(x)['p.value']))

      # Create Z-map
      zmap <- sign(ci) * abs(qnorm(pmap / 2))
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
    scaled <- min(base) +
      ((max(base) - min(base)) * (ci - min(ci[!is.na(ci)])) /
         (max(ci[!is.na(ci)]) - min(ci[!is.na(ci)])))
    # Scaling with maximum scaling factor for the given CI
  } else if (scaling == "independent") {

    # Determine the lowest possible scaling factor constant
    if (abs(range(ci[!is.na(ci)])[1]) > abs(range(ci[!is.na(ci)])[2])) {
      constant <- abs(range(ci[!is.na(ci)])[1])
    } else {
      constant <- abs(range(ci[!is.na(ci)])[2])
    }

    scaled <- (ci + constant) / (2 * constant)
    # Print warning when scaling method name is not recognized
  } else {
    warning(paste0('Scaling method \'', scaling, '\' not found. Using none.'))
    scaled <- ci
  }

  # Return the scaled CI
  return(scaled)
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
