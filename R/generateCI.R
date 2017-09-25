#' Generates classification image
#'
#' Generate classification image for for any reverse correlation task.
#'
#' This function saves the classification image as PNG to a folder and returns the CI. Your choice of scaling
#' matters. The default is \code{'matched'}, and will match the range of the intensity of the pixels to
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
#' @export
#' @import png
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom stats aggregate t.test qnorm
#' @importFrom spatstat blur as.im
#' @importFrom grDevices png dev.off
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli.
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param save_as_png Optional boolean stating whether to additionally save the CI as PNG image.
#' @param participants Optional vector specifying participant IDs. If specified, will compute the requested CIs in two steps: step 1, compute CI for each participant. Step 2, compute final CI by averaging participant CIs. If unspecified, the function defaults to averaging all data in the stimuli and responses vector.
#' @param save_individual_cis Optional boolean specifying whether individual CIs should be save as PNG images when the \code{participants} parameter is used.
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @param filename Optional string to specify a file name for the PNG image.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, or \code{independent} (default). This scaling applies to the group-level CIs if both individual-level and group-level CIs are being generated.
#' @param scaling_constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}). This scaling applies to the group-level CIs if both individual-level and group-level CIs are being generated.
#' @param individual_scaling Optional string specifying scaling method for individual CIs: \code{none}, \code{constant}, \code{independent} (default).
#' @param individual_scaling_constant Optional number specifying the value used as constant scaling factor for the noise of individual CIs (only works for \code{individual_scaling='constant'}).
#' @param mask Optional 2D matrix that defines the mask to be applied to the CI (1 = masked, 0 = unmasked). May also be a string specifying the path to a grayscale PNG image (black = masked, white = unmasked). Default: NA.
#' @param zmap Boolean specifying whether a z-map should be created (default: FALSE).
#' @param zmapmethod String specifying the method to create the z-map. Can be: \code{quick} (default), \code{t.test}.
#' @param zmapdecoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param sigma Integer specifying the amount of smoothing to apply when generating the z-maps (default: 3).
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param zmaptargetpath Optional string specifying path to save z-map PNGs to (default: ./zmaps).
#' @param n_cores Optional integer specifying the number of CPU cores to use to generate the z-map (default: detectCores()).
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined.

# Main function -----------------------------------------------------------
generateCI <- function(stimuli, responses, baseimage, rdata, participants=NA,
                       save_individual_cis=FALSE, save_as_png=TRUE, filename='',
                       targetpath='./cis', antiCI=FALSE, scaling='independent',
                       scaling_constant=0.1, individual_scaling='independent',
                       individual_scaling_constant=0.1, zmap = F,
                       zmapmethod = 'quick', zmapdecoration = T, sigma = 3,
                       threshold = 3, zmaptargetpath = './zmaps',
                       n_cores = detectCores(), mask=NA) {

  # Preprocessing -----------------------------------------------------------

  # Load parameter file (created when generating stimuli)
  load(rdata)

  # Check whether critical variables have been loaded
  if (!exists('s', envir=environment(), inherits=FALSE) &
      !exists('p', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata did not contain s or p variable.')
  }

  if (!exists('base_faces', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata did not contain base_faces variable.')
  }

  if (!exists('stimuli_params', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata did not contain stimuli_params variable.')
  }

  if (!exists('img_size', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata did not contain img_size variable.')
  }


  # Convert s to p (if rdata file originates from pre-0.3.3)
  if (exists('s', envir=environment(), inherits=FALSE)) {
    p <- list(patches=s$sinusoids, patchIdx=s$sinIdx, noise_type='sinusoid')
    rm(s)
  }

  # Get base image
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    # If no base face with the given name is found in Rdata file, throw error
    stop(paste0('File specified in rdata argument did not contain any ',
                'reference to base image label: ', baseimage, ' (NOTE: file ',
                'contains references to the following base image label(s): ',
                paste(names(base_faces), collapse=', '), ')'))
  }

  if (all(is.na(participants))) {
    # Average responses for each presented stimulus (in case stimuli have been
    # presented multiple times, or group-wise classification images are being
    # calculated, in order to reduce memory and processing load)
    aggregated <- aggregate(responses, by=list(stimuli=stimuli), FUN=mean)
    responses <- aggregated$x
    stimuli <- aggregated$stimuli
  }

  # Retrieve parameters of actually presented stimuli (this will work with
  # non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]

  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }

  # Check whether number of parameters are 4096 (this was the case in older
  # versions of rcicr) and should be truncated to 4092 to work well in this new
  # version
  if (!is.vector(params)) {
    if (ncol(params) == 4096) {
      params <- params[, 1:4092]
    }
  } else {
    # In case we only have a single trial as input
    if (length(params) == 4092) {
      params <- params[1:4092]
    }
  }

  # Generate CI(s) ----------------------------------------------------------

  # Invert parameters if antiCI is to be generated
  if (antiCI==TRUE) {
    params = -params
  }

  # If "participants" argument is not given, compute one CI based on all data
  if (all(is.na(participants))) {
    ci <- generateCINoise(params, responses, p)
  # If it is given, create a CI for each particpant and a group CI by averaging
  # across participants
  } else {

    # First generate a CI for each participant, then average across participants
    pids <- as.numeric(factor(participants))
    npids <- length(unique(pids))

    # Initialize progress bar
    pb <- txtProgressBar(min = 1, max = npids, style = 3)

    # Create cluster for parallel processing
    cl <- parallel::makeCluster(n_cores, outfile = '')
    doParallel::registerDoParallel(cl)

    # For each weighted stimulus, construct the noise pattern
    pid.cis <- foreach::foreach(obs = 1:npids,
                                .combine = 'c',
                                .packages = 'rcicr') %dopar% {

      # Update progress bar
      setTxtProgressBar(pb, obs)

      # Select only the observations of the current participant
      pid.rows <- pids == obs

      # Construct the noise pattern
      ci <- generateCINoise(params[pid.rows,], responses[pid.rows], p)

      # Check if individual CIs should be saved. If so, generate and save them
      if (save_individual_cis) {
        if (!is.na(mask)) {
          individual_ci <- applyMask(ci, mask)
        } else {
          individual_ci <- ci
        }
        scaled <- applyScaling(base, individual_ci, individual_scaling,
                               individual_scaling_constant)
        combined <- combine(scaled, base)
        saveToImage(baseimage, combined, paste0(targetpath, '/individual_cis'),
                    unique(participants)[obs], antiCI)
      }

      # Return the CI
      return(ci)
    }
    parallel::stopCluster(cl)
    dim(pid.cis) <- c(img_size, img_size, npids)

    # Average across participants for final CI and return to original variance
    ci <- apply(pid.cis, c(1, 2), mean) #* sqrt(npids)
  }

  # Check if a mask has been set. If so, apply it to the CI
  if (!is.na(mask)) {
    ci <- applyMask(ci, mask)
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
  if(zmapbool) {

    if(zmapmethod == 'quick') {
      # Blur CI
      zmap <- as.matrix(blur(as.im(ci), sigma = sigma))

      # Create z-map
      zmap <- matrix(scale(as.vector(zmap)), img_size, img_size)

      # Apply threshold
      zmap[zmap > -threshold & zmap < threshold] <- NA
    }

    if(zmapmethod == 't.test') {

      # Compute one CI in one single step based on all data
      if (all(is.na(participants))) {
        # Weigh the stimulus parameters of each trial using the given responses
        weightedparameters <- params * responses

        # Get number of observations
        n_observations <- length(responses)

        # Initialize progress bar
        pb <- txtProgressBar(min = 1, max = n_observations, style = 3)

        # Create cluster for parallel processing
        cl <- parallel::makeCluster(n_cores, outfile = '')
        doParallel::registerDoParallel(cl)

        # For each weighted stimulus, construct the complementary noise pattern
        noiseimages <- foreach::foreach(obs = 1:n_observations, .combine = 'c',
                                        .packages = 'rcicr') %dopar% {
          generateNoiseImage(weightedparameters[obs, ], p)
          # Update progress bar
          setTxtProgressBar(pb, obs)
        }
        parallel::stopCluster(cl)
        dim(noiseimages) <- c(img_size, img_size, n_observations)

      } else {
        noiseimages <- pid.cis
      }

      # Get p value for each pixel
      pmap <- apply(noiseimages, 1:2, function(x) unlist(t.test(x)['p.value']))

      # Create Z-map
      zmap <- sign(ci) * abs(qnorm(pmap/2))
    }
    # Pass zmap object to plotZmap for plotting
    plotZmap(zmap = zmap, bgimage = combined, filename = baseimage,
             sigma = sigma, threshold = threshold, size = img_size,
             decoration = zmapdecoration)
  }

  # Return data
  if (zmapbool == T) {
    return(list(ci=ci, scaled=scaled, base=base, combined=combined, zmap=zmap))
  } else {
    return(list(ci=ci, scaled=scaled, base=base, combined=combined))
  }
}

# Functions ---------------------------------------------------------------

# Apply masking to a CI
# Input: CI, mask (either a string or a matrix)
# Output: masked CI (input CI, but masked pixels are NA)
applyMask <- function(ci, mask) {
  # If mask argument is a string, treat it as a path to a bitmap and try to read
  # it into a matrix. If it is a matrix, us that. Else, throw an error.
  if (typeof(mask) == 'character') {
    mask_matrix <- png::readPNG(mask)

    # Check if the PNG uses a greyscale color palette
    if (length(dim(mask_matrix)) != 2) {
      # If the PNG uses the RGB color palette but the image itself is totally
      # greyscale (i.e. the red, green and blue color channels are identical),
      # read it in anyway
      # Thanks https://stackoverflow.com/a/30850654
      rgb_channels <- list(mask_matrix[,,1], mask_matrix[,,2], mask_matrix[,,3])
      if (all(sapply(rgb_channels, FUN = identical, mask_matrix[,,1]))) {
        mask_matrix <- mask_matrix[,,1]
      }
      # Else, throw error
      stop(paste0('This PNG is not encoded with a greyscale color palette and ',
                  'could not be converted to this encoding either. In other ',
                  'words, this is not a greyscale image.'))
    }
  } else if (typeof(mask) == 'double' && length(dim(mask)) == 2) {
    mask_matrix <- mask
  } else {
    stop('The mask argument is neither a string nor a matrix!')
  }

  # Check if mask is of the same size as the stimuli (i.e. img_size)
  if (!all(dim(mask_matrix) == 512)) {
    stop(paste0('Mask is not of the same dimensions as the stimuli! ',
                '(stimulus dimensions: ', img_size, ' x ', img_size,
                '; mask dimensions: ', dim(mask_matrix)[2],
                ' by ', dim(mask_matrix)[1], ').'))
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
    scaled <- (ci + constant) / (2*constant)
    if (max(scaled[!is.na(scaled)]) > 1.0 | min(scaled[!is.na(scaled)]) < 0) {
      warning(paste0('Chosen constant value for constant scaling made noise ',
                     'of classification image exceed possible intensity range ',
                     'of pixels (<0 or >1). Choose a lower value, or clipping ',
                     'will occur.'))
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
    }  else {
      constant <- abs(range(ci[!is.na(ci)])[2])
    }

    scaled <- (ci + constant) / (2*constant)
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
  dir.create(targetpath, recursive = T, showWarnings = F)

  # Write CI to image file
  png::writePNG(combined, paste0(targetpath, '/', filename))
}