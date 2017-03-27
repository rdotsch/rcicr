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
#' @param saveaspng Optional Boolean stating whether to additionally save the CI as PNG image.
#' @param participants Optional vector specifying participant IDs. If specified, will compute the requested CIs in two steps: step 1, compute CI for each participant. Step 2, compute final CI by averaging participant CIs. If unspecified, the function defaults to averaging all data in the stimuli and responses vector.
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @param filename Optional string to specify a file name for the PNG image.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, or \code{independent} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @param mask Optional 2D matrix that defines the mask to be applied to the CI (1 = masked, 0 = unmasked). May also be a string specifying the path to a grayscale PNG image (black = masked, white = unmasked). Default: NA.
#' @param zmap Boolean specifying whether a z-map should be created (default: FALSE).
#' @param zmapmethod String specifying the method to create the z-map. Can be: \code{quick} (default), \code{t.test}.
#' @param zmapdecoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param sigma Integer specifying the amount of smoothing to apply when generating the z-maps (default: 3).
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param zmaptargetpath Optional string specifying path to save z-map PNGs to (default: ./zmaps).
#' @param n_cores Optional integer specifying the number of CPU cores to use to generate the z-map (default: detectCores()).
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined.
generateCI <- function(mask=NA, stimuli, responses, baseimage, rdata, participants=NA, saveaspng=TRUE, filename='', targetpath='./cis', antiCI=FALSE, scaling='independent', constant=0.1, zmap = F, zmapmethod = 'quick', zmapdecoration = T, sigma = 3, threshold = 3, zmaptargetpath = './zmaps', n_cores = detectCores()) {

  # Rename zmap to zmapbool so we can use zmap for the actual zmap
  zmapbool <- zmap

  # Load parameter file (created when generating stimuli)
  load(rdata)

  # Check whether critical variables have been loaded
  if (!exists('s', envir=environment(), inherits=FALSE) & !exists('p', envir=environment(), inherits=FALSE) ) {
    stop('File specified in rdata argument did not contain s or p variable.')
  }

  if (!exists('base_faces', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain base_faces variable.')
  }

  if (!exists('stimuli_params', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain stimuli_params variable.')
  }

  if (!exists('img_size', envir=environment(), inherits=FALSE)) {
    stop('File specified in rdata argument did not contain img_size variable.')
  }


  # Convert s to p (if rdata file originates from pre-0.3.3)
  if (exists('s', envir=environment(), inherits=FALSE)) {
    p <- list(patches=s$sinusoids, patchIdx=s$sinIdx, noise_type='sinusoid')
    rm(s)
  }

  # Get base image
  base <- base_faces[[baseimage]]
  if (is.null(base)) {
    stop(paste0('File specified in rdata argument did not contain any reference to base image label: ', baseimage, ' (NOTE: file contains references to the following base image label(s): ', paste(names(base_faces), collapse=', '), ')'))
  }

  if (all(is.na(participants))) {
    # Average responses for each presented stimulus (in case stimuli have been presented multiple times,
    # or group-wise classification images are being calculated, in order to reduce memory and processing
    # load)
    aggregated <- aggregate(responses, by=list(stimuli=stimuli), FUN=mean)
    responses <- aggregated$x
    stimuli <- aggregated$stimuli
  }

  # Retrieve parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]

  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }

  # Check whether number of parameters are 4096 (this was the case in older versions of rcicr)
  # and should be truncated to 4092 to work well in this new version
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

  # Compute classification image #
  if (antiCI==TRUE) {
    params = -params
  }

  if (all(is.na(participants))) {

    # Compute one CI in one single step based on all data
    ci <- generateCINoise(params, responses, p)

  } else {

    # First CI by participant, then average across participants
    pids <- as.numeric(factor(participants))
    npids <- length(unique(pids))

    # Initialize progress bar
    pb <- txtProgressBar(min = 1, max = npids, style = 3)

    # Create cluster for parallel processing
    cl <- parallel::makeCluster(n_cores, outfile = '')
    doParallel::registerDoParallel(cl)

    # For each weighted stimulus, construct the complementary noise pattern
    pid.cis <- foreach::foreach(obs = 1:npids, .combine = 'c', .packages = 'rcicr') %dopar% {
      setTxtProgressBar(pb, obs)
      pid.rows <- pids == obs
      generateCINoise(params[pid.rows,], responses[pid.rows], p)
    }
    parallel::stopCluster(cl)
    dim(pid.cis) <- c(img_size, img_size, npids)

    # Average across participants for final CI and return to original variance
    ci <- apply(pid.cis, c(1,2), mean) #* sqrt(npids)
  }



  # Mask #

  # Check if a mask has been set
  if (!is.na(mask)) {
    # If mask argument is a string, treat it as a path to a bitmap and try to read it into a matrix
    # If mask is a matrix, use it as is
    # Else, throw error
    if (typeof(mask) == 'character') {
      mask_matrix <- png::readPNG(mask)

      # Check if the PNG uses a greyscale color palette
      if (length(dim(mask_matrix)) != 2) {
        # If the PNG uses the full color palette but the image itself is totally greyscale, read it in anyway
        # Thanks https://stackoverflow.com/a/30850654
        if (all(sapply(list(mask_matrix[,,1], mask_matrix[,,2], mask_matrix[,,3]), FUN = identical, mask_matrix[,,1]))) {
          mask_matrix <- mask_matrix[,,1]
        }
        # Else, throw error
        stop('This PNG is not encoded with a greyscale color palette and could not be converted to this encoding either. In other words, this is not a greyscale image.')
      }
    } else if (typeof(mask) == 'double' && length(dim(mask)) == 2) {
      mask_matrix <- mask
    } else {
      stop('The mask argument is neither a string nor a matrix!')
    }

    # Check if mask is of the same size as the stimuli (i.e. img_size)
    if (!all(dim(mask_matrix) == 512)) {
      stop(paste0('Mask is not of the same dimensions as the stimuli! (stimulus dimensions: ', img_size, ' x ', img_size, '; mask dimensions: ', dim(mask_matrix)[2], ' by ', dim(mask_matrix)[1], ').'))
    }

    # Check if the mask is binary
    if (length(mask_matrix) != sum(mask_matrix %in% c(0, 1))) {
      stop('This mask contains values other than 0 or 1!')
    }

    # Convert mask to boolean matrix (black == 0 == masked)
    mask <- mask_matrix == 0

    # Apply the mask to the CI. This replaces all the masked pixels with NA
    ci[mask] <- NA
  }



  # Scale
  if (scaling == 'none') {
    scaled <- ci
  } else if (scaling == 'constant') {
    scaled <- (ci + constant) / (2*constant)
    if (max(scaled[!is.na(scaled)]) > 1.0 | min(scaled[!is.na(scaled)]) < 0) {
      warning('Chosen constant value for constant scaling made noise of classification image exceed possible intensity range of pixels (<0 or >1). Choose a lower value, or clipping will occur.')
    }
  } else if (scaling == 'matched') {
    scaled <- min(base) + ((max(base) - min(base)) * (ci - min(ci[!is.na(ci)])) / (max(ci[!is.na(ci)]) - min(ci[!is.na(ci)])))

  } else if (scaling == "independent") {

    # Determine the lowest possible scaling factor constant
    if (abs(range(ci[!is.na(ci)])[1]) > abs(range(ci[!is.na(ci)])[2])) {
      constant <- abs(range(ci[!is.na(ci)])[1])
    }  else {
      constant <- abs(range(ci[!is.na(ci)])[2])
    }

    scaled <- (ci + constant) / (2*constant)

  } else {
    warning(paste0('Scaling method \'', scaling, '\' not found. Using none.'))
    scaled <- ci
  }

  # Combine with base image
  combined <- (scaled + base) / 2

  # Save to file
  if (saveaspng) {
    if (filename == '') {
      filename <- paste0(baseimage, '.png')
    }

    if (antiCI) {
      filename <- paste0('antici_', filename)
    } else {
      filename <- paste0('ci_', filename)
    }

    dir.create(targetpath, recursive = T, showWarnings = F)

    png::writePNG(combined, paste0(targetpath, '/', filename))

  }

  # Compute Z-map
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
        noiseimages <- foreach::foreach(obs = 1:n_observations, .combine = 'c', .packages = 'rcicr') %dopar% {
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
      pmap <- apply(noiseimages, c(1,2), function(x) unlist(t.test(x)['p.value']))

      # Create Z-map
      zmap <- sign(ci) * abs(qnorm(pmap/2))

    }

    # Pass zmap object to plotZmap for plotting
    plotZmap(zmap = zmap, bgimage = combined, filename = baseimage, sigma = sigma, threshold = threshold, size = img_size, decoration = zmapdecoration)

  }

  # Return list
  if (zmapbool == T) {
    return(list(ci=ci, scaled=scaled, base=base, combined=combined, zmap=zmap))
  } else {
    return(list(ci=ci, scaled=scaled, base=base, combined=combined))
  }
}
