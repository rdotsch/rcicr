# Script by Ron Dotsch, based on Matlab code by Oliver Langner and Python code by Ron Dotsch
# r.dotsch@psych.ru.nl

#' Generate single sinusoid patch
#'
#' @export
#' @import aspace
#' @import matlab
#' @param img_size Integer specifying size of sinusoid patch in number of pixels.
#' @param cycles Integer specifying number of cycles sinusoid should span.
#' @param angle Value specifying the angle (rotation) of the sinusoid.
#' @param phase Value specifying phase of sinusoid.
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid.
#' @return The sinusoid image with size \code{img_size}.
#' @examples
#' generateSinusoid(512, 2, 90, pi/2, 1.0)
generateSinusoid <- function(img_size, cycles, angle, phase, contrast) {

  # Generates an image matrix containing a sinusoid, angle (in degrees) of 0 will give vertical, 90 horizontally oriented sinusoid
  angle <- aspace::as_radians(angle)
  sinepatch = matlab::repmat(matlab::linspace(0, cycles, img_size), img_size, 1)
  sinusoid <- (sinepatch * cos(angle) + t(sinepatch) * sin(angle)) * 2 * pi
  sinusoid <- contrast * sin(sinusoid + phase)
  return(sinusoid)
}



#' Generate single gabor patch
#'
#' @export
#' @import matlab
#' @import scales
#' @param img_size Integer specifying size of gabor patch in number of pixels.
#' @param cycles Integer specifying number of cycles the sinusoid should span.
#' @param angle Value specifying the angle (rotation) of the sinusoid.
#' @param phase Value specifying phase of sinusoid.
#' @param sigma of guassian mask on top of sinusoid.
#' @param contrast Value between -1.0 and 1.0 specifying contrast of sinusoid.
#' @return The sinusoid image with size \code{img_size}.
#' @examples
#' generateSinusoid(512, 2, 90, pi/2, 1.0)
generateGabor <- function(img_size, cycles, angle, phase, sigma, contrast) {

  s <- generateSinusoid(img_size, cycles, angle, phase, contrast)
  x0 <- scales::rescale(1:img_size, to= c(-.5,.5))
  gauss <- matlab::meshgrid(x0, x0)
  gauss_mask = exp( -(((gauss$x^2)+(gauss$y^2)) / (2* (sigma/img_size)^2)) )
  return(gauss_mask * s)

}

#' Generate sinusoid noise pattern
#'
#' @export
#' @import matlab
#' @param img_size Integer specifying size of the noise pattern in number of pixels.
#' @param nscales Integer specifying the number of incremental spatial scales. Defaults to 5. Higher numbers will add higher spatial frequency scales.
#' @param noise_type String specifying noise pattern type (defaults to \code{sinusoid}; other options: \code{gabor}).
#' @param sigma Number specifying the sigma of the Gabor patch if noise_type is set to \code{gabor} (defaults to 25).
#' @param pre_0.3.0 Boolean specifying whether the noise pattern should be created in a way compatible with older versions of rcicr (< 0.3.0). If you are starting a new project, you should keep this at the default setting (FALSE). There is no reason to set this to TRUE, with the sole exception to recreate behavior of rcicr prior to version 0.3.0.
#' @return List with two elements: the 3D noise matrix with size \code{img_size}, and an indexing
#' matrix with the same size to easily change contrasts.
#' @examples
#' generateNoisePattern(256)
generateNoisePattern <- function(img_size=512, nscales=5, noise_type='sinusoid', sigma=25, pre_0.3.0=FALSE) {
  # Settings of sinusoids
  orientations <- c(0, 30, 60, 90, 120, 150)
  phases <- c(0, pi/2)
  scales <- 2^(0:(nscales-1))

  # Size of patches per scale
  mg <- matlab::meshgrid(1:img_size, 1:img_size,1:length(scales))
  patchSize = mg$x / mg$y

  # Number of patch layers needed
  nrPatches = length(scales) * length(orientations) * length(phases)

  # Preallocate memory
  patches = matlab::zeros(c(img_size, img_size, nrPatches))
  patchIdx = matlab::zeros(c(img_size, img_size, nrPatches))

  # Counters
  if (pre_0.3.0) {
    co = 0 # patch layer counter
    idx = 0 # contrast index counter
  } else {
    co = 1 # patch layer counter
    idx = 1 # contrast index counter
  }

  for (scale in scales) {
    for (orientation in orientations) {
      for (phase in phases) {
        # Generate single patch
        size <- patchSize[scale, img_size]

        if (noise_type=='gabor') {
          p <- generateGabor(size, 1.5, orientation, phase, sigma, 1)
        } else {
          p <- generateSinusoid(size, 2, orientation, phase, 1)
        }

        # Repeat to fill scale
        patches[,,co] <- matlab::repmat(p, scale)

        # Create index matrix
        for (col in 1:scale) {
          for (row in 1:scale) {

            # Insert absolute index for later contrast weighting
            patchIdx[(size * (row-1) + 1) : (size * row), (size * (col-1) + 1) : (size * col), co] = idx

            # Update contrast counter
            idx = idx + 1

          }
        }

        # Update layer counter
        co = co + 1

      }
    }
  }

  return(list(patches=patches, patchIdx=patchIdx, noise_type=noise_type, generator_version=utils::packageVersion('rcicr')))
}


#' Generate single noise image based on parameter vector
#'
#' @export
#' @param params Vector with each value specifying the contrast of each patch in noise.
#' @param p 3D patch matrix (generated using \code{generateNoisePattern()}).
#' @return The noise pattern as pixel matrix.
#' @examples
#' #params <- rnorm(4092) # generates 4092 normally distributed random values
#' #s <- generateNoisePattern(img_size=256)
#' #noise <- generateNoiseImage(params, p)
generateNoiseImage <- function(params, p) {

  # Abort stimulus generation if number of params doesn't equal number of patches
  if (length(params) != max(p$patchIdx)) {
    stop("Stimulus generation aborted: number of parameters doesn't equal number of patches!")
  }

  if ('sinusoids' %in% names(p)) {
    # Pre 0.3.3 noise pattern, rename for appropriate use
    p <- list(patches=p$sinusoids, patchIdx=p$sinIdx, noise_type='sinusoid')
  }

  noise <- apply(p$patches * array(params[p$patchIdx], dim(p$patches)), 1:2, mean)
  return(noise)

}

#' Generate classification image noise pattern based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid
#'
#' @export
#' @param stimuli Matrix with one row per trial, each row containing the 4092 parameters for the original stimulus.
#' @param responses Vector containing the response to each trial (1 if participant selected original, -1 if participant selected inverted;
#' this can be changed into a scale).
#' @param p 3D patch matrix (generated using \code{generateNoisePattern()}).
#' @return The classification image as pixel matrix.
generateCINoise <- function(stimuli, responses, p) {

  weighted <- stimuli * responses

  # Only aggregate if more than one stimulus/response row
  if(is.null(dim(weighted))) {
    params <- weighted
  } else{
    params <- colMeans(weighted)
  }

  return(generateNoiseImage(params, p))
}

#' Determines optimal scaling constant for a list of CIs
#'
#' @export
#' @import matlab
#' @import jpeg
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to JPEGs, also the base image (\code{$base}).
#' @param saveasjpegs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as JPEGs (using the key of the list as name).
#' @param targetpath Optional string specifying path to save JPEGs to (default: ./cis).
#' @return List of scaled noise patterns and determined scaling factor.
autoscale <- function(cis, saveasjpegs=TRUE, targetpath='./cis') {

  # Get range of each ci
  ranges <- matlab::zeros(length(names(cis)), 2)
  for (ciname in names(cis)) {
    ranges[which(ciname==names(cis)), ] <- range(cis[[ciname]]$ci)
  }

  # Determine the lowest possible scaling factor constant
  if (abs(min(ranges[,1])) > max(ranges[,2])) {
    constant <- abs(min(ranges[,1]))
  }  else {
    constant <- max(ranges[,2])
  }

  write(paste0("Using scaling factor constant:", constant), stdout())

  # Scale all noise patterns
  for (ciname in names(cis)) {
    cis[[ciname]]$scaled <-  (cis[[ciname]]$ci + constant) / (2*constant)

    # Combine and save to JPEG if necessary
    if (saveasjpegs) {
      ci <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2

      dir.create(targetpath, recursive=T, showWarnings = F)

      jpeg::writeJPEG(ci, paste0(targetpath, '/', ciname, '_autoscaled.jpg'), quality=1.0)
    }

  }

  return(cis)
}

#' Generates classification image
#'
#' Generate classification image for for any reverse correlation task.
#'
#' This function saves the classification image as JPEG to a folder and returns the CI. Your choice of scaling
#' matters. The default is \code{'matched'}, and will match the range of the intensity of the pixels to
#' the range of the base image pixels. This scaling is nonlinear and depends on the range of both base image
#' and noise pattern. It is truly suboptimal, because it shifts the 0 point of the noise (that is, pixels that would
#' have not changed base image at all before scaling may change the base image after scaling and vice versa). It is
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
#' @import jpeg
#' @import parallel
#' @import doParallel
#' @import foreach
#' @importFrom stats aggregate t.test
#' @importFrom spatstat blur as.im
#' @importFrom grDevices png dev.off
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli.
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveasjpeg Boolean stating whether to additionally save the CI as JPEG image.
#' @param targetpath Optional string specifying path to save JPEGs to (default: ./cis).
#' @param filename Optional string to specify a file name for the JPEG image.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, or \code{independent} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @param zmap Boolean specifying whether a z-map should be created (default: TRUE).
#' @param zmapmethod String specifying the method to create the z-map. Can be: \code{quick} (default), \code{dotsch}, \code{mangini}, or \code{ethier}.
#' @param sigma Integer specifying the amount of smoothing to apply when generating the z-maps (default: 3).
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param zmaptargetpath Optional string specifying path to save z-map PNGs to (default: ./zmaps).
#' @param n_cores Optional integer specifying the number of CPU cores to use to generate the z-map (default: detectCores()).
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined.
generateCI <- function(stimuli, responses, baseimage, rdata, saveasjpeg=TRUE, filename='', targetpath='./cis', antiCI=FALSE, scaling='independent', constant=0.1, zmap = T, zmapmethod = 'dotsch', sigma = 3, threshold = 3, zmaptargetpath = './zmaps', n_cores = detectCores()) {

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


  # Average responses for each presented stimulus (in case stimuli have been presented multiple times,
  # or group-wise classification images are being calculated, in order to reduce memory and processing
  # load)
  aggregated <- aggregate(responses, by=list(stimuli=stimuli), FUN=mean)
  responses <- aggregated$x
  stimuli <- aggregated$stimuli

  # Retrieve parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]

  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }

  # Check whether number of parameters are 4096 (this was the case in older versions of rcicr)
  # and should be truncated to 4092 to work well in this new version
  if (ncol(params) == 4096) {
    params <- params[, 1:4092]
  }

  # Compute classification image
  if (antiCI) {
    params = -params
  }
  ci <- generateCINoise(params, responses, p)

  # Scale
  if (scaling == 'none') {
    scaled <- ci
  } else if (scaling == 'constant') {
    scaled <- (ci + constant) / (2*constant)
    if (max(scaled) > 1.0 | min(scaled) < 0) {
      warning('Chosen constant value for constant scaling made noise of classification image exceed possible intensity range of pixels (<0 or >1). Choose a lower value, or clipping will occur.')
    }
  } else if (scaling == 'matched') {
    scaled <- min(base) + ((max(base) - min(base)) * (ci - min(ci)) / (max(ci) - min(ci)))

  } else if (scaling == "independent") {

    # Determine the lowest possible scaling factor constant
    if (abs(range(ci)[1]) > abs(range(ci)[2])) {
      constant <- abs(range(ci)[1])
    }  else {
      constant <- abs(range(ci)[2])
    }

    scaled <- (ci + constant) / (2*constant)

  } else {
    warning(paste0('Scaling method \'', scaling, '\' not found. Using none.'))
    scaled <- ci
  }

  # Combine with base image
  combined <- (scaled + base) / 2

  # Save to file
  if (saveasjpeg) {
    if (filename == '') {
      filename <- paste0(baseimage, '.jpg')
    }

    if (antiCI) {
      filename <- paste0('antici_', filename)
    } else {
      filename <- paste0('ci_', filename)
    }

    dir.create(targetpath, recursive=T, showWarnings = F)

    jpeg::writeJPEG(combined, paste0(targetpath, '/', filename))

  }

  # Compute Z-map
  if(zmap) {

    if(zmapmethod == 'simple') {
      # Blur CI
      zmap <- as.matrix(blur(as.im(ci), sigma = sigma))

      # Create z-map
      zmap <- matrix(scale(as.vector(zmap)), img_size, img_size)

      # Apply threshold
      zmap[zmap > -threshold & zmap < threshold] <- 0
    }

    if(zmapmethod == 'dotsch') {
      # Weigh the stimulus parameters of each trial using the given responses
      weightedparameters <- stimuli_params$gender_neutral * responses

      # Get number of observations
      n_observations <- length(responses)

      # Create cluster for parallel processing
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)

      # For each weighted stimulus, construct the resulting noise pattern
      noiseimages <- foreach(obs = 1:n_observations, .combine = 'c',
                             .packages = 'rcicr') %dopar% {
                               generateNoiseImage(weightedparameters[obs, ], p)
                             }
      stopCluster(cl)
      dim(noiseimages) <- c(img_size, img_size, n_observations)

      # Apply a T-test to each pixel and register t and p
      pmap <- apply(noiseimages, c(1,2), function(x) unlist(t.test(x)['p.value']))
      tmap <- apply(noiseimages, c(1,2), function(x) unlist(t.test(x)['statistic']))

      # Plot
      png('pmap.png', width = img_size, height = img_size)
      grid::grid.raster((pmap-min(pmap))/max(pmap-min(pmap)))
      dev.off()
      png('tmap.png', width = img_size, height = img_size)
      grid::grid.raster((tmap-min(tmap))/max(tmap-min(tmap)))
      dev.off()
    }
  }

  # Return list
  return(list(ci=ci, scaled=scaled, base=base, combined=combined, zmap=zmap))
}

#' Generates a Z-map
#'
#' Generates a Z-map given a matrix of z-scores that maps onto a specified base image.
#'
#' This function takes in a matrix of z-scores (as returned by generateCI) and an Rdata file containing a base image. It returns a Z-map image in PNG format.
#'
#' @export
#' @import dplyr
#' @param zmap A matrix containing z-scores that map onto a given base image. zmap and baseimage must have the same dimensions.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param targetpath String specifying path to save the Z-map PNG to.
#' @return Nothing. It writes a Z-map image.
generateZmap <- function(zmap, baseimage, rdata) {
}

#' Generates multiple classification images by participant or condition
#'
#' Generate classification image for any reverse correlation task that displays independently generated alternatives.
#'
#' This function saves the classification images by participant or condition as JPEG to a folder and returns the CIs.
#'
#' @export
#' @import dplyr
#' @param data Data frame
#' @param by String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for.
#' @param stimuli String specifying column name in data frame that contains the stimulus numbers of the presented stimuli.
#' @param responses String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveasjpeg Boolean stating whether to additionally save the CI as JPEG image.
#' @param targetpath Optional string specifying path to save JPEGs to (default: ./cis).
#' @param label Optional string to insert in file names of JPEGs to make them easier to identify.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant},  \code{independent} or \code{autoscale} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @return List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined).
batchGenerateCI <- function(data, by, stimuli, responses, baseimage, rdata, saveasjpeg=TRUE, targetpath='./cis', label='', antiCI=FALSE, scaling='autoscale', constant=0.1) {

  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }

  pb <- dplyr::progress_estimated(length(unique(data[,by])))
  cis <- list()

  for (unit in unique(data[,by])) {

    # Update progress bar
    pb$tick()$print()

    # Get subset of data
    unitdata <- data[data[,by] == unit, ]

    # Specify filename for CI JPEG
    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[1,by])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[1,by])
    }

    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI(unitdata[,stimuli], unitdata[,responses], baseimage, rdata, saveasjpeg, paste0(filename, '.jpg'), targetpath, antiCI, scaling, constant)
  }

  if (doAutoscale) {
    cis <- autoscale(cis, saveasjpegs=saveasjpeg, targetpath=targetpath)
  }

  pb$stop()
  return(cis)

}



#' Computes cumulative trial CIs correlations with final/target CI
#'
#' Computes cumulative trial CIs correlations with final/target CI.
#'
#' Use for instance for plotting curves of trial-final/target CI correlations to estimate how many trials are necessary in your task
#'
#' @export
#' @import dplyr
#' @importFrom stats cor
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli.
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param targetci List Target CI object generated with rcicr functions to correlate cumulative CIs with.
#' @param step Step size in sequence of trials to compute correlations with.
#' @return Vector containing correlation between cumulative CI and final/target CI.
computeCumulativeCICorrelation <- function(stimuli, responses, baseimage, rdata, targetci=list(), step=1) {

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


  # Retrieve parameters of actually presented stimuli (this will work with non-consecutive stims as well)
  params <- stimuli_params[[baseimage]][stimuli,]

  # Check whether parameters were found in this .rdata file
  if (length(params) == 0) {
    stop(paste0('No parameters found for base image: ', base))
  }

  # Compute final classification image if necessary
  if (length(targetci) == 0) {
    finalCI <- generateCINoise(params, responses, p)
  } else {
    finalCI <- targetci$ci
  }

  # Compute correlations with final CI with cumulative CI
  pb <- dplyr::progress_estimated(length(responses))

  correlations <- vector()
  corcounter <- 1
  for (trial in seq(1,length(responses), step)) {
    pb$tick()$print()

    cumCI <- generateCINoise(params[1:trial,], responses[1:trial], p)
    correlations[corcounter] <- cor(as.vector(cumCI), as.vector(finalCI))
    corcounter <- corcounter + 1
  }
  pb$stop()

  # Return correlations
  return(correlations)
}



#' Generates Z-map
#'
#' Generates Z-map for a given classification image.
#'
#' This function saves the Z-map (superimposed on the scaled classification image) to a folder.
#'
#' @export
#' @import parallel
#' @import doParallel
#' @import png
#' @importFrom spatstat blur as.im
#' @param cilist A list object (as returned by generateCI) containing at least the "ci" and "combined" element.
#' @param img_size Integer specifying size of the stimuli used in number of pixels.
#' @param sigma Optional number specifying the level of smoothing applied in calculating the z-map (default: 10).
#' @param threshold Optional number specifying the threshold absolute z-score value used in calculating the z-map (default: 3).
#' @param directional Optional boolean specifying whether different colors should be used for negative and positive z-scores (default: F).
#' @param targetpath Optional string specifying path to save JPEGs to (default: ./zmaps).
#' @param saveasjpeg Optional boolean specifying whether to write the Z-map to a JPEG image (default: TRUE).
#' @return Array representing an RGB image.
generateZmap <- function(cilist, img_size = 512, sigma = 10, threshold = 3, directional = F, saveasjpeg = T, targetpath = './zmaps') {
  # Generate Z-map
  zmap <- as.matrix(blur(as.im(cilist$ci), sigma = sigma))
  zmap <- matrix(scale(as.vector(zmap)), img_size, img_size)
  zmap[zmap > -threshold & zmap < threshold] <- 0

  # Normalize
  zmap_normalized <- zmap/max(abs(zmap))

  # Write as JPEG
  if (saveasjpeg) {
    # Import classification image
    img <- rep(cilist$combined, 3)
    dim(img) <- c(512, 512, 3)
    # Add Z-map
    # Zero out the pixels with z-scores above threshold
    img[,,1][zmap_normalized != 0] <- 0
    img[,,2][zmap_normalized != 0] <- 0
    img[,,3][zmap_normalized != 0] <- 0
    if (directional == T) {
      print('hello')
      # Add absolute values of positive normalized z-scores to red channel and
      # negative normalized z-scores to blue channel
      img[,,1][zmap_normalized > 0] <- zmap_normalized[zmap_normalized > 0]
      img[,,3][zmap_normalized < 0] <- abs(zmap_normalized[zmap_normalized < 0])
    } else {
      print('else')
      zmap_normalized <- abs(zmap_normalized)
      img[,,1][zmap_normalized != 0] <- zmap_normalized[zmap_normalized != 0]
    }

    dir.create(targetpath, recursive=T, showWarnings = F)

    # Write
    writePNG(img, paste0(targetpath, '/', 'zmap.png'))
  }

  # Return Z-map
  return(zmap)
}


# Suppress checking notes for variables loaded at runtime from .RData files
if(getRversion() >= "2.15.1")  utils::globalVariables(c("p", "s", "base_faces", "stimuli_params"))
