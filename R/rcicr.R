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
#' @import png
#' @param cis List of cis, each of which are a list containing the pixel matrices of at least the noise pattern (\code{$ci}) and if the noise patterns need to be written to PNGs, also the base image (\code{$base}).
#' @param saveaspngs Boolean, when set to true, the autoscaled noise patterns will be combined with their respective base images and saved as PNGs (using the key of the list as name).
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @return List of scaled noise patterns and determined scaling factor.
autoscale <- function(cis, saveaspngs=TRUE, targetpath='./cis') {

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

    # Combine and save to PNG if necessary
    if (saveaspngs) {
      ci <- (cis[[ciname]]$scaled + cis[[ciname]]$base) / 2

      dir.create(targetpath, recursive=T, showWarnings = F)

      png::writePNG(ci, paste0(targetpath, '/', ciname, '_autoscaled.png'))
    }

  }

  return(cis)
}

#' Generates classification image
#'
#' Generate classification image for for any reverse correlation task.
#'
#' This function saves the classification image as PNG to a folder and returns the CI. Your choice of scaling
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
#' @param zmap Boolean specifying whether a z-map should be created (default: TRUE).
#' @param zmapmethod String specifying the method to create the z-map. Can be: \code{quick} (default), \code{t.test}.
#' @param zmapdecoration Optional boolean specifying whether the Z-map should be plotted with margins, text (sigma, threshold) and a scale (default: TRUE).
#' @param sigma Integer specifying the amount of smoothing to apply when generating the z-maps (default: 3).
#' @param threshold Integer specifying the threshold z-score (default: 3). Z-scores below the threshold will not be plotted on the z-map.
#' @param zmaptargetpath Optional string specifying path to save z-map PNGs to (default: ./zmaps).
#' @param n_cores Optional integer specifying the number of CPU cores to use to generate the z-map (default: detectCores()).
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined.
generateCI <- function(stimuli, responses, baseimage, rdata, participants=NA, saveaspng=TRUE, filename='', targetpath='./cis', antiCI=FALSE, scaling='independent', constant=0.1, zmap = F, zmapmethod = 'quick', zmapdecoration = T, sigma = 3, threshold = 3, zmaptargetpath = './zmaps', n_cores = detectCores()) {

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
  if (ncol(params) == 4096) {
    params <- params[, 1:4092]
  }

  # Compute classification image #
  if (antiCI) {
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
    cl <- makeCluster(n_cores, outfile = '')
    registerDoParallel(cl)

    # For each weighted stimulus, construct the complementary noise pattern
    pid.cis <- foreach(obs = 1:npids, .combine = 'c', .packages = 'rcicr') %dopar% {
      setTxtProgressBar(pb, obs)
      pid.rows <- pids == obs
      generateCINoise(params[pid.rows,], responses[pid.rows], p)
    }
    stopCluster(cl)
    dim(pid.cis) <- c(img_size, img_size, npids)

    # Average across participants for final CI
    ci <- apply(pid.cis, c(1,2), mean)
  }

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
  if(zmap) {

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
        cl <- makeCluster(n_cores, outfile = '')
        registerDoParallel(cl)

        # For each weighted stimulus, construct the complementary noise pattern
        noiseimages <- foreach(obs = 1:n_observations, .combine = 'c', .packages = 'rcicr') %dopar% {
                                 setTxtProgressBar(pb, obs)
                                 generateNoiseImage(weightedparameters[obs, ], p)
                               }
        stopCluster(cl)
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

#' Plots a Z-map
#'
#' Plots a Z-map given a matrix of z-scores that maps onto a specified base image.
#'
#' This function takes in a matrix of z-scores (as returned by generateCI) and an Rdata file containing a base image. It returns a Z-map image in PNG format.
#' Unlisted additional arguments will be passed to raster::plot. For example, a different color palette can be specified using the \code{col} argument. See raster::plot for details.
#'
#' @export
#' @import dplyr
#' @import viridis
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
plotZmap <- function(zmap, bgimage = '', sigma, threshold = 3, mask, decoration = T, targetpath = 'zmaps', filename = 'zmap', size = 512, ...) {

  # Create target directory
  dir.create(targetpath, recursive = T, showWarnings = F)

  # If a mask is specified, import and check it
  if (exists('mask')) {
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
                 col = viridis(100), ...)
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

#' Generates multiple classification images by participant or condition
#'
#' Generate classification image for any reverse correlation task that displays independently generated alternatives.
#'
#' This function saves the classification images by participant or condition as PNG to a folder and returns the CIs.
#'
#' @export
#' @import dplyr
#' @param data Data frame
#' @param by String specifying column name that specifies the smallest unit (participant, condition) to subset the data on and calculate CIs for.
#' @param stimuli String specifying column name in data frame that contains the stimulus numbers of the presented stimuli.
#' @param responses String specifying column name in data frame that contains the responses coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param saveaspng Boolean stating whether to additionally save the CI as PNG image.
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @param label Optional string to insert in file names of PNGs to make them easier to identify.
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant},  \code{independent} or \code{autoscale} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @return List of classification image data structures (which are themselves lists of pixel matrix of classification noise only, scaled classification noise only, base image only and combined).
batchGenerateCI <- function(data, by, stimuli, responses, baseimage, rdata, saveaspng=TRUE, targetpath='./cis', label='', antiCI=FALSE, scaling='autoscale', constant=0.1) {

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

    # Specify filename for CI PNG
    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[1,by])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[1,by])
    }

    # Compute CI with appropriate settings for this subset (Optimize later so rdata file is loaded only once)
    cis[[filename]] <- generateCI(unitdata[,stimuli], unitdata[,responses], baseimage, rdata, saveaspng, paste0(filename, '.png'), targetpath, antiCI, scaling, constant)
  }

  if (doAutoscale) {
    cis <- autoscale(cis, saveaspngs=saveaspng, targetpath=targetpath)
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

# Suppress checking notes for variables loaded at runtime from .RData files
if(getRversion() >= "2.15.1")  utils::globalVariables(c("p", "s", "base_faces", "stimuli_params"))
