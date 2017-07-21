#' Generates 2IFC stimuli
#'
#' Generate stimuli for 2 images forced choice reverse correlation task.
#'
#' Will save the stimuli as
#' PNGs to a folder, including .Rdata file needed for analysis of data after data collection. This
#' .Rdata file contains the parameters that were used to generate each stimulus.
#'
#' @export
#' @import matlab
#' @import dplyr
#' @import jpeg
#' @import png
#' @import foreach
#' @import doParallel
#' @importFrom stats setNames runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param base_face_files List containing base face file names used as base images for stimuli. Accepts JPEG and PNG images.
#' @param n_trials Number specifying how many trials the task will have (function will generate two images for each trial per base image: original and inverted/negative noise).
#' @param img_size Number specifying the number of pixels that the stimulus image will span horizontally and vertically (will be square, so only one integer needed).
#' @param stimulus_path Path to save stimuli and .Rdata file to.
#' @param label Label to prepend to each file for your convenience.
#' @param use_same_parameters Boolean specifying whether for each base image the same set of parameters is used (TRUE) or a unique set is created for each base image (FALSE).
#' @param seed Integer seeding the random number generator (for reproducibility).
#' @param maximize_baseimage_contrast Boolean specifying wheter the pixel values of the base image should be rescaled to maximize its contrast.
#' @param noise_type String specifying noise pattern type (defaults to \code{sinusoid}; other options: \code{gabor}).
#' @param nscales Integer specifying the number of incremental spatial scales. Defaults to 5. Higher numbers will add higher spatial frequency scales.
#' @param sigma Number specifying the sigma of the Gabor patch if noise_type is set to \code{gabor} (defaults to 25).
#' @param ncores Number of CPU cores to use (default: detectCores()-1).
#' @param returnAsDataframe Boolean specifying whether to return a data frame with the raw noise of the stimuli that were generated (default: FALSE). Data frame columns represent pixel values, data frame rows represent stimuli.
#' @return Nothing, everything is saved to files, unless returnAsDataframe is set to TRUE.
generateStimuli2IFC <- function(base_face_files, n_trials=770, img_size=512, stimulus_path='./stimuli', label='rcic', use_same_parameters=TRUE, seed=1, maximize_baseimage_contrast=TRUE, noise_type='sinusoid', nscales=5, sigma=25, ncores=parallel::detectCores()-1, returnAsDataframe=FALSE) {

  # Initialize #
  p <- generateNoisePattern(img_size, noise_type=noise_type, nscales=nscales, sigma=sigma)
  dir.create(stimulus_path, recursive=T, showWarnings = F)
  set.seed(seed)

  stimuli_params <- list()
  base_faces <- list()

  if (!is.list(base_face_files)) {
    write("Please provide base face file name as named list, e.g. base_face_files=list(aName='baseface.jpg')", stderr())
    stop()
  }

  for (base_face in names(base_face_files)) {
    # Read base face
    filename <- base_face_files[[base_face]]
    if (grepl('png|PNG', filename)) {
      img <- png::readPNG(filename)
    } else if (grepl('jpeg|JPEG|jpg|JPG', filename)) {
      img <- jpeg::readJPEG(filename)
    } else {
      stop(paste0('Error in reading base image file ',
                  filename, ': must be a PNG or JPEG file.'))
    }

    # Check if base face is square. If not, throw an error
    if (dim(img)[1] != dim(img)[2]) {
      stop(paste0('Base face is not square! It\'s ', dim(img)[1], ' by ',
                  dim(img)[2], ' pixels. Please use a square base face.'))
    }

    # Change base face to greyscale if necessary
    if (length(dim(img)) == 3) {
      img <- apply(img, c(1, 2), mean)
    }

    # Adjust size of base face
    #base_faces[[base_face]] <- biOps::imgMedianShrink(img, x_scale=img_size/ncol(img), y_scale=img_size/nrow(img))

    # If necessary, rescale to maximize contrast
    if (maximize_baseimage_contrast) {
      img <- (img - min(img)) / (max(img) - min(img))
    }

    # Save base image to list
    base_faces[[base_face]] <- img
  }

  # Compute number of parameters needed  #
  nparams <- sum(6*2*(2^(0:(nscales-1)))^2)

  # Generate parameters #
  if (use_same_parameters) {

    # Generate stimuli parameters, one set for all base faces
    params <- matlab::zeros(n_trials, nparams)
    for (trial in 1:n_trials) {
      params[trial,] <- (runif(nparams) * 2) - 1
    }

    # Assign to each base face the same set
    for (base_face in names(base_faces)) {
      stimuli_params[[base_face]] <- params
    }

    rm(params)
  } else {
    for (base_face in names(base_faces)) {
      # Generate stimuli parameters, unique to each base face
      stimuli_params[[base_face]] <- matlab::zeros(n_trials, nparams)
      for (trial in 1:n_trials) {
        stimuli_params[[base_face]][trial,] <- (runif(nparams) * 2) - 1
      }
    }

  }

  # Generate stimuli
  pb <- txtProgressBar(min = 1, max = n_trials, style = 3)

  stimuli <- matlab::zeros(img_size, img_size, n_trials)

  cl <- parallel::makeCluster(ncores, outfile = "")
  doParallel::registerDoParallel(cl)

  stims <- foreach::foreach(
    trial = 1:n_trials, .packages = 'rcicr', .final = function(x) setNames(as.data.frame(x), as.character(1:n_trials)), .combine = 'cbind', .multicombine = TRUE) %dopar% {
    if (use_same_parameters) {
      # compute noise pattern, can be used for all base faces
      stimuli[,,trial] <- generateNoiseImage(stimuli_params[[base_face]][trial,], p)
    }

    for (base_face in names(base_faces)) {
      if (!use_same_parameters) {
        # compute noise pattern unique to this base face
        stimuli[,,trial] <- generateNoiseImage(stimuli_params[[base_face]][trial,], p)
      }

      # Scale noise (based on simulations, most values fall within this range [-0.3, 0.3], test
      # for yourself with simulateNoiseIntensities())
      stimulus <- ((stimuli[,,trial] + 0.3) / 0.6)

      # add base face
      combined <- (stimulus + base_faces[[base_face]]) / 2

      # write to file
      png::writePNG(combined, paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_ori.png", trial), sep="_"), sep='/'))

      # compute inverted stimulus
      stimulus <- ((-stimuli[,,trial] + 0.3) / 0.6)

      # add base face
      combined <- (stimulus + base_faces[[base_face]]) / 2

      # write to file
      png::writePNG(combined, paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_inv.png", trial), sep="_"), sep='/'))

      # Return CI
      if (returnAsDataframe) {
        return(as.vector(stimuli[,,trial]))
      }
    }

    # Update progress bar
    setTxtProgressBar(pb, trial)
  }
  parallel::stopCluster(cl)

  # Save all to image file (IMPORTANT, this file is necessary to analyze your data later and create classification images)
  generator_version <- '0.3.3'

  save(base_face_files, base_faces, img_size, label, n_trials, noise_type, p, seed, stimuli_params, stimulus_path, trial, use_same_parameters, generator_version, file=paste(stimulus_path, paste(label, "seed", seed, "time", format(Sys.time(), format="%b_%d_%Y_%H_%M.Rdata"), sep="_"), sep='/'), envir=environment())

  # Return CIs
  if (returnAsDataframe) {
    return(stims)
  }
}
