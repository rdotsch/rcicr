#' Generates 2IFC stimuli
#'
#' Generate stimuli for 2 images forced choice reverse correlation task.
#'
#' Will save the stimuli as PNGs to a folder, including .Rdata file needed for analysis of data
#' after data collection. This .Rdata file contains the parameters that were used to generate each stimulus.
#'
#' @export
#' @import jpeg
#' @import png
#' @import foreach
#' @import doSNOW
#' @importFrom stats setNames runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param base_face_files Named list of base face file names used as base images for stimuli, e.g. \code{list(aName = 'baseface.jpg')}. Accepts JPEG and PNG images, recognised by a \code{.png}, \code{.jpg} or \code{.jpeg} extension. Each name labels that base image's stimulus files and indexes the .Rdata file that \code{\link{generateCI}} reads back, so every element must be named, and named uniquely. Each image must be square and exactly \code{img_size} by \code{img_size} pixels: rcicr does not resize base images. All of this is checked before any stimuli are generated, and the message names the offending entry.
#' @param n_trials Number specifying how many trials the task will have (function will generate two images for each trial per base image: original and inverted/negative noise).
#' @param img_size Number specifying the number of pixels that the stimulus image will span horizontally and vertically (will be square, so only one integer needed).
#' @param stimulus_path String specifying the directory to save the stimuli and the .Rdata file to. Required unless both \code{save_as_png} and \code{save_rdata} are FALSE; there is no default path. It is created if it does not exist. Use \code{tempdir()} if you only want to try the function out.
#' @param label Label to prepend to each file for your convenience.
#' @param use_same_parameters Boolean specifying whether for each base image the same set of parameters is used (TRUE) or a unique set is created for each base image (FALSE).
#' @param seed Integer seeding the random number generator (for reproducibility).
#' @param maximize_baseimage_contrast Boolean specifying whether the pixel values of the base image should be rescaled to maximize its contrast. A base image with no contrast at all — every pixel the same value — has nothing to rescale, and is rejected with an error rather than silently turned into an all-\code{NaN} base image. Such an image is still usable with \code{maximize_baseimage_contrast = FALSE}.
#' @param noise_type String specifying noise pattern type (defaults to \code{sinusoid}; other options: \code{gabor}).
#' @param nscales Integer specifying the number of incremental spatial scales. Defaults to 5. Higher numbers will add higher spatial frequency scales.
#' @param sigma Number specifying the sigma of the Gabor patch if noise_type is set to \code{gabor} (defaults to 25).
#' @param ncores Number of CPU cores to use (default: \code{detectCores()-1}; 2 under \code{R CMD check}, per CRAN policy).
#' @param return_as_dataframe Boolean specifying whether to return a data frame with the raw noise of the stimuli that were generated (default: FALSE). Data frame columns represent pixel values, data frame rows represent stimuli. The frame holds one noise image per trial, so it is meaningful only under the default \code{use_same_parameters = TRUE}, where every base image shares a single parameter set and one noise image per trial is all there is. With \code{use_same_parameters = FALSE} and more than one base image, only the first base image's noise is returned; the frame's shape cannot represent trial x base image. Stimuli are still written to disk for every base image either way, and \code{save_rdata = TRUE} records the full parameter set, so nothing is lost from the files themselves.
#' @param save_as_png Boolean specifying whether to write the stimuli as images to disk (default: TRUE).
#' @param save_rdata Boolean specifying whether .RData file with stimulus parameters will be saved (default: TRUE). Note: you always need to save the .RData file so that you can retrieve the stimulus parameters to compute classification images. This function argument exists primarily for internal rcicr use.
#' @return Nothing, everything is saved to files, unless return_as_dataframe is set to TRUE.
#' @examples
#' # a synthetic square grayscale image stands in for a real base face photo
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#'
#' generateStimuli2IFC(
#'   base_face_files = list(face = base_face),
#'   n_trials = 4,
#'   img_size = 32,
#'   stimulus_path = tempdir(),
#'   seed = 1,
#'   ncores = 1,
#'   nscales = 1
#' )
generateStimuli2IFC <- function(base_face_files, n_trials = 770, img_size = 512, stimulus_path, label = 'rcic', use_same_parameters = TRUE, seed = 1, maximize_baseimage_contrast = TRUE, noise_type = 'sinusoid', nscales = 5, sigma = 25, ncores = default_ncores(), return_as_dataframe = FALSE, save_as_png = TRUE, save_rdata = TRUE) {

  # stimulus_path is required, not defaulted: a default path writes to the
  # user's filespace uninvited, which CRAN policy does not allow.
  writes_to_disk <- save_as_png || save_rdata
  if (writes_to_disk && missing(stimulus_path)) {
    stop(paste0('No stimulus_path was given. Supply stimulus_path = <a ',
      'directory> to say where the stimuli and the .Rdata file ',
      'should go. Use tempdir() if you only want to try the ',
      'function out.'
    ))
  }

  validateBaseFaceFiles(base_face_files)

  # Before the noise basis, which is slow at the default 512px, and before the
  # directory is created.
  base_faces <- readBaseImages(
    base_face_files, img_size, maximize_baseimage_contrast, 'generateStimuli2IFC'
  )

  # Initialize #
  p <- generateNoisePattern(img_size, noise_type = noise_type, nscales = nscales, sigma = sigma)

  # Only create the directory when something is written to it.
  # generateReferenceDistribution2IFC() calls this with both save flags FALSE
  # and used to leave a stray ./stimuli directory behind.
  if (writes_to_disk) {
    dir.create(stimulus_path, recursive = TRUE, showWarnings = FALSE)
  } else if (missing(stimulus_path)) {
    # Bind it anyway. stimulus_path appears in the %dopar% body below, and
    # foreach exports every free variable of that body to the workers -- which
    # means get()ting it, and a missing argument aborts there even though the
    # branch that uses it cannot run. Never read.
    stimulus_path <- NA_character_
  }

  # More depends on this call than the stimuli. generateReferenceDistribution2IFC()
  # re-generates stimuli through this function and then draws its simulated
  # responses with runif(), *after* this set.seed() has run - which is what makes
  # every Informational Value reproducible from the stimulus file alone, and is
  # documented as a guarantee in ?generateReferenceDistribution2IFC.
  #
  # So moving or removing this line, or reseeding after it, changes every InfoVal
  # ever computed with this package without touching computeInfoVal2IFC() at all.
  set.seed(seed)

  stimuli_params <- list()

  # Compute number of parameters needed  #
  nparams <- sum(6 * 2 * (2^(0:(nscales - 1)))^2)

  # Generate parameters #
  if (use_same_parameters) {

    # Generate stimuli parameters, one set for all base faces
    params <- matlab::zeros(n_trials, nparams)
    for (trial in 1:n_trials) {
      params[trial, ] <- (runif(nparams) * 2) - 1
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
        stimuli_params[[base_face]][trial, ] <- (runif(nparams) * 2) - 1
      }
    }

  }

  # Generate stimuli
  pb <- txtProgressBar(min = 1, max = n_trials, style = 3)

  # NULL when ncores == 1: the loop below then runs in this process instead of
  # in a one-worker cluster. See startBackend() in parallel.R.
  cl <- startBackend(ncores)
  if (!is.null(cl)) {
    on.exit(stopClusterSafely(cl), add = TRUE)
  }

  stims <- foreach::foreach(
    trial = 1:n_trials, .packages = 'rcicr', .final = function(x) setNames(as.data.frame(x), as.character(1:n_trials)), .combine = 'cbind', .multicombine = TRUE,
    .options.snow = progressOption(pb, cl)
  ) %dopar% {
    # Each iteration only ever needs the noise for its own trial, so this is a
    # plain matrix. It used to write into a preallocated
    # zeros(img_size, img_size, n_trials) array declared before the cluster was
    # created - at the defaults that is a 1.5 GB object (512 x 512 x 770), and
    # because it existed in the parent environment foreach exported a full copy
    # to *every* worker. Each worker then wrote one slice into its own private
    # copy and discarded it, so the memory was pure overhead. See issue #12.
    if (use_same_parameters) {
      # One parameter set is shared by every base face, so any key gives the
      # same values; take the first explicitly rather than relying on `base_face`
      # still holding a value left over from the base-image loop above.
      trial_noise <- generateNoiseImage(stimuli_params[[names(base_faces)[1]]][trial, ], p)
    }

    for (base_face in names(base_faces)) {
      if (!use_same_parameters) {
        # compute noise pattern unique to this base face
        trial_noise <- generateNoiseImage(stimuli_params[[base_face]][trial, ], p)
      }

      # Scale noise (based on simulations, most values fall within this range [-0.3, 0.3], test
      # for yourself with simulateNoiseIntensities())
      stimulus <- ((trial_noise + 0.3) / 0.6)

      # add base face
      combined <- (stimulus + base_faces[[base_face]]) / 2

      # write to file
      if (save_as_png) {
        png::writePNG(combined, paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_ori.png", trial), sep = "_"), sep = '/'))
      }

      # compute inverted stimulus
      stimulus <- ((-trial_noise + 0.3) / 0.6)

      # add base face
      combined <- (stimulus + base_faces[[base_face]]) / 2

      # write to file
      if (save_as_png) {
        png::writePNG(combined, paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_inv.png", trial), sep = "_"), sep = '/'))
      }

      # Return CI
      if (return_as_dataframe) {
        # Advance the bar here too. This return exits the entire foreach body,
        # not just this loop, so the update below it never ran on this path and
        # the bar sat at zero for the whole run (issue #82). It is duplicated
        # rather than hoisted above the loop deliberately: `trial_noise` is
        # reassigned per base face when use_same_parameters is FALSE, so moving
        # this return would change *which* base face's noise is returned.
        if (is.null(cl)) setTxtProgressBar(pb, trial)
        return(as.vector(trial_noise))
      }
    }

    # Serial path only; in parallel the bar is driven from the parent by
    # .options.snow, because this assignment would land in a worker's copy.
    if (is.null(cl)) setTxtProgressBar(pb, trial)

    # The body's value feeds .combine/.final even when it is discarded (it is,
    # unless return_as_dataframe). Return it explicitly: the guard above is
    # NULL-valued when parallel, and cbind()ing NULLs collapses the frame that
    # .final then tries to setNames() to n_trials columns.
    trial
  }
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }
  cl <- NULL

  # Save all to image file (IMPORTANT, this file is necessary to analyze your data later and create classification images)
  #
  # This records which rcicr wrote the file. It was a hardcoded '0.4.0' string
  # from 2016 until 1.2.0, so *every* .Rdata written by 0.4.0 through 1.1.0
  # claims to come from 0.4.0 no matter what actually wrote it. Anything reading
  # this field must therefore treat '0.4.0' as "unknown, somewhere in that
  # range" rather than as a real version, and must accept both a character
  # string (old files) and the package_version object written here (new ones) --
  # note that comparing versions as strings is wrong anyway, since '0.10.0' sorts
  # below '0.4.0'. p$generator_version has always held the real version and is
  # the more trustworthy of the two on any file that has it.
  generator_version <- utils::packageVersion('rcicr')

  if (save_rdata) {
    # nscales and sigma are saved so that anything re-generating this stimulus
    # set later (notably generateReferenceDistribution2IFC(), which builds the
    # infoVal null distribution) reproduces the same noise basis. They were
    # previously omitted, so re-generation silently fell back to the defaults.
    save(base_face_files, base_faces, img_size, label, n_trials, noise_type, nscales, sigma, p, seed, stimuli_params, stimulus_path, use_same_parameters, generator_version, file = paste(stimulus_path, paste(label, "seed", seed, "time", format(Sys.time(), format = "%b_%d_%Y_%H_%M.Rdata"), sep = "_"), sep = '/'), envir = environment())
  }

  # Return CIs
  if (return_as_dataframe) {
    return(stims)
  }
}

# Check base_face_files up front and name the offending entry.
#
# Each of these used to surface far from its cause: a bare stop() carrying an
# empty message, or "attempt to select less than one element in get1index" from
# inside a parallel worker. See issues #124 and #180.
validateBaseFaceFiles <- function(base_face_files) {
  example <- 'e.g. base_face_files = list(aName = "baseface.jpg")'

  if (!is.list(base_face_files)) {
    stop(paste0('base_face_files must be a named list, ', example,
      '. It is of class ',
      paste(class(base_face_files), collapse = '/'), '.'
    ))
  }

  if (length(base_face_files) == 0) {
    stop(paste0('base_face_files is empty. Supply at least one base image, ',
      example, '.'
    ))
  }

  nms <- names(base_face_files)
  if (is.null(nms) || any(is.na(nms) | nms == '')) {
    stop(paste0('Every element of base_face_files must be named. The names ',
      'label the stimulus files and index the .Rdata file that ',
      'generateCI() reads back, so they cannot be left off: ',
      example, '.'
    ))
  }

  if (anyDuplicated(nms)) {
    stop(paste0('base_face_files has duplicate names (',
      paste(unique(nms[duplicated(nms)]), collapse = ', '),
      '). Only the first entry under each name would be used, and ',
      'the rest would be silently skipped, so give every base image ',
      'its own name.'
    ))
  }

  for (base_face in nms) {
    filename <- base_face_files[[base_face]]

    if (!is.character(filename) || length(filename) != 1 || is.na(filename)) {
      stop(paste0('Base image "', base_face, '" must be a single file name, ',
        'but it is of class ', paste(class(filename), collapse = '/'),
        ' and length ', length(filename), '. ', example, '.'
      ))
    }

    if (is.na(baseImageFormat(filename))) {
      stop(paste0('Base image "', base_face, '" (', filename, ') must be a ',
        'PNG or JPEG file, named with a .png, .jpg or .jpeg ',
        'extension.'
      ))
    }

    if (!file.exists(filename)) {
      stop(paste0('Base image "', base_face, '" does not exist: ', filename,
        '. Paths are resolved relative to the working directory, ',
        getwd(), '.'
      ))
    }
  }

  invisible(TRUE)
}
