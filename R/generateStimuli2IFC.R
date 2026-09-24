#' Generates 2IFC stimuli
#'
#' Generate stimuli for a two-image forced-choice reverse correlation task.
#'
#' Saves the stimuli as PNGs, together with an \code{.Rdata} file holding the parameters used to
#' generate each stimulus. Analysing the responses later requires that file.
#'
#' @export
#' @import jpeg
#' @import png
#' @import foreach
#' @import doSNOW
#' @importFrom stats setNames runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param base_face_files Named list of base image files, e.g. \code{list(aName = 'baseface.jpg')}. JPEG and PNG images are accepted, recognised by a \code{.png}, \code{.jpg} or \code{.jpeg} extension. Each name labels that base image's stimulus files and is the key \code{\link{generateCI}} uses to find it in the \code{.Rdata} file, so every element needs a unique name. Each image must be square and exactly \code{img_size} pixels wide: rcicr does not resize base images. All of this is checked before any stimuli are generated, and an error names the offending entry.
#' @param n_trials Number of trials. Each trial gets two images per base image: one with the noise added (original) and one with it subtracted (inverted).
#' @param img_size Width and height of the square stimulus images, in pixels. It must be divisible by \code{2^(nscales - 1)}, because the finest scale tiles the image with that many patches per side.
#' @param stimulus_path Directory to save the stimuli and the \code{.Rdata} file to. Required unless both \code{save_as_png} and \code{save_rdata} are FALSE; there is no default. The directory is created if it does not exist; to just try the function out, use \code{tempdir()}.
#' @param label Label put at the start of each file name.
#' @param use_same_parameters Boolean: all base images share one set of noise parameters (\code{TRUE}) or each gets its own (\code{FALSE}).
#' @param seed Seed for the random number generator, for reproducibility. It is saved in the \code{.Rdata} file, where the default InfoVal reference replays it. With \code{seed = NULL} there is nothing to replay, so InfoVal references for that file need an explicit \code{response_seed}.
#' @param maximize_baseimage_contrast Boolean: rescale the base image's pixel values to maximize its contrast. A base image with no contrast at all, every pixel the same value, cannot be rescaled and is rejected with an error. It can still be used with \code{maximize_baseimage_contrast = FALSE}.
#' @param noise_type Noise pattern type: \code{sinusoid} (default) or \code{gabor}.
#' @param nscales Number of spatial scales (default: 5). Each additional scale adds a higher spatial frequency. \code{img_size} must be divisible by \code{2^(nscales - 1)}.
#' @param sigma Sigma of the Gabor patches when \code{noise_type = 'gabor'} (default: 25).
#' @param ncores Number of CPU cores to use (default: \code{detectCores() - 1}; 2 under \code{R CMD check}, per CRAN policy).
#' @param return_as_dataframe Boolean: return a data frame with the raw noise of the generated stimuli (default: \code{FALSE}), one row per pixel and one column per trial. With the default \code{use_same_parameters = TRUE} every base image shares the same noise, so that is all of it. With \code{use_same_parameters = FALSE} and more than one base image, only the first base image's noise is returned, because one column per trial cannot hold several. The stimuli are still written for every base image, and \code{save_rdata = TRUE} records every parameter set, so nothing is missing from the files.
#' @param save_as_png Boolean: write the stimuli to disk as PNG images (default: \code{TRUE}).
#' @param save_rdata Boolean: save the \code{.Rdata} file with the stimulus parameters (default: \code{TRUE}). Computing classification images needs that file, so keep this \code{TRUE}; the argument exists mainly for internal use. The file is named \code{<label>_seed_<seed>_time_<month>_<day>_<year>_<hour>_<minute>.Rdata}, for the minute the call started. An existing file of that name is never overwritten: the call stops before generating anything. So does a call into the same folder with the same seed, started in the same minute, while another is still running.
#' @return Nothing: everything is saved to files. With \code{return_as_dataframe = TRUE}, the data frame described there.
#' @examples
#' # a synthetic square grayscale image stands in for a real base face photo
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#'
#' generateStimuli2IFC(
#'   base_face_files = list(face = base_face),
#'   n_trials = 4,
#'   img_size = 32,
#'   stimulus_path = tempfile("stimuli"),
#'   seed = 1,
#'   ncores = 1,
#'   nscales = 1
#' )
generateStimuli2IFC <- function(base_face_files, n_trials = 770, img_size = 512, stimulus_path, label = 'rcic', use_same_parameters = TRUE, seed = 1, maximize_baseimage_contrast = TRUE, noise_type = 'sinusoid', nscales = 5, sigma = 25, ncores = default_ncores(), return_as_dataframe = FALSE, save_as_png = TRUE, save_rdata = TRUE) {

  # Read before any setup, which can cross a minute boundary: two calls started
  # in the same minute must reach the same lock.
  started <- stimulusTime()

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
  base_faces <- list()

  for (base_face in names(base_face_files)) {
    # Read base face
    filename <- base_face_files[[base_face]]
    img_format <- baseImageFormat(filename)

    img <- tryCatch(
      if (img_format == 'png') png::readPNG(filename) else jpeg::readJPEG(filename),
      error = function(e) {
        stop(paste0('Base image "', base_face, '" (', filename, ') could not ',
            'be read as ', img_format, ': ', conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    # Check if base face is square. If not, throw an error
    if (dim(img)[1] != dim(img)[2]) {
      stop(paste0('Base image "', base_face, '" (', filename, ') is not ',
        'square! It\'s ', dim(img)[1], ' by ', dim(img)[2],
        ' pixels. Please use a square base face.'
      ))
    }

    # Change base face to greyscale if necessary.
    #
    # Alpha is not a colour, so it is dropped rather than averaged in. Two
    # channels is grey plus alpha and keeps the first; three or four is RGB or
    # RGBA and averages the first three. Averaging every channel turned an
    # opaque black pixel into 0.5 for grey-plus-alpha and 0.25 for RGBA.
    if (length(dim(img)) == 3) {
      channels <- dim(img)[3]
      keep <- if (channels == 2) 1L else seq_len(min(3, channels))
      img <- apply(img[, , keep, drop = FALSE], c(1, 2), mean)
    }

    # Check that the base face matches the requested stimulus size. Automatic
    # resizing used to happen here via biOps, but that dependency was dropped
    # and never replaced, so a mismatch would otherwise surface much later as
    # an opaque "non-conformable arrays" error from inside a parallel worker
    # (when the noise is added to the base image).
    if (nrow(img) != img_size) {
      stop(paste0('Base image "', base_face, '" (', filename, ') is ',
        nrow(img), ' by ', ncol(img), ' pixels, but img_size is ',
        img_size, '. rcicr does not resize base images: please ',
        'either resize the image to ', img_size, ' by ', img_size,
        ' pixels, or call generateStimuli2IFC() with img_size = ',
        nrow(img), '.'
      ))
    }

    # If necessary, rescale to maximize contrast
    if (maximize_baseimage_contrast) {
      # (img - min) / (max - min) is 0/0 on a uniform image, and the all-NaN
      # base face went into the .Rdata unremarked (#176). Only an error here:
      # with the rescale off a flat base image is usable and produces valid
      # stimuli, so rejecting it outright would break a legitimate call.
      if (max(img) == min(img)) {
        stop(paste0('Base image "', base_face, '" (', filename, ') has no ',
          'contrast: every pixel is ', min(img), '. Contrast cannot ',
          'be maximized on a uniform image, and doing so would make ',
          'the base image entirely NaN. Use a base image with some ',
          'variation, or call generateStimuli2IFC() with ',
          'maximize_baseimage_contrast = FALSE.'
        ))
      }

      img <- (img - min(img)) / (max(img) - min(img))
    }

    # Save base image to list
    base_faces[[base_face]] <- img
  }

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

  if (save_rdata) {
    rdata_file <- stimulusRdataPath(stimulus_path, label, seed, started)
    rdata_lock <- acquireStimulusLock(rdata_file, seed, started)
    on.exit(unlink(rdata_lock, recursive = TRUE), add = TRUE)
  }

  # Reference generation replays these parameter draws to preserve the historical
  # response stream. Changing the seeding or draw count here requires revisiting
  # seedResponseStream(); reproducibility also requires the same RNGkind().
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

    # Nothing past the first base face is written when save_as_png is FALSE, and
    # the frame below carries only the first base face's noise either way. This
    # keeps generateReferenceDistribution2IFC()'s no-PNG re-generation at one
    # generateNoiseImage() call per trial however many base faces there are.
    trial_bases <- if (save_as_png) names(base_faces) else names(base_faces)[1]
    returned_noise <- NULL

    for (base_face in trial_bases) {
      if (!use_same_parameters) {
        # compute noise pattern unique to this base face
        trial_noise <- generateNoiseImage(stimuli_params[[base_face]][trial, ], p)
      }

      # The frame holds one noise image per trial, so it can carry only the
      # first base face's; ?generateStimuli2IFC documents that. Captured inside
      # the loop because after it trial_noise holds the *last* base face's noise
      # when use_same_parameters is FALSE.
      if (return_as_dataframe && is.null(returned_noise)) {
        returned_noise <- as.vector(trial_noise)
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
    }

    # Serial path only; in parallel the bar is driven from the parent by
    # .options.snow, because this assignment would land in a worker's copy.
    if (is.null(cl)) setTxtProgressBar(pb, trial)

    # The body's value feeds .combine/.final even when it is discarded (it is,
    # unless return_as_dataframe), and must never be NULL: cbind()ing NULLs
    # collapses the frame that .final then tries to setNames() to n_trials
    # columns.
    if (return_as_dataframe) returned_noise else trial
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
    save(base_face_files, base_faces, img_size, label, n_trials, noise_type, nscales, sigma, p, seed, stimuli_params, stimulus_path, use_same_parameters, generator_version, file = rdata_file, envir = environment())
  }

  # Return CIs
  if (return_as_dataframe) {
    return(stims)
  }
}

# Which reader a base image needs: 'png', 'jpeg', or NA.
#
# Anchored to the extension. The old grepl('png|PNG', filename) matched anywhere
# in the path, so a JPEG under a directory named "png" went to png::readPNG().
baseImageFormat <- function(filename) {
  if (grepl('\\.png$', filename, ignore.case = TRUE)) {
    return('png')
  }
  if (grepl('\\.jpe?g$', filename, ignore.case = TRUE)) {
    return('jpeg')
  }
  NA_character_
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
