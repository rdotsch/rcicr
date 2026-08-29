#' Generates 2IFC stimuli in a generative model's latent space
#'
#' Generate stimuli for a 2 images forced choice reverse correlation task by
#' perturbing a point in the latent space of a generative face model, rather
#' than by adding sinusoid noise to a base image.
#'
#' Each trial draws one random perturbation of the base latent and renders two
#' faces from it: the base latent plus the perturbation, and the base latent
#' minus it. That mirrors the original and inverted stimuli of
#' \code{\link{generateStimuli2IFC}}, and keeps the same response coding, so an
#' existing task script needs no reshaping: \code{1} means the participant chose
#' the original, \code{-1} the inverted one. The pair is also antithetic, which
#' halves the variance of the classification image computed from it.
#'
#' Perturbations are drawn as \code{latent_sigma} times the generator's own
#' \code{latent_sd}, so \code{latent_sigma = 1} moves a face by about as much as
#' real faces differ from each other along each dimension. Smaller values make a
#' harder task with subtler differences between the two images.
#'
#' Like the pixel-noise pipeline, this writes an .Rdata file that later analysis
#' needs, and nothing about the stimuli is recoverable without it. For a
#' generator built by \code{\link{latentGeneratorPCA}} that file carries the
#' generator too, so it is self-contained. For a generator that is too large to
#' store, the file records a fingerprint instead and
#' \code{\link{generateLatentCI}} needs the generator handed back to it.
#'
#' @export
#' @param generator A generator built by \code{\link{latentGeneratorPCA}} or
#'   another of rcicr's latentGenerator functions.
#' @param n_trials Number specifying how many trials the task will have (the
#'   function renders two images for each trial: original and inverted).
#' @param stimulus_path String specifying the directory to save the stimuli and
#'   the .Rdata file to. Required unless both \code{save_as_png} and
#'   \code{save_rdata} are FALSE; there is no default path. It is created if it
#'   does not exist. Use \code{tempdir()} if you only want to try the function
#'   out.
#' @param base_latent The point in latent space the stimuli vary around.
#'   Defaults to NULL, which uses the generator's \code{latent_mean}: the centre
#'   of the space, which for an eigenface generator is the mean face. Supply a
#'   latent of your own to run the task around a particular face.
#' @param label Label to prepend to each file for your convenience.
#' @param latent_sigma Number scaling the perturbations, in units of the
#'   generator's \code{latent_sd}. Named \code{latent_sigma} rather than
#'   \code{sigma} because the pixel-noise pipeline already saves a \code{sigma}
#'   with an unrelated meaning.
#' @param seed Integer seeding the random number generator (for
#'   reproducibility).
#' @param batch_size Number of stimuli rendered per call into the generator.
#'   Larger batches are faster and use more memory.
#' @param save_as_png Boolean specifying whether to write the stimuli as images
#'   to disk (default: TRUE).
#' @param save_rdata Boolean specifying whether the .Rdata file with the
#'   stimulus parameters will be saved (default: TRUE). You always need this
#'   file to compute classification images afterwards.
#' @return Invisibly, a list with the path of the .Rdata file written (or NA),
#'   the matrix of latent perturbations, and the base latent.
#' @examples
#' # This function is part of the experimental latent-space module.
#' options(rcicr.experimental = TRUE)
#'
#' faces <- replicate(6, tempfile(fileext = ".png"))
#' for (i in seq_along(faces)) {
#'   png::writePNG(matrix(runif(16 * 16), 16, 16), faces[i])
#' }
#' generator <- latentGeneratorPCA(faces, n_components = 3, img_size = 16)
#'
#' # Its own directory, not tempdir() itself: the .Rdata file written here
#' # would otherwise sit among any others already there.
#' stimulus_path <- file.path(tempdir(), "rcicr_latent_example")
#'
#' generateStimuliLatent2IFC(
#'   generator = generator,
#'   n_trials = 4,
#'   stimulus_path = stimulus_path,
#'   seed = 1
#' )
generateStimuliLatent2IFC <- function(generator, n_trials = 300, stimulus_path, base_latent = NULL, label = 'rcic_latent', latent_sigma = 1, seed = 1, batch_size = 32, save_as_png = TRUE, save_rdata = TRUE) {

  requireExperimental('generateStimuliLatent2IFC')

  # stimulus_path is required, not defaulted: a default path writes to the
  # user's filespace uninvited, which CRAN policy does not allow.
  writes_to_disk <- save_as_png || save_rdata
  if (writes_to_disk && missing(stimulus_path)) {
    msg <- paste0(
      'No stimulus_path was given. Supply stimulus_path = <a directory> to say ',
      'where the stimuli and the .Rdata file should go. Use tempdir() if you ',
      'only want to try the function out.'
    )
    stop(msg, call. = FALSE)
  }

  validateGenerator(generator)

  if (!is.numeric(n_trials) || length(n_trials) != 1 || n_trials < 1) {
    stop('n_trials must be a single positive number.', call. = FALSE)
  }

  if (is.null(base_latent)) {
    base_latent <- generator$latent_mean
  }
  base_latent <- as.numeric(base_latent)

  if (length(base_latent) != generator$latent_dim) {
    msg <- paste0(
      'base_latent has ', length(base_latent), ' dimension(s), but this ',
      'generator\'s latent_dim is ', generator$latent_dim, '.'
    )
    stop(msg, call. = FALSE)
  }

  if (writes_to_disk) {
    dir.create(stimulus_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Drawn after this line and nowhere else, so the whole stimulus set follows
  # from the seed, the generator and these parameters -- the same guarantee
  # generateStimuli2IFC() makes.
  set.seed(seed)

  n_trials <- as.integer(n_trials)
  latent_params <- matrix(
    stats::rnorm(n_trials * generator$latent_dim), nrow = n_trials
  )
  latent_params <- sweep(latent_params, 2, latent_sigma * generator$latent_sd, '*')

  run_id <- newSearchRunId()

  if (save_as_png) {
    refuseToOverwriteStimuli(stimulus_path, label, seed, 'generateStimuliLatent2IFC')
    writeLatentStimuli(generator, base_latent, latent_params, stimulus_path,
                       label, seed, batch_size)
  }

  rdata_file <- NA_character_
  if (save_rdata) {
    rdata_file <- writeLatentRdata(generator, base_latent, latent_params,
                                   stimulus_path, label, latent_sigma, seed, run_id)
  }

  return(invisible(list(
    rdata = rdata_file,
    latent_params = latent_params,
    base_latent = base_latent
  )))
}

# Rendering runs in this process rather than through the foreach/doSNOW backend
# the pixel pipeline uses. A generator is a single external resource -- a GPU,
# a Python interpreter, a loaded model -- and fanning it out across workers
# either duplicates the model in memory or contends for the device. The speed
# is in render()'s batching instead, which is why the contract makes it take a
# matrix of latents rather than one at a time.
writeLatentStimuli <- function(generator, base_latent, latent_params, stimulus_path, label, seed, batch_size) {
  renderStimuliInBatches(generator, base_latent, latent_params, batch_size,
    function(trials, originals, inverted) {
      writeStimulusBatch(trials, originals, inverted, stimulus_path, label, seed)
    }
  )

  return(invisible(NULL))
}

# One pass over the trials, rendering each batch exactly once and handing it to
# the caller. What a batch is then used for -- written to disk, shown to a
# respond callback, or both -- is the caller's business.
#
# Sharing this rather than rendering separately per destination is not only
# speed. A generator is often an external program, so rendering the same
# stimulus twice doubles the cost of the slowest step; and a generator with any
# nondeterminism would otherwise show a participant pixels that the archived PNG
# does not match, which is unrecoverable after the fact.
renderStimuliInBatches <- function(generator, base_latent, latent_params, batch_size, on_batch, progress = TRUE) {
  n_trials <- nrow(latent_params)
  batch_size <- max(1L, as.integer(batch_size))
  starts <- seq(1, n_trials, by = batch_size)
  collected <- vector('list', length(starts))

  pb <- NULL
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = n_trials, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (b in seq_along(starts)) {
    trials <- seq(starts[b], min(starts[b] + batch_size - 1L, n_trials))
    deltas <- latent_params[trials, , drop = FALSE]

    base <- matrix(base_latent, nrow = length(trials), ncol = length(base_latent), byrow = TRUE)
    originals <- renderUnchecked(generator, base + deltas, validate = FALSE)
    inverted <- renderUnchecked(generator, base - deltas, validate = FALSE)

    collected[[b]] <- on_batch(trials, originals, inverted)

    if (!is.null(pb)) {
      utils::setTxtProgressBar(pb, max(trials))
    }
  }

  return(invisible(collected))
}

writeStimulusBatch <- function(trials, originals, inverted, stimulus_path, label, seed) {
  for (i in seq_along(trials)) {
    png::writePNG(originals[i, , ], latentStimulusFile(stimulus_path, label, seed, trials[i], 'ori'))
    png::writePNG(inverted[i, , ], latentStimulusFile(stimulus_path, label, seed, trials[i], 'inv'))
  }

  return(invisible(NULL))
}

# The %05d trial number is the stimulus index into latent_params, matching the
# pixel pipeline's naming so a task script reads either set the same way.
# Two calls at the same label and seed write the same image names, so the second
# replaces the first participant's stimuli while the first .Rdata survives,
# pointing at perturbations whose faces have been swapped underneath it.
#
# Refused rather than renamed. The trial number in a stimulus filename is part of
# what this module documents -- it is how a task script finds image 42 and how a
# response file is matched back to a perturbation -- so making the names
# unpredictable to keep them unique would break the thing they exist for. An
# error names the collision and the two ways out.
refuseToOverwriteStimuli <- function(stimulus_path, label, seed, caller) {
  existing <- list.files(stimulus_path,
                         pattern = paste0('^', sprintf('%s_%s_', label, seed), '[0-9]{5}_(ori|inv)[.]png$'))

  if (length(existing) > 0) {
    msg <- paste0(
      caller, '() would overwrite ', length(existing), ' stimulus image(s) ',
      'already in ', stimulus_path, ', because they carry the same label ("',
      label, '") and seed (', seed, '). Their .Rdata would survive and point at ',
      'perturbations whose images had been replaced. Give this set its own ',
      'label, or write it to an empty directory.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

latentStimulusFile <- function(stimulus_path, label, seed, trial, suffix) {
  name <- sprintf('%s_%s_%05d_%s.png', label, seed, trial, suffix)

  return(file.path(stimulus_path, name))
}

writeLatentRdata <- function(generator, base_latent, latent_params, stimulus_path, label, latent_sigma, seed, run_id) {
  generator_spec <- generatorSpec(generator)
  img_size <- generator$img_size
  n_trials <- nrow(latent_params)
  generator_version <- utils::packageVersion('rcicr')

  # The pixel pipeline's timestamp resolves to the minute, which is enough there
  # because generating a stimulus set at 512px takes longer than that. Here a set
  # can be made in under a second, so two calls sharing a label and seed within
  # one minute wrote the same file and the first call's returned path then
  # pointed at the second's perturbations. The run id makes the name unique
  # without consuming the random stream that seed reproduces the set from.
  name <- paste(
    label, 'seed', seed, 'time',
    paste0(format(Sys.time(), format = '%b_%d_%Y_%H_%M'), '_', run_id, '.Rdata'),
    sep = '_'
  )
  file <- file.path(stimulus_path, name)

  save(latent_params, base_latent, latent_sigma, generator_spec, n_trials,
       img_size, seed, label, stimulus_path, generator_version,
       file = file, envir = environment())

  return(file)
}
