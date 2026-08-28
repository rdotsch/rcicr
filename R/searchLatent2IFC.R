#' Searches latent space across generations of 2IFC trials
#'
#' Run the latent-space 2IFC task repeatedly, moving the point the stimuli vary
#' around towards whatever the participant is choosing, instead of sampling once
#' around a fixed point.
#'
#' \code{\link{generateLatentCI}} estimates the direction a participant prefers.
#' It cannot know how far along it to go, because a single round of trials
#' carries no information about distance: the scaling constant that turns the
#' direction into a face is the researcher's choice. Running the estimate again
#' from the new position answers that. Each generation re-centres on the
#' participant's running optimum, the estimated direction shrinks as the centre
#' approaches what they are looking for, and the search settles there.
#'
#' The update is the response-weighted mean of a generation's perturbations,
#' rescaled to a step of \code{alpha} standard deviations. That is the
#' evolution-strategy estimate of the gradient of the participant's implicit
#' preference, and it is the same estimator \code{\link{generateLatentCI}}
#' computes, used as a direction of travel rather than as an answer.
#'
#' \subsection{Two ways to run it}{
#'
#' With \code{respond}, the whole search runs in one call: the function is
#' handed each trial's two rendered faces and returns \code{1} or \code{-1}.
#' That is what a live experiment harness or a simulation uses.
#'
#' Without \code{respond}, one generation runs per call. The stimuli and a state
#' file are written and the function returns; you run the trials, then call it
#' again with \code{state} pointing at that file and \code{responses} holding
#' what the participant said. This is the mode that fits a lab, where the trials
#' happen in another program on another day.
#' }
#'
#' \subsection{Reading the diagnostics}{
#'
#' The \code{generations} data frame reports, per generation, the length of the
#' estimated direction and its cosine with the previous generation's. Once the
#' centre reaches what the participant is looking for their responses stop
#' carrying a direction, so successive estimates become unrelated and the cosine
#' falls to around zero. A cosine that is near zero from the first generation
#' means something different: they were never discriminating.
#' }
#'
#' @export
#' @param generator A generator built by \code{\link{latentGeneratorPCA}} or
#'   another of rcicr's latentGenerator functions.
#' @param n_generations Number of generations to run. In resumable mode this is
#'   the number of generations the search will have in total, used to tell you
#'   when it is finished.
#' @param n_trials Number of trials per generation.
#' @param stimulus_path String specifying the directory to save the stimuli and
#'   state files to. Required unless \code{save_as_png} is FALSE and
#'   \code{respond} was given; there is no default path.
#' @param base_latent The point the first generation varies around. Defaults to
#'   NULL, which uses the generator's \code{latent_mean}.
#' @param respond Function called once per trial with a list holding the two
#'   rendered images, the two latents, and the generation and trial numbers,
#'   returning 1 when the original is chosen and -1 when the inverted one is.
#'   Defaults to NULL, which runs the search one generation per call.
#' @param state String pointing at the state file from the previous generation,
#'   in resumable mode.
#' @param responses Vector of responses to the previous generation's trials, in
#'   resumable mode.
#' @param latent_sigma Number scaling the perturbations of the first generation,
#'   in units of the generator's \code{latent_sd}.
#' @param sigma_decay Number multiplying \code{latent_sigma} at each generation.
#'   Defaults to 1, which holds the perturbation size constant. Values below 1
#'   shrink it as the search goes on, which sounds like it should sharpen the
#'   answer and measures worse: shrinking the perturbations shrinks the
#'   difference between the two faces while the participant's own uncertainty
#'   stays where it is, so the responses get noisier just as the search needs
#'   them to get finer. Against simulated observers across a range of internal
#'   noise, no annealing beat every decay tried, and by the largest margin for
#'   the noisiest observers.
#' @param alpha How far the centre moves on the first generation, in standard
#'   deviations. Later generations move \code{alpha / generation^step_decay}.
#' @param step_decay How quickly the step shrinks. The default of 1 gives the
#'   Robbins-Monro schedule \code{alpha / generation}, whose steps sum to
#'   infinity, so the search can reach an answer however far away it starts,
#'   while their squares converge, so the noise in each generation's estimate
#'   averages out instead of accumulating. Values below 1 keep the steps larger
#'   for longer, which reaches a distant answer sooner and settles less
#'   precisely. The schedule is what supplies the distance: see the note in the
#'   source of \code{searchStepSize} for why the responses themselves cannot.
#' @param seed Integer seeding the random number generator (for
#'   reproducibility).
#' @param batch_size Number of stimuli rendered per call into the generator.
#' @param save_as_png Boolean specifying whether to write each generation's
#'   stimuli as images.
#' @return In callback mode, a list with the final latent, the trajectory of
#'   centres, the per-generation diagnostics, and the rendered first and last
#'   faces. In resumable mode, a list with the state file to pass back, the
#'   generation just written, and how many remain.
#' @examples
#' # This function is part of the experimental latent-space module.
#' options(rcicr.experimental = TRUE)
#'
#' faces <- replicate(8, tempfile(fileext = ".png"))
#' for (i in seq_along(faces)) {
#'   png::writePNG(matrix(runif(16 * 16), 16, 16), faces[i])
#' }
#' generator <- latentGeneratorPCA(faces, n_components = 4, img_size = 16)
#'
#' # A simulated observer that prefers brighter faces.
#' brighter <- function(trial) if (mean(trial$original) > mean(trial$inverted)) 1 else -1
#'
#' search <- searchLatent2IFC(
#'   generator, n_generations = 3, n_trials = 10,
#'   respond = brighter, save_as_png = FALSE
#' )
#' search$generations
searchLatent2IFC <- function(generator, n_generations = 10, n_trials = 30, stimulus_path, base_latent = NULL, respond = NULL, state = NULL, responses = NULL, latent_sigma = 1, sigma_decay = 1, alpha = 2, step_decay = 1, seed = 1, batch_size = 32, save_as_png = TRUE) {

  requireExperimental('searchLatent2IFC')

  writes_to_disk <- save_as_png || is.null(respond)
  if (writes_to_disk && missing(stimulus_path)) {
    msg <- paste0(
      'No stimulus_path was given. Supply stimulus_path = <a directory> to say ',
      'where the stimuli and the state file should go, or run with a respond ',
      'function and save_as_png = FALSE to search without writing anything.'
    )
    stop(msg, call. = FALSE)
  }

  validateGenerator(generator)

  if (missing(stimulus_path)) {
    stimulus_path <- NA_character_
  } else {
    dir.create(stimulus_path, recursive = TRUE, showWarnings = FALSE)
  }

  if (is.null(respond)) {
    return(searchGenerationStep(generator, n_generations, n_trials, stimulus_path,
                                base_latent, state, responses, latent_sigma,
                                sigma_decay, alpha, step_decay, seed, batch_size,
                                save_as_png))
  }

  if (!is.null(state)) {
    stop('Pass either respond or state, not both.', call. = FALSE)
  }

  set.seed(seed)

  centre <- if (is.null(base_latent)) generator$latent_mean else as.numeric(base_latent)
  checkBaseLatent(centre, generator)

  first_image <- renderLatent(generator, centre, validate = FALSE)[1, , ]
  trajectory <- matrix(centre, nrow = 1)
  diagnostics <- NULL
  previous <- NULL

  for (generation in seq_len(n_generations)) {
    sigma_now <- latent_sigma * sigma_decay^(generation - 1)
    deltas <- drawLatentDeltas(n_trials, generator, sigma_now)

    answers <- collectResponses(generator, centre, deltas, respond, generation,
                                stimulus_path, seed, batch_size, save_as_png)

    direction <- weightedLatentMean(deltas, answers)
    step <- searchStep(direction, generator$latent_sd,
                       searchStepSize(alpha, generation, step_decay))

    diagnostics <- rbind(diagnostics, searchDiagnostics(generation, sigma_now,
                                                        direction, previous,
                                                        generator$latent_sd))
    previous <- direction
    centre <- centre + step
    trajectory <- rbind(trajectory, centre)
  }

  rownames(trajectory) <- NULL

  return(list(
    latent = centre,
    latent_image = renderLatent(generator, centre, validate = FALSE)[1, , ],
    base_latent = trajectory[1, ],
    base_image = first_image,
    trajectory = trajectory,
    generations = diagnostics,
    generator_fingerprint = generator$fingerprint
  ))
}

# One generation, written to disk and handed back for the trials to be run
# elsewhere. The state file carries everything the next call needs, so a search
# survives the weeks between sessions that the .Rdata contract exists for.
searchGenerationStep <- function(generator, n_generations, n_trials, stimulus_path, base_latent, state, responses, latent_sigma, sigma_decay, alpha, step_decay, seed, batch_size, save_as_png) {
  if (is.null(state)) {
    if (!is.null(responses)) {
      stop('responses were given without a state file to apply them to.', call. = FALSE)
    }

    set.seed(seed)
    centre <- if (is.null(base_latent)) generator$latent_mean else as.numeric(base_latent)
    checkBaseLatent(centre, generator)

    previous <- list(generation = 0L, centre = centre, trajectory = matrix(centre, nrow = 1),
                     diagnostics = NULL, direction = NULL)
  } else {
    previous <- advanceSearchState(generator, state, responses, alpha, step_decay)
  }

  generation <- previous$generation + 1L

  if (generation > n_generations) {
    return(list(state = NULL, generation = previous$generation, remaining = 0L,
                latent = previous$centre,
                latent_image = renderLatent(generator, previous$centre, validate = FALSE)[1, , ],
                trajectory = previous$trajectory,
                generations = previous$diagnostics,
                generator_fingerprint = generator$fingerprint))
  }

  sigma_now <- latent_sigma * sigma_decay^(generation - 1)
  deltas <- drawLatentDeltas(n_trials, generator, sigma_now)

  if (save_as_png) {
    writeLatentStimuli(generator, previous$centre, deltas, stimulus_path,
                       searchLabel(generation), seed, batch_size)
  }

  search_state <- list(
    generation = generation,
    centre = previous$centre,
    latent_params = deltas,
    latent_sigma = sigma_now,
    trajectory = previous$trajectory,
    diagnostics = previous$diagnostics,
    direction = previous$direction,
    generator_spec = generatorSpec(generator),
    n_generations = n_generations,
    alpha = alpha,
    seed = seed
  )

  file <- file.path(stimulus_path,
                    sprintf('search_gen%03d_seed_%s.Rdata', generation, seed))
  save(search_state, file = file)

  return(list(
    state = file,
    generation = generation,
    remaining = n_generations - generation,
    latent_params = deltas,
    centre = previous$centre
  ))
}

# The previous generation's state plus its responses, turned into the centre the
# next generation varies around.
advanceSearchState <- function(generator, state, responses, alpha, step_decay) {
  env <- new.env(parent = emptyenv())
  load(state, envir = env)

  if (!exists('search_state', envir = env, inherits = FALSE)) {
    msg <- paste0(
      'The file given as state does not hold a search state. Pass the file ',
      'searchLatent2IFC() returned as $state.'
    )
    stop(msg, call. = FALSE)
  }

  stored <- get('search_state', envir = env, inherits = FALSE)
  matchGenerator(generator, stored$generator_spec, 'searchLatent2IFC')

  if (is.null(responses)) {
    msg <- paste0(
      'Generation ', stored$generation, ' has stimuli waiting but no responses. ',
      'Pass responses = <the ', nrow(stored$latent_params), ' answers>, coded 1 ',
      'for the original image and -1 for the inverted one.'
    )
    stop(msg, call. = FALSE)
  }

  responses <- unlist(responses, use.names = FALSE)
  if (length(responses) != nrow(stored$latent_params)) {
    msg <- paste0(
      'Generation ', stored$generation, ' had ', nrow(stored$latent_params),
      ' trials but ', length(responses), ' responses were given.'
    )
    stop(msg, call. = FALSE)
  }

  direction <- weightedLatentMean(stored$latent_params, responses)
  step <- searchStep(direction, generator$latent_sd,
                     searchStepSize(alpha, stored$generation, step_decay))
  centre <- stored$centre + step

  latest <- searchDiagnostics(stored$generation, stored$latent_sigma, direction,
                              stored$direction, generator$latent_sd)
  diagnostics <- rbind(stored$diagnostics, latest)

  return(list(
    generation = stored$generation,
    centre = centre,
    trajectory = rbind(stored$trajectory, centre),
    diagnostics = diagnostics,
    direction = direction
  ))
}

# How far the centre moves on generation g.
#
# A 2IFC response is a sign, and a sign carries no magnitude: an observer who is
# barely sure and one who is certain answer identically. Working out what the
# response-weighted mean converges to shows the consequence exactly. For
# perturbations of covariance S and a preference direction u, it converges to
# sqrt(2/pi) * S %*% u / sqrt(t(u) %*% S %*% u), whose length in standard
# deviations is the same whatever the length of u. So the estimate says which
# way to go and nothing at all about how far, and a step proportional to it
# moves the same distance however close the centre already is. Measured, that
# leaves the search circling the answer at a fixed radius: steps of 0.93, 1.05,
# 1.15, 1.05, 0.85, 1.06, 1.05, 1.28 standard deviations across eight
# generations, with the distance to the target bouncing between 0.3 and 1.0
# instead of settling.
#
# The distance has to come from the schedule instead. alpha / g^step_decay is
# the Robbins-Monro form: the steps still sum to infinity, so the search can
# reach an answer however far away it starts, while their squares converge, so
# the noise in each estimate averages out rather than accumulating.
searchStepSize <- function(alpha, generation, step_decay) {
  return(alpha / generation^step_decay)
}

searchStep <- function(direction, latent_sd, alpha) {
  return(applyLatentScaling(direction, 'sd', alpha, latent_sd))
}

drawLatentDeltas <- function(n_trials, generator, latent_sigma) {
  deltas <- matrix(stats::rnorm(n_trials * generator$latent_dim), nrow = n_trials)

  return(sweep(deltas, 2, latent_sigma * generator$latent_sd, '*'))
}

checkBaseLatent <- function(base_latent, generator) {
  if (length(base_latent) != generator$latent_dim) {
    msg <- paste0(
      'base_latent has ', length(base_latent), ' dimension(s), but this ',
      'generator\'s latent_dim is ', generator$latent_dim, '.'
    )
    stop(msg, call. = FALSE)
  }

  return(invisible(TRUE))
}

# The callback sees the two faces rather than the two latents, because that is
# what a participant sees and what a real harness has to work from.
collectResponses <- function(generator, centre, deltas, respond, generation, stimulus_path, seed, batch_size, save_as_png) {
  if (save_as_png) {
    writeLatentStimuli(generator, centre, deltas, stimulus_path,
                       searchLabel(generation), seed, batch_size)
  }

  n_trials <- nrow(deltas)
  base <- matrix(centre, nrow = n_trials, ncol = length(centre), byrow = TRUE)
  originals <- renderLatent(generator, base + deltas, validate = FALSE)
  inverted <- renderLatent(generator, base - deltas, validate = FALSE)

  answers <- numeric(n_trials)
  for (i in seq_len(n_trials)) {
    answers[i] <- respond(list(
      original = originals[i, , ],
      inverted = inverted[i, , ],
      latent_original = centre + deltas[i, ],
      latent_inverted = centre - deltas[i, ],
      generation = generation,
      trial = i
    ))
  }

  if (!all(answers %in% c(-1, 1))) {
    msg <- paste0(
      'respond must return 1 (the original was chosen) or -1 (the inverted ',
      'image was) for every trial. Generation ', generation, ' returned ',
      paste(utils::head(unique(setdiff(answers, c(-1, 1))), 3), collapse = ', '), '.'
    )
    stop(msg, call. = FALSE)
  }

  return(answers)
}

searchLabel <- function(generation) {
  return(sprintf('rcic_search_gen%03d', generation))
}

# A cosine near zero after several generations says the centre has arrived: the
# participant's responses no longer carry a direction, so successive estimates
# are unrelated to each other.
searchDiagnostics <- function(generation, latent_sigma, direction, previous, latent_sd) {
  cosine <- NA_real_
  if (!is.null(previous)) {
    denominator <- sqrt(sum(previous^2)) * sqrt(sum(direction^2))
    if (denominator > 0) {
      cosine <- sum(previous * direction) / denominator
    }
  }

  return(data.frame(
    generation = generation,
    latent_sigma = latent_sigma,
    direction_norm = latentNorm(direction, latent_sd),
    cosine_with_previous = cosine
  ))
}
