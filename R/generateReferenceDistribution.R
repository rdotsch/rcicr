#' Generates reference distribution
#'
#' Generates reference distribution of norms for a particular set of task parameters.
#'
#' In order to compute the Informational Value metric. Saves its results in the supplied rdata file for later reuse.
#'
#' @section Reproducibility:
#' With the default \code{response_seed = NULL}, the reference distribution is determined by
#' the stimulus \code{.Rdata} file and the session's \code{\link{RNGkind}}. It does not
#' depend on the ambient random number state, and it does not depend on \code{ncores}. Two
#' researchers who compute InfoVal from the same stimulus file therefore get the same number,
#' and the same reference distribution, on different machines and in different sessions --
#' provided both sessions use the same RNG kind.
#'
#' \code{set.seed()} keeps whatever kind the session already has rather than restoring one,
#' and no stimulus file records which was in force, so a session that has changed
#' \code{RNGkind()} changes the simulated response draws and therefore the null. The saved
#' noise basis and stimulus parameters remain unchanged. To reproduce an earlier reference,
#' use the RNG kind that built it; see \url{https://github.com/rdotsch/rcicr/issues/315}.
#'
#' The noise is reconstructed from the basis and parameters the stimulus file saved, so the
#' base images themselves are never reopened and an archived or moved experiment can still be
#' scored. Responses are seeded from the state following one parameter matrix's draws at the
#' saved stimulus seed, preserving the historical default response stream.
#'
#' Pass an explicit \code{response_seed} to draw a *different* null from the same stimuli --
#' for instance to check how much Monte Carlo error a given \code{iter} leaves in your
#' InfoVal. This changes only the simulated responses; the stimuli themselves, and so the
#' noise basis the null is built on, are unaffected.
#'
#' @export
#' @importFrom stats runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param iter Number of iterations for the simulation (i.e., the number of norms generated with classification images based on random responding).
#' @param ncores Number of CPU cores to use when rebuilding the saved noise (default: \code{detectCores()-1}; 2 under \code{R CMD check}, per CRAN policy).
#' @param response_seed Optional seed for the simulated random responses. The default
#' (\code{NULL}) draws them from the state the stimulus generator left behind, which is the
#' reproducible behaviour described under Reproducibility. Supply a number to obtain an
#' independent draw of the null from the same stimuli.
#' @param save_rdata Boolean specifying whether the reference distribution should be written
#' back into the \code{rdata} file (default \code{TRUE}). Set to \code{FALSE} to compute a
#' distribution without changing what later calls to \code{\link{computeInfoVal2IFC}} will
#' use -- worth doing whenever \code{response_seed} is set, so a one-off null does not become
#' the file's permanent reference.
#' @param baseimage Saved base-image label, using the same key as \code{generateCI()}.
#' Required when the saved base images have different noise parameters. With a single base
#' or identical parameter matrices, \code{NULL} retains the shared reference behavior.
#' @section Independent base images:
#' When saved parameter matrices differ, supply \code{baseimage} explicitly to say which
#' base's noise to use. Cached distributions are stored in \code{reference_norms_by_base},
#' keyed by base label, with \code{norms} and \code{response_seed} in each entry. Old unscoped
#' \code{reference_norms} are neither reused nor overwritten for independent bases.
#' Existing shared-parameter files continue using their unscoped cache and reconstruction.
#' @return The reference distribution, invisibly, as a numeric vector of \code{iter} norms.
#' Unless \code{save_rdata = FALSE}, it is also added to the supplied \code{rdata} file as
#' \code{reference_norms} (alongside \code{reference_norms_seed}, recording the
#' \code{response_seed} it was generated with), so a later call to
#' \code{\link{computeInfoVal2IFC}} using the same file can reuse it instead of re-simulating.
#' Independent-base references instead use \code{reference_norms_by_base}, as described above.
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
#' # iter is kept tiny here for a fast example; in practice use iter >= 10000.
#' suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3, ncores = 1))
generateReferenceDistribution2IFC <- function(rdata, iter = 10000, ncores = default_ncores(), response_seed = NULL, save_rdata = TRUE, baseimage = NULL) { # nolint: object_length_linter.

  reference_selection <- selectReferenceBase(rdata, baseimage)
  if (reference_selection$independent) {
    return(invisible(generateBaseReference(reference_selection, rdata, iter,
                                           ncores, response_seed, save_rdata)))
  }

  # load() assigns straight into this function's frame, so any object stored in
  # the .Rdata file silently overwrites an argument of the same name. This
  # function re-saves its frame at the end, so files it has already written
  # contain `rdata` (and, since ncores was added, `ncores`) - meaning a second
  # call on the same file would ignore the ncores the caller passed and write
  # back to the path recorded during the first call. Keep private copies and
  # restore them after loading.
  #
  # This is also why the response seed is called `response_seed` and not `seed`:
  # `seed` is the stimulus seed stored in the file, so an argument of that name
  # would overwrite it here and then be written back, corrupting the record of
  # how the stimuli were generated.
  .args <- list(rdata = rdata, iter = iter, ncores = ncores,
    response_seed = response_seed, save_rdata = save_rdata, baseimage = baseimage
  )

  # Load parameter file (created when generating stimuli)
  load(rdata)

  rdata <- .args$rdata
  iter <- .args$iter
  ncores <- .args$ncores
  response_seed <- .args$response_seed
  save_rdata <- .args$save_rdata
  baseimage <- .args$baseimage

  # Shared parameters need only the first base. Inline arguments avoid locals
  # leaking into the frame that is re-saved below.
  write("Building the reference from the saved noise, please wait...", stdout())
  stimuli <- referenceNoise(environment(), names(stimuli_params)[1], ncores)

  # Simulate random responding in 2IFC task with ntrials trials across iter iterations
  write("Computing reference distribution, please wait...", stdout())

  # Seed the *responses* only. A response_seed replaces the stimulus stream
  # entirely; without one, the draws the generator spent on the parameters are
  # replayed so the responses continue from where they always did. Handing the
  # stimulus seed a different value instead would describe stimuli the
  # participants never saw.
  seedResponseStream(environment(), names(stimuli_params)[1], response_seed)

  if (iter < 10000) {
    warning("You should set iter >= 10000 for InfoVal statistic to be reliable")
  }

  # Initialize progressbar (dplyr::progress_estimated() is deprecated)
  pb <- txtProgressBar(min = 0, max = iter, style = 3)

  # Run simulation
  reference_norms <- vector(length = iter)

  for (i in 1:iter) {
    setTxtProgressBar(pb, i)

    # Generate random responses for this iteration.
    # This is exactly what the deprecated purrr::rbernoulli(n, p) did
    # internally. It is spelled out rather than swapped for rbinom() on
    # purpose: rbinom() consumes the random stream differently, so it would
    # silently change every reference distribution - and therefore every
    # infoVal - computed from a given seed.
    responses <- ((runif(n_trials) > 0.5) * 2) - 1

    # Compute classification image for this iteration
    ci <- (stimuli %*% as.matrix(responses)) / ncol(stimuli)

    # Save norm for this iteration
    reference_norms[i] <- norm(ci, "f")
  }

  close(pb)

  if (save_rdata) {

    # Save reference norms to rdata file
    write("\nSaving simulated reference distribution to rdata file...", stdout())

    # Provenance belongs to the saved norms; function arguments and scratch state do not.
    reference_norms_seed <- response_seed # nolint: object_usage_linter.
    reference_norms_source <- "saved_noise" # nolint: object_usage_linter.
    reference_norms_fingerprint <- referenceFingerprint(reference_norms) # nolint: object_usage_linter.
    outfile <- rdata
    internals <- c("stimuli", "responses", "pb", "ci", "i", ".args",
      "rdata", "iter", "ncores", "response_seed", "save_rdata",
      "outfile", "internals", "reference_selection", "baseimage"
    )
    save(list = setdiff(ls(all.names = TRUE), internals), file = outfile,
      envir = environment()
    )

  }

  invisible(reference_norms)

}
