#' Generates reference distribution
#'
#' Generates reference distribution of norms for a particular set of task parameters.
#'
#' In order to compute the Informational Value metric. Saves its results in the supplied rdata file for later reuse.
#'
#' @section Reproducibility:
#' With the default \code{response_seed = NULL}, the reference distribution is determined by
#' the stimulus \code{.Rdata} file alone. It does not depend on the ambient random number
#' state of the calling session, and it does not depend on \code{ncores}. Two researchers
#' who compute InfoVal from the same stimulus file therefore get the same number, and the
#' same reference distribution, on different machines and in different sessions.
#'
#' This is a guarantee, not a coincidence, and it is relied upon: the function re-generates
#' the stimuli through \code{\link{generateStimuli2IFC}}, whose internal \code{set.seed()}
#' call uses the seed stored in the \code{.Rdata} file and lands before the random responses
#' below are drawn.
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
#' @param ncores Number of CPU cores to use when re-generating the stimuli (default: \code{detectCores()-1}; 2 under \code{R CMD check}, per CRAN policy).
#' @param response_seed Optional seed for the simulated random responses. The default
#' (\code{NULL}) draws them from the state left by the stimulus re-generation, which is the
#' reproducible behaviour described under Reproducibility. Supply a number to obtain an
#' independent draw of the null from the same stimuli.
#' @param save_rdata Boolean specifying whether the reference distribution should be written
#' back into the \code{rdata} file (default \code{TRUE}). Set to \code{FALSE} to compute a
#' distribution without changing what later calls to \code{\link{computeInfoVal2IFC}} will
#' use -- worth doing whenever \code{response_seed} is set, so a one-off null does not become
#' the file's permanent reference.
#' @return The reference distribution, invisibly, as a numeric vector of \code{iter} norms.
#' Unless \code{save_rdata = FALSE}, it is also added to the supplied \code{rdata} file as
#' \code{reference_norms} (alongside \code{reference_norms_seed}, recording the
#' \code{response_seed} it was generated with), so a later call to
#' \code{\link{computeInfoVal2IFC}} using the same file can reuse it instead of re-simulating.
#' @examples
#' \donttest{
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
#' # Run from a temp working directory: this function re-generates stimuli via
#' # generateStimuli2IFC() without forwarding a stimulus_path, so it always
#' # creates a ./stimuli directory relative to the working directory.
#' withr::with_dir(tempdir(), {
#'   suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3, ncores = 1))
#' })
#' }
generateReferenceDistribution2IFC <- function(rdata, iter=10000, ncores=default_ncores(), response_seed=NULL, save_rdata=TRUE) {

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
                response_seed = response_seed, save_rdata = save_rdata)

  # Load parameter file (created when generating stimuli)
  load(rdata)

  rdata         <- .args$rdata
  iter          <- .args$iter
  ncores        <- .args$ncores
  response_seed <- .args$response_seed
  save_rdata    <- .args$save_rdata

  # Recover the noise-basis parameters used for the real stimuli. These were
  # not saved before this version, so .Rdata files written by older rcicr lack
  # them. Falling back to the defaults silently would rebuild the reference
  # distribution on a *different* noise basis than participants actually saw,
  # producing a wrong infoVal - so warn loudly rather than guessing quietly.
  if (!exists('nscales', envir=environment(), inherits=FALSE)) {
    nscales <- 5
    warning(paste0('This .Rdata file was written by a version of rcicr that ',
                   'did not save `nscales`, so the default (5) is assumed for ',
                   'the reference distribution. If the stimuli were generated ',
                   'with a different nscales, the resulting infoVal will be ',
                   'wrong - regenerate the stimulus set with this version of ',
                   'rcicr to fix this.'))
  }
  if (!exists('sigma', envir=environment(), inherits=FALSE)) {
    sigma <- 25
  }

  # Re-generate stimuli based on rdata parameters in matrix form
  write("Re-generating stimuli based on rdata file, please wait...", stdout())
  stimuli <- generateStimuli2IFC(base_face_files, n_trials, img_size, seed=seed, noise_type=noise_type, nscales=nscales, sigma=sigma, ncores=ncores, return_as_dataframe=TRUE, save_as_png=FALSE, save_rdata=FALSE)

  # Simulate random responding in 2IFC task with ntrials trials across iter iterations
  write("Computing reference distribution, please wait...", stdout())

  # Seed the *responses* only, and only when asked. This deliberately sits after
  # the stimulus re-generation above rather than being passed into it: handing a
  # different seed to generateStimuli2IFC() would rebuild a different stimulus
  # set, so the null would describe stimuli the participants never saw. What
  # varies here is the simulated responding, on the same stimuli.
  #
  # NULL means no set.seed() call at all, leaving the stream exactly as the
  # stimulus re-generation left it - so the default path is byte-identical to
  # what previous versions produced, not merely equivalent.
  if (!is.null(response_seed)) {
    set.seed(response_seed)
  }

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
      ci <- (as.matrix(stimuli) %*% as.matrix(responses)) / ncol(stimuli)

      # Save norm for this iteration
      reference_norms[i] <- norm(ci, "f")
  }

  close(pb)

  if (save_rdata) {

    # Save reference norms to rdata file
    write("\nSaving simulated reference distribution to rdata file...", stdout())

    # Record which seed produced these norms, so a file carrying a deliberately
    # varied null is distinguishable from one carrying the default. NULL records
    # the default stream. Files written before this version simply lack the
    # field, which is why every read of it must be guarded with exists().
    reference_norms_seed <- response_seed

    # Save everything that came from (or belongs in) the stimulus file, but none
    # of this function's own arguments or scratch variables. Writing `rdata` and
    # `ncores` back into the file is what causes the clobbering described at the
    # top of this function, so they are excluded at the source rather than only
    # worked around on read. `response_seed` and `save_rdata` are excluded for
    # the same reason - `reference_norms_seed` is the field that records the
    # seed, and it is a description of the norms rather than an input.
    outfile <- rdata
    internals <- c("stimuli", "responses", "pb", "ci", "i", ".args",
                   "rdata", "iter", "ncores", "response_seed", "save_rdata",
                   "outfile", "internals")
    save(list=setdiff(ls(all.names=TRUE), internals), file=outfile,
         envir=environment())

  }

  invisible(reference_norms)

}
