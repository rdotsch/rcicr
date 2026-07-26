#' Generates reference distribution
#'
#' Generates reference distribution of norms for a particular set of task parameters.
#'
#' In order to compute the Informational Value metric. Saves its results in the supplied rdata file for later reuse.
#'
#' @export
#' @importFrom stats runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param iter Number of iterations for the simulation (i.e., the number of norms generated with classification images based on random responding).
#' @param ncores Number of CPU cores to use when re-generating the stimuli (default: detectCores()-1).
#' @return Nothing. The reference distribution (\code{reference_norms}) is added to the supplied
#' \code{rdata} file, so a later call to \code{\link{computeInfoVal2IFC}} using the same file can
#' reuse it instead of re-simulating.
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
generateReferenceDistribution2IFC <- function(rdata, iter=10000, ncores=parallel::detectCores()-1) {

  # Load parameter file (created when generating stimuli)
  load(rdata)

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

  # Save reference norms to rdata file
  write("\nSaving simulated reference distribution to rdata file...", stdout())
  close(pb)
  rm(stimuli, responses, pb, iter, ci)
  save(list=ls(all.names=TRUE), file=rdata, envir=environment())

}
