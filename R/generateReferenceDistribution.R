#' Generates reference distribution
#'
#' Generates reference distribution of norms for a particular set of task parameters.
#'
#' In order to compute the Informational Value metric. Saves its results in the supplied rdata file for later reuse.
#'
#' @export
#' @importFrom purrr rbernoulli
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param iter Number of iterations for the simulation (i.e., the number of norms generated with classification images based on random responding).
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
#' # iter is kept tiny here for a fast example; in practice use iter >= 10000
#' # (run from a temp working directory: generateReferenceDistribution2IFC()
#' # internally re-generates stimuli via generateStimuli2IFC() without
#' # forwarding a stimulus_path, so it always creates a ./stimuli directory
#' # relative to the working directory)
#' withr::with_dir(tempdir(), {
#'   suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3))
#' })
#' }
generateReferenceDistribution2IFC <- function(rdata, iter=10000) {

  # Load parameter file (created when generating stimuli)
  load(rdata)

  # Re-generate stimuli based on rdata parameters in matrix form
  write("Re-generating stimuli based on rdata file, please wait...", stdout())
  stimuli <- generateStimuli2IFC(base_face_files, n_trials, img_size, seed=seed, noise_type=noise_type,ncores=parallel::detectCores()-1, return_as_dataframe=TRUE, save_as_png=FALSE, save_rdata=FALSE)

  # Simulate random responding in 2IFC task with ntrials trials across iter iterations
  write("Computing reference distribution, please wait...", stdout())

  if (iter < 10000) {
    warning("You should set iter >= 10000 for InfoVal statistic to be reliable")
  }

  # Initialize progressbar
  pb <- progress_estimated(iter)

  # Run simulation
  reference_norms <- vector(length = iter)

  for (i in 1:iter) {
      pb$tick()$print()

      # Generate random responses for this iteration
      responses <- (purrr::rbernoulli(n_trials, p=0.5) * 2) - 1

      # Compute classification image for this iteration
      ci <- (as.matrix(stimuli) %*% as.matrix(responses)) / ncol(stimuli)

      # Save norm for this iteration
      reference_norms[i] <- norm(ci, "f")
  }

  # Save reference norms to rdata file
  write("\nSaving simulated reference distribution to rdata file...", stdout())
  rm(stimuli, responses, pb, iter, ci)
  save(list=ls(all.names=TRUE), file=rdata, envir=environment())

}
