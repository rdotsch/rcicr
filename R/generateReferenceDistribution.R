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
#' @return Nothing.

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

      # Compute classication image for this iteration
      ci <- (as.matrix(stimuli) %*% as.matrix(responses)) / ncol(stimuli)

      # Save norm for this iteration
      reference_norms[i] <- norm(ci, "f")
  }

  # Save reference norms to rdata file
  write("\nSaving simulated reference distribution to rdata file...", stdout())
  rm(stimuli, responses, pb, iter, ci)
  save(list=ls(all.names=TRUE), file=rdata, envir=environment())

}
