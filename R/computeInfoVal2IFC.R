#' Computes Informational Value
#'
#' Computes Informational Value for a single CI in a 2IFC task.
#'
#' The Informational Value metric can be considered as a z-score that quantifies the signal
#' present in a classification image. The higher the Informational Value, the more signal. It is
#' possible to use a cut-off such as z = 1.96 to select classification images with significant
#' signal under alpha = 0.05.
#'
#' Informational Value is computed by simulating random responding under identical task parameters to
#' an empirical dataset (called the reference distribution). The metric quantifies how unlikely it is
#' to observe these data under the null-hypothesis that there is no signal (i.e., that there is only random responding).
#'
#' The simulation to compute the reference distribution takes a long time, and is only run locally when
#' pre-computed values for the reference distribution matching the stimulus set in the .Rdata file have
#' not been supplied by the rcicr package.
#'
#' For more information see Brinkman, Goffin, Aarts, van Haren, & Dotsch (in prep).
#'
#' @export
#' @importFrom stats mad median
#' @importFrom tibble tribble
#' @importFrom dplyr filter
#' @import yesno
#' @param target_ci A classification image object (list-type) as returned by generateCI
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli and possibly its corresponding reference distribution generated with generateReferenceDistribution().
#' @param iter Number of iterations for the simulation of the reference distribution (only used if reference distribution is not already pre-generated and present in rdata file)
#' @param force_gen_ref_dist Boolean specifying whether to override the default behavior to use pre-computed values for the reference distribution for specific task parameters and instead force to recompute the reference distribution (default: FALSE).
#' @return Informational value (z-score)

computeInfoVal2IFC <- function(target_ci, rdata, iter = 10000, force_gen_ref_dist = FALSE) {

  # RD: To supress notes from R CMD CHECK, but thise should not be necessary -- debug
  ref_seed <- NA
  ref_img_size <- NA
  ref_n_trials <- NA

  # Load parameter file (created when generating stimuli)
  load(rdata)

  # Check whether reference norms are present or can be looked up from table. If not, re-generate.
  if (!force_gen_ref_dist & !exists("reference_norms", envir=environment(), inherits=FALSE)) {

    # Pre-computed reference distribution table (TODO: read from external file)
    ref_lookup <- tribble(
      ~ref_seed, ~ref_img_size, ~ref_iter, ~ref_n_trials, ~ref_median,   ~ref_mad,
    #  1,         512,           10000,     100,           1097.7394,     52.54232,
    #  1,         512,           10000,     300,           634.0318,      30.51781,
    #  1,         512,           10000,     500,           490.4709,      23.71276,
    #  1,         512,           10000,     1000,          347.2960,      16.64761
    )

    # Check whether we have a perfect match
    ref_values <- ref_lookup %>%
      filter(ref_seed==seed, ref_img_size==img_size, ref_n_trials==n_trials, ref_iter==iter)

    if (ref_values %>% count() == 1) {
      # We have a match, use the values
      write("Pre-computed reference values matching your exact parameters found.", stdout())

      ref_median <- ref_values$ref_median
      ref_mad    <- ref_values$ref_mad
      ref_iter   <- ref_values$ref_iter

    } else {
      # Check whether at least seed, img_size, and n_trials match
      ref_values <- ref_lookup %>%
        filter(ref_seed==seed, ref_img_size==img_size, ref_n_trials==n_trials)

      if (ref_values %>% count() > 0) {
        write("I found pre-computed reference values that matched seed, image size, and number of trials, but not the number of reference distribution iterations.", stdout())
        max_ref_iter <- as.numeric(ref_values %>% summarise(max(ref_iter)))
        user_response <- yesno::yesno(paste0("I did find pre-computed values for ", max_ref_iter, " iterations matching all other parameters. Do you want to use those instead?"))

        if (user_response) {
          write(paste0("Using pre-computed reference values for ", max_ref_iter, " instead of ", iter, " iterations."), stdout())
          ref_values <- ref_lookup %>%
            filter(ref_seed==seed, ref_img_size==img_size, ref_n_trials==n_trials, ref_iter==max_ref_iter)

          ref_median <- ref_values$ref_median
          ref_mad    <- ref_values$ref_mad
          ref_iter   <- ref_values$ref_iter
        }
      }
    }
  }

  if (!exists("ref_median", envir=environment(), inherits=FALSE)) {

    if (!exists("reference_norms", envir=environment(), inherits=FALSE)) {

      # Reference norms not present in rdata file, re-generate
      generateReferenceDistribution2IFC(rdata, iter=iter)

      # Re-load rdata file
      load(rdata)

      write("Note that now that this simulated reference distribution has been saved to the .Rdata file, the next time you call computeInfoVal2IFC(), it will not need to be computed again.")

    } else {

      write("Using reference distribution found in rdata file.", stdout())

    }

    # Compute reference values
    ref_median <- median(reference_norms)
    ref_mad    <- mad(reference_norms)
    ref_iter   <- length(reference_norms)
  }

  # Compute informational value metric
  cinorm <- norm(matrix(target_ci[["ci"]]), "f")
  infoVal <- (cinorm - ref_median ) / (ref_mad)

  write( paste0("Informational value: z = ", infoVal, " (ci norm = ", cinorm,"; reference median = ", ref_median, "; MAD = ", ref_mad, "; iterations = ", ref_iter,  ")"), stdout() )

  return(infoVal)

}
