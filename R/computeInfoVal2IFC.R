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
#' @importFrom dplyr filter count summarise %>%
#' @import yesno
#' @param target_ci A classification image object (list-type) as returned by generateCI
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli and possibly its corresponding reference distribution generated with generateReferenceDistribution().
#' @param iter Number of iterations for the simulation of the reference distribution (only used if reference distribution is not already pre-generated and present in rdata file)
#' @param force_gen_ref_dist Boolean specifying whether to override the default behavior to use pre-computed values for the reference distribution for specific task parameters and instead force to recompute the reference distribution (default: FALSE).
#' @param response_seed Optional seed for the simulated random responses used to build the
#' reference distribution. The default (\code{NULL}) uses the reference distribution stored in
#' the \code{rdata} file, or generates the reproducible default one described under
#' Reproducibility in \code{\link{generateReferenceDistribution2IFC}}. Supplying a number
#' forces a fresh reference distribution to be simulated from an independent draw, which is
#' how you check how much Monte Carlo error \code{iter} leaves in this Informational Value.
#' The result is deliberately \emph{not} written back to the \code{rdata} file, so a one-off
#' check cannot change the number every later analysis of that stimulus set reports.
#' @param baseimage Saved base-image label used to generate \code{target_ci}. Required when
#' saved base images have different noise parameters; use the same label passed to
#' \code{generateCI()}. With a single base or identical parameter matrices, the default
#' \code{NULL} retains the shared reference behavior. Independent-base caches are separate
#' for each label; old unscoped reference norms are ignored for those files.
#' @return Informational value (z-score)
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
#' # compute (and cache in rdata_file) a reference distribution; iter is kept
#' # tiny here for a fast example, in practice use iter >= 10000.
#' suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3, ncores = 1))
#'
#' responses <- sample(c(1, -1), 6, replace = TRUE)
#' target_ci <- generateCI(
#'   stimuli = 1:6, responses = responses, baseimage = "face",
#'   rdata = rdata_file, save_as_png = FALSE
#' )
#'
#' computeInfoVal2IFC(target_ci = target_ci, rdata = rdata_file)
computeInfoVal2IFC <- function(target_ci, rdata, iter = 10000, force_gen_ref_dist = FALSE, response_seed = NULL, baseimage = NULL) {

  reference_selection <- selectReferenceBase(rdata, baseimage)
  if (reference_selection$independent) {
    return(computeBaseInfoVal(target_ci, rdata, iter, force_gen_ref_dist,
      response_seed, reference_selection))
  }

  # RD: To supress notes from R CMD CHECK, but thise should not be necessary -- debug
  ref_seed <- NA
  ref_img_size <- NA
  ref_n_trials <- NA

  # Load parameter file (created when generating stimuli)
  # load() assigns into this frame, so an .Rdata file written by an older
  # generateReferenceDistribution2IFC() - which saved its own `rdata` argument
  # into the file - would overwrite the path we were called with. This function
  # still uses `rdata` further down (to regenerate and re-load), so it would
  # then operate on whatever path that file happened to record.
  #
  # Keep every argument rather than the three that were known to collide: the
  # hazard has bitten from the .Rdata side twice now (a `sigma` field added to
  # the file captured generateCI()'s z-map argument; fixed in #146), so the
  # guard has to hold for fields that do not exist yet. `target_ci` is the one
  # that would hurt most here - it is read at the very end to compute the CI
  # norm, so a file carrying that name would silently score somebody else's
  # classification image and return a plausible number rather than an error.
  .args <- captureArgs(environment())
  load(rdata)
  list2env(.args, envir = environment())

  # Asking for a specific response seed is asking for a specific reference
  # distribution, so it has to imply regeneration. Without this the argument
  # would be silently ignored on every file that already has reference_norms -
  # which is every file after the first call, since generating one writes it
  # back. That is the same shape of bug as force_gen_ref_dist being ignored
  # (see the comment further down) and the documented-but-unapplied `mask`
  # argument of plotZmap(): accepted, documented, and doing nothing.
  if (!is.null(response_seed)) {
    force_gen_ref_dist <- TRUE
  }

  # Check whether reference norms are present or can be looked up from table. If not, re-generate.
  if (!force_gen_ref_dist && !exists("reference_norms", envir = environment(), inherits = FALSE)) {

    # Pre-computed reference distribution table (TODO: read from external file).
    #
    # THIS TABLE IS EMPTY, AND THAT IS CORRECT. The four rows below were
    # measured under the pre-2018 infoVal formula and were commented out by
    # commit 01e547e, which adopted the Euclidean norm and scaling factor k from
    # the erratum to Schmitz et al. (2019). That change redefined the norms these
    # medians and MADs summarise, so reusing them would silently score every CI
    # against a null from the wrong formula. Emptying the table was the right
    # call; the numbers were never re-measured.
    #
    # The consequence is that every lookup below misses and the reference
    # distribution is always regenerated -- correct, just slow. The matching and
    # prompting machinery is kept rather than deleted so that repopulating the
    # table is a matter of measuring four numbers, not rebuilding the feature.
    # To repopulate: run generateReferenceDistribution2IFC() at each parameter
    # combination and record median(reference_norms) and mad(reference_norms).
    ref_lookup <- tribble(
      ~ref_seed, ~ref_img_size, ~ref_iter, ~ref_n_trials, ~ref_median, ~ref_mad,
      #  1,         512,           10000,     100,           1097.7394,     52.54232,
      #  1,         512,           10000,     300,           634.0318,      30.51781,
      #  1,         512,           10000,     500,           490.4709,      23.71276,
      #  1,         512,           10000,     1000,          347.2960,      16.64761
    )

    # Check whether we have a perfect match
    ref_values <- ref_lookup %>%
      filter(ref_seed == seed, ref_img_size == img_size, ref_n_trials == n_trials, ref_iter == iter)

    if (ref_values %>% count() == 1) {
      # We have a match, use the values
      write("Pre-computed reference values matching your exact parameters found.", stdout())

      ref_median <- ref_values$ref_median
      ref_mad <- ref_values$ref_mad
      ref_iter <- ref_values$ref_iter

    } else {
      # Check whether at least seed, img_size, and n_trials match
      ref_values <- ref_lookup %>%
        filter(ref_seed == seed, ref_img_size == img_size, ref_n_trials == n_trials)

      if (ref_values %>% count() > 0) {
        write("I found pre-computed reference values that matched seed, image size, and number of trials, but not the number of reference distribution iterations.", stdout())
        max_ref_iter <- as.numeric(ref_values %>% summarise(max(ref_iter)))

        # Only ask when there is somebody to answer. In a non-interactive
        # session -- a batch script, knitr, R CMD check -- there is no way to
        # respond, so decline and regenerate rather than silently substituting a
        # reference distribution built with a different number of iterations than
        # the caller asked for. Regenerating is slower but is what was requested.
        if (interactive()) {
          user_response <- yesno::yesno(paste0("I did find pre-computed values for ", max_ref_iter, " iterations matching all other parameters. Do you want to use those instead?"))
        } else {
          write(paste0("Not running interactively, so I cannot ask -- regenerating the reference distribution with the requested ", iter, " iterations rather than reusing the pre-computed ", max_ref_iter, ". Pass the pre-computed value explicitly if you want it."), stdout())
          user_response <- FALSE
        }

        if (user_response) {
          write(paste0("Using pre-computed reference values for ", max_ref_iter, " instead of ", iter, " iterations."), stdout())
          ref_values <- ref_lookup %>%
            filter(ref_seed == seed, ref_img_size == img_size, ref_n_trials == n_trials, ref_iter == max_ref_iter)

          ref_median <- ref_values$ref_median
          ref_mad <- ref_values$ref_mad
          ref_iter <- ref_values$ref_iter
        }
      }
    }
  }

  if (!exists("ref_median", envir = environment(), inherits = FALSE)) {

    # Regenerate when there is no cached distribution, or when the caller
    # explicitly asked for one. force_gen_ref_dist previously only skipped the
    # lookup-table branch above and never reached here, so it was silently
    # ignored whenever reference_norms already existed in the .Rdata file.
    if (force_gen_ref_dist || !exists("reference_norms", envir = environment(), inherits = FALSE)) {

      # Reference norms not present in rdata file (or regeneration forced).
      #
      # A caller-supplied response_seed is a one-off check, not a redefinition
      # of this stimulus set's null, so it is never cached: saving it would
      # silently change what every later InfoVal computed from this file means,
      # for anyone using it, with nothing in the call to say so.
      cache_ref_dist <- is.null(response_seed)

      reference_norms <- generateReferenceDistribution2IFC(
        rdata, iter = iter, response_seed = response_seed, save_rdata = cache_ref_dist
      )

      if (cache_ref_dist) {

        # Re-load rdata file, to pick up the reference_norms just written to it.
        # This load carries the same hazard as the one above and needs the same
        # restore: `target_ci` is read after this point (line ~230, for the CI
        # norm), so without it a file containing that name would replace the
        # caller's classification image between here and the computation of it.
        load(rdata)
        list2env(.args, envir = environment())

        # NB: write() defaults to file = "data", so omitting stdout() here did
        # not print this message - it silently created a file called "data" in
        # the working directory. Every other write() in the package passes
        # stdout(); this one was missed.
        write("Note that now that this simulated reference distribution has been saved to the .Rdata file, the next time you call computeInfoVal2IFC(), it will not need to be computed again.", stdout())

      } else {

        write(paste0("Reference distribution simulated with response_seed = ", response_seed,
                ". This is an independent draw of the null, not the reference distribution ",
                "stored in the .Rdata file, and it has deliberately not been saved there."
              ), stdout())

      }

    } else {

      write("Using reference distribution found in rdata file.", stdout())

    }

    # Compute reference values
    ref_median <- median(reference_norms)
    ref_mad <- mad(reference_norms)
    ref_iter <- length(reference_norms)
  }

  # Compute informational value metric
  cinorm <- norm(matrix(target_ci[["ci"]]), "f")
  infoVal <- (cinorm - ref_median) / (ref_mad)

  write(paste0("Informational value: z = ", infoVal, " (ci norm = ", cinorm, "; reference median = ", ref_median, "; MAD = ", ref_mad, "; iterations = ", ref_iter, ")"), stdout())

  return(infoVal)

}
