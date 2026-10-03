#' Computes Informational Value
#'
#' Computes the Informational Value of a single CI from a 2IFC task.
#'
#' The Informational Value is a z-score for the signal in a classification image: the higher it
#' is, the more signal. A cut-off such as z = 1.96 selects classification images with significant
#' signal at alpha = 0.05.
#'
#' It is computed against a reference distribution: classification images simulated from random
#' responses under the same task parameters as the real data. The Informational Value expresses
#' how unlikely the observed CI is under the null hypothesis that the responses were random.
#'
#' Simulating the reference distribution takes a long time. It is simulated whenever the
#' \code{rdata} file does not already hold one, and then stored there for reuse. If the file
#' cannot be written, the simulated reference is used for this call only: a note names the file
#' (and, for independent base images, the base) that could not be updated, and the next call
#' simulates it again. An archive on read-only media can therefore still be scored, at the cost of
#' simulating each time.
#'
#' @section Matching the reference to the classification image:
#' The reference must be built from the same stimuli as the classification image (Brinkman et
#' al., 2019, Part I). By default it uses every saved stimulus, with one response each. Matching
#' it to your design is your responsibility, because the default cannot know which trials you
#' dropped. A CI built from fewer stimuli (after removing missed trials, say, or from part of the
#' set) has a larger norm under random responding, and the default reference inflates its
#' InfoVal. At 512 pixels, in the stimulus sets measured, pure-noise CIs had a median InfoVal of
#' 0.09 to 0.11 with 1\% of trials missing and 0.51 to 0.58 with 5\%
#' (\url{https://github.com/rdotsch/rcicr/blob/main/analyses/infoval-design-mismatch.md}).
#'
#' Pass the stimuli the CI was built from as \code{reference_stimuli}. \code{\link{generateCI}}
#' records them on its result, so
#' \code{computeInfoVal2IFC(ci, rdata, reference_stimuli = attr(ci, "trial_design")$stimuli)}
#' does it, and a message says when a CI is scored against a reference over different stimuli.
#' The message never changes the number returned. To score many classification images, one per
#' participant say, use \code{\link{batchComputeInfoVal2IFC}}, which takes \code{reference_stimuli}
#' per image and computes each distinct reference once.
#'
#' No reference is defined for a CI that averages repeated presentations of a stimulus or several
#' participants: every reference here assumes one response per stimulus from one responder. For
#' participants who each saw every stimulus once, compute the InfoVal of each participant's own
#' CI instead.
#'
#' @section Masked classification images:
#' A classification image computed with \code{mask} (see \code{\link{generateCI}}) holds
#' \code{NA} in its masked pixels. Its InfoVal is computed over the unmasked pixels only, against
#' a reference built over the same pixels from the same stimuli: Brinkman et al.'s (2019)
#' statistic for the region analysed. The mask is read from the classification image itself.
#' Each mask gets its own reference, simulated once and stored in the \code{rdata} file apart
#' from the default one; \code{generateReferenceDistribution2IFC(mask = )} stores one in
#' advance. An InfoVal over part of the image is not comparable with one over the whole image.
#'
#' For the method, see Brinkman, L., Goffin, S., van de Schoot, R., van Haren, N. E. M.,
#' Dotsch, R., & Aarts, H. (2019). Quantifying the informational value of classification
#' images. \emph{Behavior Research Methods}, \emph{51}, 2059-2073.
#' \doi{10.3758/s13428-019-01232-2}
#'
#' @export
#' @importFrom stats mad median
#' @param target_ci A classification image, as the list returned by \code{\link{generateCI}}.
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus and, once computed, the reference distribution (see \code{\link{generateReferenceDistribution2IFC}}).
#' @param iter Number of simulated classification images in the reference distribution. Used only when the reference distribution has to be simulated.
#' @param force_gen_ref_dist Boolean: simulate the reference distribution again even if the \code{rdata} file already holds one (default: \code{FALSE}).
#' @param response_seed Optional seed for the simulated random responses behind the reference
#' distribution. The default, \code{NULL}, uses the reference distribution stored in the
#' \code{rdata} file, or simulates the reproducible default one described under Reproducibility
#' in \code{\link{generateReferenceDistribution2IFC}}, which needs the stimulus seed saved in the
#' file. For a file without one, pass a number. A number forces a fresh reference
#' distribution from an independent draw; use it to check how much Monte Carlo error \code{iter}
#' leaves in the Informational Value. That result is \emph{not} written back to the \code{rdata}
#' file, so a one-off check cannot change the number every later analysis of the stimulus set
#' reports.
#' @param baseimage The base-image label \code{target_ci} was computed for, the same one passed
#' to \code{generateCI()}. Required when the base images have different noise parameters; each
#' base then gets its own stored reference distribution, and a shared one left by an older version
#' is ignored. With a single base image, or base images sharing one parameter set, leave it at
#' \code{NULL}.
#' @param reference_stimuli Optional stimulus numbers the classification image was built from,
#' each once, when that is not every saved stimulus. The reference is then built over exactly
#' those stimuli and stored in the \code{rdata} file apart from the default one, as described in
#' \code{\link{generateReferenceDistribution2IFC}}. The default, \code{NULL}, uses every saved
#' stimulus, as does passing all of them, with the same result in every case. A subset too small
#' for random responses to give distinct norms (one stimulus, and usually two) is refused, since
#' its MAD is 0 and the InfoVal would not be a number. See "Matching the reference to the classification image" below.
#' @param reference_method \code{"gram"} (the default) or \code{"images"}: how the reference is
#' computed, when it has to be. \code{"images"} reproduces rcicr 1.5.0 and earlier bit for bit; a
#' reference already stored in \code{rdata} is reused whichever is given. See "Reference method" in
#' \code{\link{generateReferenceDistribution2IFC}}.
#' @return The Informational Value, a z-score.
#' @seealso \code{vignette("reverse-correlation-walkthrough")}, section "Is there actually signal?"; \code{vignette("recipes")}, "Matching the InfoVal reference to the design", for dropped trials, masks, batches and group averages.
#' @examples
#' # a synthetic square grayscale image stands in for a real base face photo
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#'
#' stimulus_path <- tempfile("stimuli")
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
computeInfoVal2IFC <- function(target_ci, rdata, iter = 10000, force_gen_ref_dist = FALSE, response_seed = NULL, baseimage = NULL, reference_stimuli = NULL, reference_method = c("gram", "images")) {

  reference_method <- match.arg(reference_method)
  selection <- selectReferenceBase(rdata, baseimage)
  report <- function(n_trials, reference_stimuli) reportTrialDesign(target_ci, n_trials, reference_stimuli)
  subset <- subsetReferenceFor(rdata, selection, reference_stimuli, ciMask(target_ci))
  reference <- if (!is.null(subset)) {
    subsetReference(rdata, iter, force_gen_ref_dist, response_seed, subset$source, selection,
                    subset$reference_stimuli, reference_method, report, subset$masked)
  } else if (selection$independent) {
    baseReference(rdata, iter, force_gen_ref_dist, response_seed, selection, reference_method,
                  report)
  } else {
    sharedReference(rdata, iter, force_gen_ref_dist, response_seed, reference_method, report)
  }

  cinorm <- ciNorm(target_ci)
  info_val <- (cinorm - reference$median) / reference$mad
  write(paste0('Informational value: z = ', info_val, ' (', reference$note, 'ci norm = ', cinorm,
               '; reference median = ', reference$median, '; MAD = ', reference$mad,
               '; iterations = ', reference$iter, ')'), stdout())

  return(info_val)
}

# The default reference of a file whose bases share one parameter matrix, as
# list(median, mad, iter, note). report(n_trials, NULL) runs once the file is
# loaded, before anything is simulated.
sharedReference <- function(rdata, iter, force_gen_ref_dist, response_seed, reference_method, report) {
  # Old files may contain argument names, including a stale rdata path.
  .args <- captureArgs(environment())
  loadRdata(rdata, environment())
  list2env(.args, envir = environment())

  report(get0('n_trials', envir = environment(), inherits = FALSE), NULL)

  if (!is.null(response_seed)) force_gen_ref_dist <- TRUE
  cached_reference <- if (exists('reference_norms', envir = environment(), inherits = FALSE)) {
    list(
      norms = reference_norms,
      response_seed = get0('reference_norms_seed', envir = environment(), inherits = FALSE),
      source = get0('reference_norms_source', envir = environment(), inherits = FALSE),
      fingerprint = get0('reference_norms_fingerprint', envir = environment(), inherits = FALSE)
    )
  } else {
    NULL
  }

  seed <- get0('seed', envir = environment(), inherits = FALSE)
  if (!force_gen_ref_dist && is.null(cached_reference)) requireStimulusSeed(seed, rdata)

  reference_norms <- resolveReferenceNorms(cached_reference, rdata, iter,
                                           force_gen_ref_dist, response_seed,
                                           seedless = is.null(seed),
                                           reference_method = reference_method)
  referenceSummary(reference_norms, 'in the stimulus file', 'in the stimulus file', '')
}
