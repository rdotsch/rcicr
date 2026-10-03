#' Generates reference distribution
#'
#' Generates the reference distribution of norms for a stimulus set.
#'
#' \code{\link{computeInfoVal2IFC}} scores a classification image against this distribution. By
#' default the result is saved in the \code{rdata} file for later reuse.
#'
#' @section Reproducibility:
#' With the default \code{response_seed = NULL}, the reference distribution depends only on the
#' stimulus \code{.Rdata} file: not on the current random state, not on \code{ncores}, and,
#' for files that record the RNG kind, not on the session's \code{\link{RNGkind}}. Two
#' researchers computing InfoVal from the same stimulus file get the same reference distribution,
#' and the same number, on any machine and in any session.
#'
#' This needs the stimulus seed saved in the file. A file without one (made with
#' \code{generateStimuli2IFC(seed = NULL)}, or with the field removed) has no stream to continue,
#' so the default stops with an error; pass a \code{response_seed} instead.
#'
#' Stimulus files record the \code{\link{RNGkind}} the stimuli were drawn under, as
#' \code{rng_kind}, and the simulated responses are drawn under it, with or without a
#' \code{response_seed}. When it differs from the session's kind, the caller's random stream is
#' restored afterwards, kind and position, so the session is not left on the file's generator. In
#' the session's own kind the stream is left where the draws end, as it always was. Files written
#' before the kind was recorded draw under the session's kind: a session with a different
#' \code{RNGkind()} draws different responses, and so a different null. To reproduce an earlier
#' reference from such a file, use the RNG kind that built it; see
#' \url{https://github.com/rdotsch/rcicr/issues/315}.
#'
#' The noise is rebuilt from the basis and parameters saved in the stimulus file, so the base
#' images are never reopened and a moved or archived experiment can still be scored. The simulated
#' responses continue the random stream from the saved stimulus seed, after the draws for one
#' parameter matrix, which reproduces the historical default reference.
#'
#' Pass a \code{response_seed} to draw a \emph{different} null from the same stimuli, for
#' instance to check how much Monte Carlo error a given \code{iter} leaves in your InfoVal. This
#' changes only the simulated responses, not the stimuli or the noise basis the null is built on.
#'
#' @section Reference method:
#' \code{reference_method = "gram"}, the default, computes the norms from the stimulus Gram
#' matrix, without calling \code{generateNoiseImage()} or multiplying the noise for every draw.
#' For a small stimulus set the Gram matrix is built from the noise rendered through the sparse
#' basis; otherwise from the basis's cross-product, kept sparse. Either way the basis is built a
#' block of image columns at a time, so no full-size copy of it, or of the noise, is held.
#' \code{"images"} is the calculation of rcicr 1.5.0 and earlier: every noise image rendered, in
#' parallel over \code{ncores}, and multiplied for every draw. It reproduces references from those
#' versions bit for bit.
#'
#' The two agree to rounding. No configuration measured was bit-identical, and none differed by
#' more than a relative 5e-14 in a single norm (reference BLAS; see
#' \url{https://github.com/rdotsch/rcicr/blob/main/analyses/gram-reference-accuracy.md}).
#'
#' Computing a reference with \code{"gram"} prints a message saying so, and naming
#' \code{"images"} as the way to reproduce earlier references. Silence it with
#' \code{suppressMessages()}. The method is never switched for you.
#'
#' A reference already stored in the file is reused as stored, whichever method is asked for. To
#' rebuild one with a particular method, add \code{force_gen_ref_dist = TRUE} in
#' \code{\link{computeInfoVal2IFC}}, or call this function with \code{save_rdata = TRUE}. The
#' method a reference was computed with is stored beside it, as \code{reference_norms_method} or
#' as \code{method} in each \code{reference_norms_by_base} and \code{reference_norms_by_stimuli}
#' entry.
#'
#' @export
#' @importFrom stats runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus.
#' @param iter Number of simulated classification images, each built from random responses; the distribution holds one norm per image.
#' @param ncores Number of CPU cores used to render the saved noise with \code{reference_method = "images"} (default: \code{detectCores() - 1}; 2 under \code{R CMD check}, per CRAN policy). \code{"gram"} does not use it.
#' @param response_seed Optional seed for the simulated random responses. The default,
#' \code{NULL}, continues from the state the stimulus generator left behind, as described under
#' Reproducibility; it needs the stimulus seed saved in the file. A number gives an independent
#' draw of the null from the same stimuli.
#' @param save_rdata Boolean: write the reference distribution into the \code{rdata} file
#' (default \code{TRUE}). With \code{FALSE}, later calls to \code{\link{computeInfoVal2IFC}}
#' keep using what the file already holds. Set it to \code{FALSE} whenever you set
#' \code{response_seed}, so a one-off null does not become the file's permanent reference. The
#' exception is a file without a stimulus seed: there the seeded reference is the one to keep, as
#' the file's reproducible reference.
#' @param baseimage Base-image label, the same key passed to \code{generateCI()}. Required when
#' the base images have different noise parameters. With a single base image, or base images
#' sharing one parameter set, leave it at \code{NULL}.
#' @param reference_stimuli Optional stimulus numbers, each once, to build the reference over
#' instead of every saved stimulus: the stimuli a classification image was built from, when that
#' is not all of them. See "Matching the reference to the classification image" in
#' \code{\link{computeInfoVal2IFC}}. The simulated responses continue the same stream as the
#' default, so the reference is reproducible from the file. Passing every saved stimulus is the
#' same as the default, \code{NULL}.
#' @param reference_method \code{"gram"} (the default) or \code{"images"}: how the reference
#' norms are computed. \code{"images"} reproduces rcicr 1.5.0 and earlier bit for bit. See
#' "Reference method" below.
#' @param mask Optional mask, in any form \code{\link{generateCI}} accepts: a 0/1 matrix or the
#' path to a PNG, black (0) where masked. The reference is then built over the unmasked pixels
#' only, as \code{\link{computeInfoVal2IFC}} needs for a classification image computed with
#' that mask, and stored apart from the default one. The default, \code{NA}, uses every pixel.
#' See "Masked classification images" in \code{\link{computeInfoVal2IFC}}.
#' @section Independent base images:
#' When the base images have different parameter matrices, \code{baseimage} says whose noise to
#' use. The distributions are then stored in \code{reference_norms_by_base}, one entry per base
#' label, each holding \code{norms} and \code{response_seed}. A shared \code{reference_norms}
#' in such a file is neither used nor overwritten. Files whose base images share one parameter
#' matrix keep using \code{reference_norms}.
#' @return The reference distribution, invisibly, as a numeric vector of \code{iter} norms.
#' Unless \code{save_rdata = FALSE}, it is also added to the \code{rdata} file as
#' \code{reference_norms}, with \code{reference_norms_seed} recording the \code{response_seed}
#' it was drawn with. A later \code{\link{computeInfoVal2IFC}} call on the same file then reuses
#' it instead of simulating again. For independent base images it goes in
#' \code{reference_norms_by_base} instead, as described above. A reference over
#' \code{reference_stimuli} goes in \code{reference_norms_by_stimuli}, a list with one entry per
#' stimulus set (and base image, where the bases have different noise), leaving the default
#' reference untouched.
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
#' # iter is kept tiny here for a fast example; in practice use iter >= 10000.
#' suppressWarnings(generateReferenceDistribution2IFC(rdata_file, iter = 3, ncores = 1))
generateReferenceDistribution2IFC <- function(rdata, iter = 10000, ncores = default_ncores(), response_seed = NULL, save_rdata = TRUE, baseimage = NULL, reference_stimuli = NULL, reference_method = c("gram", "images"), mask = NA) { # nolint: object_length_linter.

  reference_method <- match.arg(reference_method)
  validateIter(iter)
  reference_selection <- selectReferenceBase(rdata, baseimage)
  # Only a proper subset leaves here; an explicit full set takes the default path.
  subset <- subsetReferenceFor(rdata, reference_selection, reference_stimuli, mask)
  if (!is.null(subset)) {
    return(invisible(generateSubsetReference(subset$source, rdata, reference_selection,
                                             subset$reference_stimuli, iter, ncores,
                                             response_seed, save_rdata, reference_method,
                                             subset$masked)))
  }
  if (reference_selection$independent) {
    return(invisible(generateBaseReference(reference_selection, rdata, iter,
                                           ncores, response_seed, save_rdata, reference_method)))
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
    response_seed = response_seed, save_rdata = save_rdata, baseimage = baseimage,
    reference_method = reference_method
  )

  # Removed before load() so the frame re-saved at the end holds the file's own
  # objects of these names, if any, and nothing else; the method is read from
  # .args below.
  rm(reference_stimuli, subset, reference_method, mask)

  # Load parameter file (created when generating stimuli)
  loadRdata(rdata, environment())

  rdata <- .args$rdata
  iter <- .args$iter
  ncores <- .args$ncores
  response_seed <- .args$response_seed
  save_rdata <- .args$save_rdata
  baseimage <- .args$baseimage

  # Shared parameters need only the first base. Inline arguments avoid locals
  # leaking into the frame that is re-saved below.
  if (is.null(response_seed)) requireStimulusSeed(get0("seed", envir = environment(), inherits = FALSE), rdata)
  write("Building the reference from the saved noise, please wait...", stdout())
  write("Computing reference distribution, please wait...", stdout())
  if (iter < 10000) {
    warning("You should set iter >= 10000 for InfoVal statistic to be reliable")
  }

  # Seed the *responses* only. A response_seed replaces the stimulus stream
  # entirely; without one, the draws the generator spent on the parameters are
  # replayed so the responses continue from where they always did. Handing the
  # stimulus seed a different value instead would describe stimuli the
  # participants never saw.
  reference_norms <- referenceNorms(environment(), names(stimuli_params)[1], NULL, iter, ncores,
                                    response_seed, .args$reference_method)

  if (save_rdata) {

    # Save reference norms to rdata file
    write("\nSaving simulated reference distribution to rdata file...", stdout())

    # Provenance belongs to the saved norms; function arguments and scratch state do not.
    reference_norms_seed <- response_seed # nolint: object_usage_linter.
    reference_norms_method <- .args$reference_method # nolint: object_usage_linter.
    reference_norms_source <- "saved_noise" # nolint: object_usage_linter.
    reference_norms_fingerprint <- referenceSnapshot(reference_norms) # nolint: object_usage_linter.
    outfile <- rdata
    internals <- c("stimuli", "responses", "pb", "ci", "i", ".args",
      "rdata", "iter", "ncores", "response_seed", "save_rdata",
      "outfile", "internals", "reference_selection", "baseimage"
    )
    saveRdataSafely(setdiff(ls(all.names = TRUE), internals), outfile, environment())

  }

  invisible(reference_norms)

}
