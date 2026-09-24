#' Generates reference distribution
#'
#' Generates the reference distribution of norms for a stimulus set.
#'
#' \code{\link{computeInfoVal2IFC}} scores a classification image against this distribution. By
#' default the result is saved in the \code{rdata} file for later reuse.
#'
#' @section Reproducibility:
#' With the default \code{response_seed = NULL}, the reference distribution depends only on the
#' stimulus \code{.Rdata} file and the session's \code{\link{RNGkind}}: not on the current
#' random state, and not on \code{ncores}. Two researchers computing InfoVal from the same
#' stimulus file get the same reference distribution, and the same number, on any machine and in
#' any session, provided both use the same RNG kind.
#'
#' The RNG kind is the one gap. \code{set.seed()} keeps whatever kind the session already has,
#' and no stimulus file records which kind was in use. A session with a different
#' \code{RNGkind()} therefore draws different simulated responses, and so a different null; the
#' saved noise basis and stimulus parameters stay the same. To reproduce an earlier reference, use
#' the RNG kind that built it; see \url{https://github.com/rdotsch/rcicr/issues/315}.
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
#' @export
#' @importFrom stats runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus.
#' @param iter Number of simulated classification images, each built from random responses; the distribution holds one norm per image.
#' @param ncores Number of CPU cores used to rebuild the saved noise (default: \code{detectCores() - 1}; 2 under \code{R CMD check}, per CRAN policy).
#' @param response_seed Optional seed for the simulated random responses. The default,
#' \code{NULL}, continues from the state the stimulus generator left behind, as described under
#' Reproducibility. A number gives an independent draw of the null from the same stimuli.
#' @param save_rdata Boolean: write the reference distribution into the \code{rdata} file
#' (default \code{TRUE}). With \code{FALSE}, later calls to \code{\link{computeInfoVal2IFC}}
#' keep using what the file already holds. Set it to \code{FALSE} whenever you set
#' \code{response_seed}, so a one-off null does not become the file's permanent reference.
#' @param baseimage Base-image label, the same key passed to \code{generateCI()}. Required when
#' the base images have different noise parameters. With a single base image, or base images
#' sharing one parameter set, leave it at \code{NULL}.
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
#' \code{reference_norms_by_base} instead, as described above.
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
    reference_norms_fingerprint <- referenceSnapshot(reference_norms) # nolint: object_usage_linter.
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
