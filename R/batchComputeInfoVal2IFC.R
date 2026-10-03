#' Computes Informational Values for several classification images
#'
#' Computes the Informational Value of every classification image in a list, from one 2IFC
#' stimulus file, as \code{\link{computeInfoVal2IFC}} would for each in turn.
#'
#' Each distinct reference distribution is found or simulated once and used for every
#' classification image scored against it. The stimulus file is read at most three times to find
#' stored references, however many images there are, and again for each reference that has to be
#' simulated. A loop over \code{computeInfoVal2IFC()} reads the file twice
#' per image and, when the reference cannot be stored (a read-only file, a \code{response_seed},
#' or \code{force_gen_ref_dist = TRUE}), simulates it again for every image. The values are
#' identical to that loop's.
#'
#' A masked classification image is scored over its unmasked pixels, as described under "Masked
#' classification images" in \code{\link{computeInfoVal2IFC}}. Images with the same stimuli and
#' the same mask share one reference.
#'
#' Messages about the trial design (see "Matching the reference to the classification image" in
#' \code{\link{computeInfoVal2IFC}}) are collected into one message per kind, naming the
#' classification images concerned, rather than one per image.
#'
#' @export
#' @param target_cis A list of classification images, each as returned by
#' \code{\link{generateCI}}; for example the result of \code{\link{batchGenerateCI2IFC}}.
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. Every
#' classification image in \code{target_cis} must come from it.
#' @param iter,force_gen_ref_dist,response_seed,reference_method As in
#' \code{\link{computeInfoVal2IFC}}, applied to every classification image.
#' @param baseimage As in \code{\link{computeInfoVal2IFC}}: the one base-image label every
#' classification image in \code{target_cis} was computed for. Score CIs for different base
#' images in separate calls.
#' @param reference_stimuli The stimulus numbers each classification image was built from, when
#' that is not every saved stimulus. \code{NULL} (the default) scores every image against the
#' full set; a vector scores every image against that subset; a list with one element per image
#' (\code{NULL} for the full set) gives each its own. To use what \code{\link{generateCI}}
#' recorded, pass \code{lapply(target_cis, function(ci) attr(ci, "trial_design")$stimuli)}.
#' @return A numeric vector of Informational Values, one per classification image, named as
#' \code{target_cis} is.
#' @seealso \code{vignette("recipes", package = "rcicr")}, "Several classification images at once".
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
#' # two "participants", three trials each
#' data <- data.frame(
#'   participant = rep(c("p1", "p2"), each = 3),
#'   stimulus = 1:6,
#'   response = sample(c(1, -1), 6, replace = TRUE)
#' )
#' cis <- suppressWarnings(batchGenerateCI2IFC(
#'   data = data, by = "participant", stimuli = "stimulus", responses = "response",
#'   baseimage = "face", rdata = rdata_file, save_as_png = FALSE, scaling = "none"
#' ))
#'
#' # iter is kept tiny here for a fast example, in practice use iter >= 10000.
#' # Each participant saw three of the six stimuli, so each is scored over those.
#' suppressWarnings(batchComputeInfoVal2IFC(
#'   cis, rdata_file, iter = 20,
#'   reference_stimuli = list(1:3, 4:6)
#' ))
batchComputeInfoVal2IFC <- function(target_cis, rdata, iter = 10000, force_gen_ref_dist = FALSE, response_seed = NULL, baseimage = NULL, reference_stimuli = NULL, reference_method = c("gram", "images")) { # nolint: object_length_linter.

  reference_method <- match.arg(reference_method)
  validateCiBatch(target_cis)
  n <- length(target_cis)
  labels <- ciLabels(target_cis)
  requested <- if (is.list(reference_stimuli)) {
    if (length(reference_stimuli) != n) {
      stop('reference_stimuli has ', length(reference_stimuli), ' elements for ', n,
           ' classification images; give one per image, or a single vector for all of them.',
           call. = FALSE)
    }
    reference_stimuli
  } else {
    rep(list(reference_stimuli), n)
  }

  selection <- selectReferenceBase(rdata, baseimage)
  source <- referenceSource(rdata, selection)
  if (!validTrialCount(source$n_trials)) stop('The stimulus file must contain a positive integer n_trials.')
  canonical <- lapply(requested, canonicalReferenceStimuli, n_trials = source$n_trials)
  masks <- lapply(target_cis, function(ci) referenceMask(ciMask(ci), source$img_size))
  keys <- vapply(seq_len(n), function(i) {
    key <- maskKey(masks[[i]])
    paste(paste(canonical[[i]], collapse = ','), paste(key$lengths, collapse = ','),
          key$values[1], sep = '|')
  }, character(1))
  groups <- split(seq_len(n), factor(keys, levels = unique(keys)))

  issues <- lapply(seq_len(n), function(i) {
    trialDesignIssue(target_cis[[i]], source$n_trials, canonical[[i]])
  })
  reportBatchTrialDesign(issues, labels)

  info_vals <- numeric(n)
  quiet <- function(n_trials, reference_stimuli) NULL
  for (members in groups) {
    ids <- canonical[[members[1]]]
    masked <- masks[[members[1]]]
    reference <- if (!is.null(ids) || !is.null(masked)) {
      subsetReference(rdata, iter, force_gen_ref_dist, response_seed, source, selection, ids,
                      reference_method, quiet, masked)
    } else if (selection$independent) {
      baseReference(rdata, iter, force_gen_ref_dist, response_seed, selection, reference_method,
                    quiet)
    } else {
      sharedReference(rdata, iter, force_gen_ref_dist, response_seed, reference_method, quiet)
    }
    write(paste0('Reference for ', length(members), ' of ', n, ' classification images (',
                 reference$note, 'reference median = ', reference$median, '; MAD = ',
                 reference$mad, '; iterations = ', reference$iter, ')'), stdout())
    for (i in members) {
      cinorm <- ciNorm(target_cis[[i]])
      info_vals[i] <- (cinorm - reference$median) / reference$mad
    }
  }

  names(info_vals) <- names(target_cis)
  return(info_vals)
}

isClassificationImage <- function(x) {
  is.list(x) && is.matrix(x[['ci']]) && is.numeric(x[['ci']])
}

# A single CI is a list too; its own ci element, not the name, gives it away,
# since list(ci = first, control = second) is a valid batch.
validateCiBatch <- function(target_cis) {
  if (isClassificationImage(target_cis)) {
    stop('target_cis is a single classification image. Score it with computeInfoVal2IFC(), ',
         'or wrap it in a list.', call. = FALSE)
  }
  if (!is.list(target_cis) || !length(target_cis)) {
    stop('target_cis must be a list of classification images, as returned by generateCI().',
         call. = FALSE)
  }
  bad <- which(!vapply(target_cis, isClassificationImage, logical(1)))
  if (length(bad)) {
    stop('Element ', paste(ciLabels(target_cis)[bad], collapse = ', '), ' of target_cis is ',
         'not a classification image: expected a list with a numeric matrix ci, as returned ',
         'by generateCI().', call. = FALSE)
  }
  invisible(NULL)
}

ciLabels <- function(target_cis) {
  labels <- names(target_cis)
  positions <- as.character(seq_along(target_cis))
  if (is.null(labels)) return(positions)
  ifelse(is.na(labels) | !nzchar(labels), positions, labels)
}

# One message per kind of issue, not one per CI: with a hundred participants,
# a hundred copies would bury the one that differs.
reportBatchTrialDesign <- function(issues, labels) {
  n <- length(issues)
  some <- function(which) {
    shown <- firstFew(labels[which])
    paste0(length(which), ' of the ', n, ' classification images (', shown, ')')
  }
  averaged <- which(vapply(issues, function(x) isTRUE(x$averaged), logical(1)))
  mismatched <- which(vapply(issues, function(x) identical(x$averaged, FALSE), logical(1)))
  if (length(averaged)) {
    message(some(averaged), ' average several participants or repeated presentations. ',
            uncalibratedDesignAdvice())
  }
  if (length(mismatched)) {
    message(some(mismatched), ' were built from different stimuli than the reference they are ',
            'scored against. ', mismatchedDesignAdvice(
              'reference_stimuli = lapply(<your CIs>, function(ci) attr(ci, "trial_design")$stimuli)'
            ))
  }
  invisible(NULL)
}
