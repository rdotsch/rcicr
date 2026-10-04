#' Check a stimulus file against stimulus PNGs
#'
#' Compares the noise a stimulus \code{.Rdata} file records with the stimulus PNGs in a folder,
#' trial by trial. An original and an inverted stimulus hold the same base image, so their
#' difference is the trial's noise alone: \code{ori - inv} equals the noise divided by 0.6
#' wherever neither image is clipped at black or white. The base images are therefore not needed.
#'
#' Use it to confirm that a regenerated \code{.Rdata} file matches an archive of stimulus PNGs whose
#' original file was lost; \code{vignette("recipes", package = "rcicr")}, "When the
#' \code{.Rdata} file is lost", shows how.
#'
#' @export
#' @param rdata Path to the \code{.Rdata} file to check, as written by
#'   \code{\link{generateStimuli2IFC}}.
#' @param png_dir Directory holding the stimulus PNGs. Required; there is no default.
#' @param label,seed The \code{label} and \code{seed} the PNGs were generated with, which their
#'   names carry. Default: the values stored in \code{rdata}. Pass the archive's own when
#'   \code{rdata} is a candidate regenerated with other settings, so the PNGs can be found.
#' @return A data frame with one row per base label and trial: \code{base}, \code{trial},
#'   \code{share} (the share of compared pixels whose \code{ori - inv} agrees with the file's noise
#'   divided by 0.6 to within 1/255), \code{compared} (the number of pixels neither image clips) and
#'   \code{missing} (\code{TRUE} when either PNG is absent, with \code{share} \code{NA}). Its
#'   \code{unchecked} attribute lists PNGs in \code{png_dir} that carry \code{label} and
#'   \code{seed} but belong to no row, such as trials beyond the file's \code{n_trials}.
#'
#'   No verdict is returned. In the configurations measured in
#'   \url{https://github.com/rdotsch/rcicr/blob/main/analyses/stimulus-png-residuals.md}, a wrong
#'   \code{nscales}, seed, noise type, base order or \code{use_same_parameters} agreed on 2.1\% to
#'   7.4\% of pixels per trial, but a Gabor \code{sigma} of 24 where the PNGs used 25 on 95\% to
#'   99.6\%: compare candidates against the same PNGs. Check every base label, since with several bases a
#'   wrong \code{use_same_parameters} shows only after the first.
#'
#'   Warns about missing and unchecked PNGs, and stops if no PNG named for \code{label} and
#'   \code{seed} exists.
#' @seealso \code{vignette("recipes", package = "rcicr")}, "When the \code{.Rdata} file is lost".
#' @examples
#' base_face <- tempfile(fileext = ".png")
#' png::writePNG(matrix(runif(32 * 32), 32, 32), base_face)
#' stimulus_path <- tempfile("stimuli")
#' generateStimuli2IFC(list(face = base_face), n_trials = 3, img_size = 32,
#'                     stimulus_path = stimulus_path, seed = 1, ncores = 1, nscales = 2)
#' rdata <- list.files(stimulus_path, pattern = "\\.Rdata$", full.names = TRUE)
#' checkStimulusPNGs(rdata, stimulus_path)
checkStimulusPNGs <- function(rdata, png_dir, label = NULL, seed = NULL) {
  if (missing(png_dir) || !is.character(png_dir) || length(png_dir) != 1L || !dir.exists(png_dir)) {
    stop('png_dir must be the path of an existing directory of stimulus PNGs.', call. = FALSE)
  }
  loaded <- loadStimulusParams(rdata, require_img_size = FALSE)
  if (is.null(label) || is.null(seed)) {
    stored <- new.env(parent = emptyenv())
    loadRdata(rdata, stored)
    if (is.null(label)) label <- get0('label', envir = stored, inherits = FALSE)
    if (is.null(seed)) seed <- get0('seed', envir = stored, inherits = FALSE)
    if (is.null(label) || is.null(seed)) {
      stop('rdata stores no ', if (is.null(label)) 'label' else 'seed',
           '. Pass the label and seed the PNGs were generated with.', call. = FALSE)
    }
  }

  rows <- list()
  checked <- character(0)
  for (base in names(loaded$stimuli_params)) {
    n_trials <- NROW(loaded$stimuli_params[[base]])
    params <- selectStimulusParams(loaded$stimuli_params, base, seq_len(n_trials))
    if (is.null(dim(params))) params <- matrix(params, nrow = 1L)
    for (trial in seq_len(n_trials)) {
      ori_file <- stimulusPngPath(png_dir, label, base, seed, trial, 'ori')
      inv_file <- stimulusPngPath(png_dir, label, base, seed, trial, 'inv')
      checked <- c(checked, basename(ori_file), basename(inv_file))
      rows[[length(rows) + 1L]] <- comparePngPair(ori_file, inv_file, params[trial, ], loaded$p, base,
                                                  trial)
    }
  }
  result <- do.call(rbind, rows)

  archived <- stimulusPngNames(png_dir, label, seed)
  if (length(archived) == 0L) {
    stop('No stimulus PNG named for label "', label, '" and seed ', seed, ' in ', png_dir,
         ' (looked for names such as ', basename(stimulusPngPath(png_dir, label, result$base[1], seed,
                                                                 1L, 'ori')),
         '). label and seed must be the ones the PNGs were generated with.', call. = FALSE)
  }
  if (any(result$missing)) {
    absent <- result[result$missing, ]
    warning(sum(result$missing), ' trial(s) have no ori or inv PNG: ',
            firstFew(paste0(absent$base, ' ', absent$trial)), '.', call. = FALSE)
  }

  unchecked <- setdiff(archived, checked)
  if (length(unchecked) > 0L) {
    warning(length(unchecked), ' PNG(s) named for label "', label, '" and seed ', seed,
            ' are outside the trials and base labels rdata holds: ', firstFew(unchecked), '.',
            call. = FALSE)
  }
  attr(result, 'unchecked') <- unchecked
  return(result)
}

comparePngPair <- function(ori_file, inv_file, params, p, base, trial) {
  if (!file.exists(ori_file) || !file.exists(inv_file)) {
    return(data.frame(base = base, trial = trial, share = NA_real_, compared = 0L, missing = TRUE))
  }
  ori <- firstChannel(png::readPNG(ori_file))
  inv <- firstChannel(png::readPNG(inv_file))
  noise <- generateNoiseImage(params, p)
  compared <- ori > 0 & ori < 1 & inv > 0 & inv < 1
  agrees <- abs(ori - inv - noise / stimulusNoiseScale)[compared] <= 1 / 255
  data.frame(base = base, trial = trial, share = if (any(compared)) mean(agrees) else NA_real_,
             compared = sum(compared), missing = FALSE)
}

# Stimulus PNGs are written grey; an RGB(A) copy holds the grey in each colour channel.
firstChannel <- function(img) {
  if (length(dim(img)) == 3L) img[, , 1] else img
}

# Every PNG in png_dir named as stimulusPngPath() names one for this label and seed. The base
# label may contain underscores, so the name is matched by its prefix and suffix; %05d is a minimum
# width, so the trial has five or more digits.
stimulusPngNames <- function(png_dir, label, seed) {
  files <- list.files(png_dir, all.files = TRUE)
  suffix <- paste0('_', escapeRegex(paste(seed)), '_[0-9]{5,}_(ori|inv)\\.png$')
  files[startsWith(files, paste0(label, '_')) & grepl(suffix, files)]
}

escapeRegex <- function(x) {
  gsub('([][{}()+*^$|\\\\.?])', '\\\\\\1', x)
}
