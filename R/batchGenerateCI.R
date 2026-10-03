#' Generates multiple classification images by participant or condition
#'
#' Generate one classification image per participant or condition, for any reverse correlation task.
#'
#' Splits \code{data} by the \code{by} column, calls \code{\link{generateCI}} for each part, and returns the CIs. By default each CI is also saved as a PNG.
#'
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @param data Data frame with one row per trial.
#' @param by Name of the column that splits the data into units, such as participants or conditions. One CI is computed per unit.
#' @param stimuli Name of the column holding the stimulus numbers of the presented stimuli.
#' @param responses Name of the column holding the responses: 1 where the original stimulus was chosen, -1 where the inverted one was. The column must be numeric with no missing or infinite values; the call stops before computing anything otherwise.
#' @param baseimage String naming the base image: not its file name, but its key in the \code{base_face_files} list passed to \code{\link{generateStimuli2IFC}}.
#' @param rdata Path to the \code{.Rdata} file written when the stimuli were generated. It holds the contrast parameters of every stimulus.
#' @param save_as_png Boolean: also save the CI as a PNG image.
#' @param targetpath Directory to save PNGs to. Required when \code{save_as_png = TRUE}; there is no default. The directory is created if it does not exist; to just try the function out, use \code{tempdir()}.
#' @param label Optional string added to the PNG file names to make them easier to identify.
#' @param antiCI Boolean: compute the anti-CI, the classification image with its sign flipped, instead of the CI.
#' @param scaling Scaling method: \code{none}, \code{constant}, \code{matched}, \code{independent} or \code{autoscale} (default). \code{autoscale} computes the CIs unscaled and then puts them on one scale with \code{\link{autoscale}}.
#' @param constant Scaling constant for the noise. Used only when \code{scaling = 'constant'}.
#' @param participants Optional name of a column holding a participant ID per trial. When given, each unit's CI is computed as \code{\link{generateCI}} does with \code{participants}: one CI per participant in the unit, then their average. With \code{by} naming a condition, that gives one CI per condition with participants nested in it, and each participant counts equally however many trials they contributed. Every row needs an ID; the call stops before computing anything if any is missing. Default: \code{NULL}, all of a unit's trials pooled into one CI.
#' @return Named list with one classification image per unit. Each is itself a list of pixel matrices: the raw noise (\code{ci}), the scaled noise (\code{scaled}), the base image (\code{base}) and the two combined (\code{combined}).
#' @seealso \code{vignette("reverse-correlation-walkthrough")}, section "Several participants at once".
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
#'
#' cis <- suppressWarnings(batchGenerateCI(
#'   data = data, by = "participant", stimuli = "stimulus", responses = "response",
#'   baseimage = "face", rdata = rdata_file, save_as_png = FALSE
#' ))
batchGenerateCI <- function(data, by, stimuli, responses, baseimage, rdata, save_as_png = TRUE, targetpath, label = '', antiCI = FALSE, scaling = 'autoscale', constant = 0.1, participants = NULL) {
  return(batchCIs(
    data = data, by = by, stimuli = stimuli, responses = responses, baseimage = baseimage,
    rdata = rdata, save_as_png = save_as_png, targetpath = targetpath, label = label,
    antiCI = antiCI, scaling = scaling, constant = constant, participants = participants
  ))
}
