#' Generates 2IFC classification image
#'
#' Generate classification image for 2 images forced choice reverse correlation task.  This function exists for backwards compatibility. You can also just use \code{generateCI()}, which this function wraps.
#'
#' This function saves the classification image as PNG to a folder and returns the CI. Your choice of scaling
#' matters. The default is \code{'matched'}, and will match the range of the intensity of the pixels to
#' the range of the base image pixels. This scaling is non linear and depends on the range of both base image
#' and noise pattern. It is truly suboptimal, because it shifts the 0 point of the noise (that is, pixels that would
#' have not changed base image at all before scaling may change the base image after scaling and vice versa). It is
#' however the quick and dirty way to see how the CI noise affects the base image.
#'
#' For more control, use \code{'constant'} scaling, where the scaling is independent of
#' the base image and noise range, but where the choice of constant is arbitrary (provided by the user with t
#' the \code{constant} parameter). The noise is then scale as follows: \code{scaled <- (ci + constant) / (2*constant)}.
#' Note that pixels can take intensity values between 0 and 1 If your scaled noise exceeds those values,
#' a warning will be given. You should pick a higher constant (but do so consistently for different classification images
#' that you want to compare). The higher the constant, the less visible the noise will be in the resulting image.
#'
#' When creating multiple classification images a good strategy is to find the lowest constant that works for all
#' classification images. This can be automatized using the \code{autoscale} function.
#'
#' @export
#' @param stimuli Vector with stimulus numbers (should be numeric) that were presented in the order of the response vector. Stimulus numbers must match those in file name of the generated stimuli.
#' @param responses Vector specifying the responses in the same order of the stimuli vector, coded 1 for original stimulus selected and -1 for inverted stimulus selected.
#' @param baseimage String specifying which base image was used. Not the file name, but the key used in the list of base images at time of generating the stimuli.
#' @param rdata String pointing to .RData file that was created when stimuli were generated. This file contains the contrast parameters of all generated stimuli.
#' @param save_as_png Boolean stating whether to additionally save the CI as PNG image.
#' @param filename Optional string to specify a file name for the PNG image.
#' @param targetpath Optional string specifying path to save PNGs to (default: ./cis).
#' @param antiCI Optional boolean specifying whether antiCI instead of CI should be computed.
#' @param scaling Optional string specifying scaling method: \code{none}, \code{constant}, \code{matched}, or \code{independent} (default).
#' @param constant Optional number specifying the value used as constant scaling factor for the noise (only works for \code{scaling='constant'}).
#' @return List of pixel matrix of classification noise only, scaled classification noise only, base image only and combined.
generateCI2IFC <- function(stimuli, responses, baseimage, rdata, save_as_png=TRUE, filename='', targetpath="./cis", antiCI=FALSE, scaling='independent', constant=0.1) {

  # For backwards compatibility
  return(generateCI(
    stimuli = stimuli,
    responses = responses,
    baseimage = baseimage,
    rdata = rdata,
    save_as_png = save_as_png,
    filename = filename,
    targetpath = targetpath,
    antiCI = antiCI,
    scaling = scaling,
    scaling_constant = constant))

}
