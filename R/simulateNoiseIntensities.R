#' Simulate pixel intensity range for noise
#'
#' @export
#' @import matlab
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats runif
#' @importFrom graphics boxplot
#' @param nrep Number of replications
#' @param img_size Size of noise pattern in pixels (one value equal for width and height)
#' @return Matrix with range of noise intensities for each replication
#' @note This function currently always errors: it references undefined \code{data}/\code{by}
#' variables when sizing its progress bar (apparently copy-pasted from \code{\link{batchGenerateCI}}),
#' and it also ignores the \code{img_size} argument internally, always simulating at 512px. Not fixed
#' here; the example below is not run.
#' @examples
#' \dontrun{
#' simulateNoiseIntensities(nrep = 10, img_size = 512)
#' }
simulateNoiseIntensities <- function(nrep=1000, img_size=512) {

  results <- matlab::zeros(nrep, 2)
  # img_size was previously ignored here, hardcoding a 512px noise pattern
  # regardless of what the caller asked for.
  s <- generateNoisePattern(img_size=img_size)

  # The progress bar used to be sized with `data[, by]`, neither of which is a
  # parameter of this function - `data` resolved to utils::data, so every call
  # failed with "object of type 'closure' is not subsettable".
  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  for (i in 1:nrep) {
    setTxtProgressBar(pb, i)

    # One contrast weight per patch index, matching the noise pattern actually
    # generated above rather than a hardcoded 4096.
    params <- (runif(max(s$patchIdx)) * 2) - 1

    noise <- generateNoiseImage(params, s)
    results[i,] <- range(noise)
  }
  close(pb)
  boxplot(results)
  return(results)
}
