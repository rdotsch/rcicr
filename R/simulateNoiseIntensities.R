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
#' @examples
#' # nrep and img_size are kept small here so the example is fast; the defaults
#' # (1000 replications at 512px) are what you want for a real estimate.
#' simulateNoiseIntensities(nrep = 10, img_size = 64)
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
