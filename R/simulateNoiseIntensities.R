#' Simulate pixel intensity range for noise
#'
#' @export
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats runif
#' @importFrom graphics boxplot
#' @param nrep Number of replications
#' @param img_size Size of noise pattern in pixels (one value equal for width and height). The pattern has five scales, so it must be a multiple of 16.
#' @return Matrix with range of noise intensities for each replication
#' @examples
#' # nrep and img_size are kept small here so the example is fast; the defaults
#' # (1000 replications at 512px) are what you want for a real estimate.
#' simulateNoiseIntensities(nrep = 10, img_size = 64)
simulateNoiseIntensities <- function(nrep = 1000, img_size = 512) {

  results <- array(0, c(nrep, 2))
  s <- renderingBasis(generateNoisePattern(img_size = img_size))

  pb <- txtProgressBar(min = 0, max = nrep, style = 3)
  for (i in 1:nrep) {
    setTxtProgressBar(pb, i)

    # One contrast weight per patch index of the pattern above.
    params <- (runif(max(s$patchIdx)) * 2) - 1

    noise <- generateNoiseImage(params, s)
    results[i, ] <- range(noise)
  }
  close(pb)
  boxplot(results)
  return(results)
}
