# The two ways generateCI() turns a classification image into a z-map. Both
# return the z-map matrix and write nothing; plotZmap() does the saving.

# The cheap z-map: blur the CI, rescale to z-scores, drop everything inside the
# threshold.
computeZmapQuick <- function(ci, sigma, threshold, img_size) {
  zmap <- as.matrix(blur(as.im(ci), sigma = sigma))
  zmap <- matrix(scale(as.vector(zmap)), img_size, img_size)

  zmap[zmap > -threshold & zmap < threshold] <- NA

  return(zmap)
}

# The t-test z-map: a one-sample t test per pixel across a stack of images.
#
# What that stack is depends on how the CI was built, and this is the one place
# in the pipeline that reaches across the two designs. With `participants` given
# the per-participant CIs are the stack -- they already exist, so building a
# noise image per trial would be both wasteful and a different test. Without it
# there is nothing to reuse and the stack is built here.
#
# That switch is `is.null(pid_cis)` rather than a `participants` argument
# because the dependency is on the CIs themselves: pid_cis is produced if and
# only if participants was given, so the two are equivalent, but only one of
# them says what this function actually needs.
computeZmapTTest <- function(ci, params, responses, p, pid_cis, img_size,
                             n_cores) {
  if (is.null(pid_cis)) {
    # Weigh the stimulus parameters of each trial using the given responses
    weightedparameters <- params * responses
    n_observations <- length(responses)

    pb <- txtProgressBar(min = 1, max = n_observations, style = 3)

    cl <- startBackend(n_cores)
    if (!is.null(cl)) {
      on.exit(stopClusterSafely(cl), add = TRUE)
    }

    # For each weighted stimulus, construct the complementary noise pattern
    noiseimages <- foreach::foreach(obs = 1:n_observations, .combine = 'c',
      .packages = 'rcicr',
      .options.snow = progressOption(pb, cl)
    ) %dopar% {
      noiseimage <- generateNoiseImage(weightedparameters[obs, ], p)
      if (is.null(cl)) setTxtProgressBar(pb, obs)
      return(noiseimage)
    }
    if (!is.null(cl)) {
      parallel::stopCluster(cl)
    }
    # Blanked so the on.exit() above does not tear down an already-stopped
    # cluster; it resolves cl at exit time, not at registration.
    cl <- NULL
    dim(noiseimages) <- c(img_size, img_size, n_observations)

  } else {
    noiseimages <- pid_cis
  }

  # Get p value for each pixel
  pmap <- apply(noiseimages, 1:2, function(x) unlist(t.test(x)['p.value']))

  return(sign(ci) * abs(qnorm(pmap / 2)))
}
