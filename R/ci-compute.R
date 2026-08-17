# Computing classification images from the selected noise parameters. The
# per-participant design lives here; the single-shot design is one call to
# generateCINoise() and stays inline in generateCI().

# One CI per participant, plus their average as the group CI.
#
# Every value the %dopar% body reads is a formal of this function or is built
# above the loop, and that is load-bearing rather than tidiness.
# foreach::getexports() scans the body for free variables and get()s each one
# from the environment the loop is evaluated in, so what is in scope decides
# whether the call runs at all -- #235 was an absent-and-required `targetpath`
# aborting every participant call at the default core count, from a branch that
# could not even execute.
#
# Returns both the group CI and the per-participant stack: the t.test z-map
# short-circuits to the latter instead of rebuilding a noise image per trial,
# so dropping it would sever that path silently.
computeParticipantCIs <- function(params, responses, participants, p, base,
                                  baseimage, img_size, mask, n_cores,
                                  save_individual_cis, targetpath,
                                  individual_scaling,
                                  individual_scaling_constant, antiCI) {
  pids <- as.numeric(factor(participants))
  npids <- length(unique(pids))

  pb <- txtProgressBar(min = 1, max = npids, style = 3)

  cl <- startBackend(n_cores)
  if (!is.null(cl)) {
    on.exit(stopClusterSafely(cl), add = TRUE)
  }

  pid.cis <- foreach::foreach(obs = 1:npids, # nolint: object_name_linter.
    .combine = 'c',
    .packages = 'rcicr',
    .options.snow = progressOption(pb, cl)
  ) %dopar% {

    # Serial path only; in parallel .options.snow ticks the bar in the parent.
    if (is.null(cl)) setTxtProgressBar(pb, obs)

    # Select only the observations of the current participant
    pid.rows <- pids == obs # nolint: object_name_linter.

    # Construct the noise pattern
    ci <- generateCINoise(params[pid.rows, ], responses[pid.rows], p)

    # Check if individual CIs should be saved. If so, generate and save them
    if (save_individual_cis) {
      if (hasMask(mask)) {
        individual_ci <- applyMask(ci, mask, img_size)
      } else {
        individual_ci <- ci
      }
      scaled <- applyScaling(base, individual_ci, individual_scaling,
        individual_scaling_constant
      )
      combined <- combine(scaled, base)
      # sort(), not unique(): obs indexes the *sorted* factor levels built
      # above, so naming the file from appearance order gave every participant
      # someone else's ID whenever the two orders disagreed. sort(unique(x))
      # is factor()'s own level order, and keeps the caller's type so numeric
      # IDs still format exactly as they did.
      saveToImage(baseimage, combined, paste0(targetpath, '/individual_cis'),
        sort(unique(participants))[obs], antiCI
      )
    }

    # Return the CI
    return(ci)
  }
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }
  # Blanked so the on.exit() above does not tear down an already-stopped
  # cluster; it resolves cl at exit time, not at registration.
  cl <- NULL
  dim(pid.cis) <- c(img_size, img_size, npids) # nolint: object_name_linter.

  # Average across participants for final CI and return to original variance
  return(list(ci = apply(pid.cis, c(1, 2), mean), pid_cis = pid.cis)) #* sqrt(npids)
}
