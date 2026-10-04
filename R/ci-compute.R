# Computing classification images from the selected noise parameters. The
# per-participant design lives here; the single-shot design is one call to
# generateCINoise() and stays inline in generateCI().

# One CI per participant, plus their average as the group CI.
#
# Every value the %dopar% body reads is a formal of this function or is built
# above the loop. foreach::getexports() get()s every free variable of the body,
# so one that is unbound aborts the call, even in a branch that cannot run (#235).
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
  # A single selected trial arrives as a parameter vector, which the per-row
  # selection below would index as a matrix.
  if (is.null(dim(params))) params <- matrix(params, nrow = 1)

  pb <- txtProgressBar(min = 0, max = npids, style = 3)

  cl <- startBackend(n_cores, npids)
  if (!is.null(cl)) {
    on.exit(stopClusterSafely(cl), add = TRUE)
  }

  # In parallel each task carries its own participant's rows, so no worker is
  # sent the whole parameter matrix. Serially the rows are taken in the loop,
  # so no second copy of the matrix is made.
  if (is.null(cl)) {
    pid_params <- pid_responses <- vector('list', npids)
  } else {
    pid_params <- lapply(seq_len(npids), function(obs) params[pids == obs, ])
    pid_responses <- lapply(seq_len(npids), function(obs) responses[pids == obs])
  }

  pid.cis <- foreach::foreach(obs = seq_len(npids), obs_params = pid_params, # nolint: object_name_linter.
    obs_responses = pid_responses,
    .combine = 'c',
    .packages = 'rcicr',
    .noexport = c('params', 'responses', 'pid_params', 'pid_responses'),
    .options.snow = progressOption(pb, cl)
  ) %dopar% {

    # Serial path only; in parallel .options.snow ticks the bar in the parent.
    if (is.null(cl)) setTxtProgressBar(pb, obs)

    if (is.null(cl)) {
      pid.rows <- pids == obs # nolint: object_name_linter.
      obs_params <- params[pid.rows, ]
      obs_responses <- responses[pid.rows]
    }

    ci <- generateCINoise(obs_params, obs_responses, p)

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
      # sort(), not unique(): obs indexes factor()'s sorted levels, and
      # sort(unique(x)) keeps the caller's type, so numeric IDs format as
      # before (#267).
      saveToImage(baseimage, combined, paste0(targetpath, '/individual_cis'),
        sort(unique(participants))[obs], antiCI
      )
    }

    return(ci)
  }
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }
  # Blanked so the on.exit() above does not tear down an already-stopped
  # cluster; it resolves cl at exit time, not at registration.
  cl <- NULL
  dim(pid.cis) <- c(img_size, img_size, npids) # nolint: object_name_linter.

  return(list(ci = apply(pid.cis, c(1, 2), mean), pid_cis = pid.cis)) #* sqrt(npids)
}

# The scaling each participant's CI got, or would get, as its own PNG. Only
# 'independent' derives a constant from the CI; the other methods record the
# caller's constant or NA, so they never scan the stack. The stack is unmasked
# (the loop masks only the copy it renders), so the mask is applied here first,
# read once rather than per participant, or a masked-out extreme would set the
# constant. Participants are in the stack's order, factor()'s sorted levels,
# which is the order the PNGs are named in.
individualScalingRecord <- function(pid_cis, participants, mask, img_size, individual_scaling,
                                    individual_scaling_constant) {
  if (scalingMethod(individual_scaling) != 'independent') {
    return(scalingRecord(NULL, individual_scaling, individual_scaling_constant))
  }
  kept <- if (hasMask(mask)) !is.na(applyMask(matrix(0, img_size, img_size), mask, img_size)) else TRUE
  constants <- vapply(seq_len(dim(pid_cis)[3]), function(i) {
    independentConstant(pid_cis[, , i][kept])
  }, numeric(1))
  list(method = 'independent', constant = stats::setNames(constants, sort(unique(participants))))
}
