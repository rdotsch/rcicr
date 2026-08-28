# The gate on the latent-space module.
#
# Every exported function in that module opens with this call. The point is not
# to hide the code -- it is exported, documented and checked -- but to make a
# user state that they accept a weaker guarantee than the rest of the package
# offers: while the module is experimental its numeric output is not covered by
# the promise that a stored analysis script reproduces its results. A latent
# classification image computed today may not reproduce under the next
# development build.
#
# The gate comes off once a golden master pins the numbers and a release has
# shipped them, at which point the module joins the ordinary guarantee.

requireExperimental <- function(caller) {
  if (isTRUE(getOption('rcicr.experimental', FALSE))) {
    return(invisible(TRUE))
  }

  msg <- paste0(
    caller, '() is part of the experimental latent-space module, ',
    'which is switched off by default. Its numeric output is not yet stable: ',
    'unlike the rest of rcicr, a classification image it computes today is ',
    'not guaranteed to reproduce under a later version. To use it anyway, run ',
    'options(rcicr.experimental = TRUE) first, and record the rcicr version ',
    'alongside any result you keep.'
  )

  stop(msg, call. = FALSE)
}
