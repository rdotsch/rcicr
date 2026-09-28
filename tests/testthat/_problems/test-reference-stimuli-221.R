# Extracted from test-reference-stimuli.R:221

# prequel ----------------------------------------------------------------------
quiet <- function(expr) {
  out <- NULL
  suppressWarnings(utils::capture.output(out <- expr))
  out
}
reference_of <- function(rdata, iter = 20, ...) {
  quiet(generateReferenceDistribution2IFC(rdata, iter = iter, ncores = 1, ...))
}
messages_of <- function(expr) {
  seen <- character()
  quiet(withCallingHandlers(expr, message = function(m) {
    seen <<- c(seen, conditionMessage(m))
    invokeRestart("muffleMessage")
  }))
  seen
}
ci_of <- function(rdata, stimuli, participants = NA) {
  responses <- rep(c(1, -1), length.out = length(stimuli))
  quiet(generateCI(stimuli, responses, "base", rdata, participants = participants,
                   save_as_png = FALSE, n_cores = 1))
}
saved <- function(rdata) {
  e <- new.env()
  load(rdata, envir = e)
  as.list(e, all.names = TRUE)
}
direct_reference <- function(rdata, ids, iter, response_seed = NULL) {
  f <- saved(rdata)
  params <- f$stimuli_params$base
  noise <- vapply(ids, function(i) as.vector(generateNoiseImage(params[i, ], f$p)),
                  numeric(length(f$p$patches[, , 1])))
  if (is.null(response_seed)) {
    set.seed(f$seed)
    for (trial in seq_len(f$n_trials)) runif(ncol(params))
  } else {
    set.seed(response_seed)
  }
  vapply(seq_len(iter), function(i) {
    r <- ((runif(length(ids)) > 0.5) * 2) - 1
    norm((noise %*% r) / length(ids), "f")
  }, numeric(1))
}

# test -------------------------------------------------------------------------
tmp <- withr::local_tempdir()
rdata <- make_fixture_rdata(tmp, n_trials = 6)
score <- function(ci, ...) messages_of(computeInfoVal2IFC(ci, rdata, iter = 5, response_seed = 1, ...))
