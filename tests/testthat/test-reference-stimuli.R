# reference_stimuli builds the InfoVal reference over the stimuli a CI was
# built from (#349), and generateCI() records its design so
# computeInfoVal2IFC() can say when the two do not match.

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

# The arithmetic of the reference, spelled out over the selected stimuli.
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

test_that("an explicit full set, in any order and storage mode, is the default reference", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 8)

  default <- reference_of(rdata, response_seed = 3, save_rdata = FALSE)
  explicit <- reference_of(rdata, response_seed = 3, save_rdata = FALSE,
                           reference_stimuli = as.numeric(c(8, 1:7)))
  expect_identical(explicit, default)

  # The default path re-saves its frame, so neither the argument nor anything
  # bound while deciding on the path may reach the file, omitted or not.
  added <- c("reference_norms", "reference_norms_fingerprint", "reference_norms_seed",
             "reference_norms_source")
  fresh <- names(saved(rdata))
  reference_of(rdata)
  expect_setequal(setdiff(names(saved(rdata)), fresh), added)
  other <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  reference_of(other, reference_stimuli = 1:8)
  expect_setequal(setdiff(names(saved(other)), fresh), added)
  expect_identical(saved(other)$reference_norms, saved(rdata)$reference_norms)
})

test_that("a subset reference is the default arithmetic over those stimuli", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 12)
  ids <- c(2L, 5L, 6L, 9L, 11L)

  seeded <- reference_of(rdata, response_seed = 4, save_rdata = FALSE, reference_stimuli = ids)
  expect_equal(seeded, direct_reference(rdata, ids, 20, response_seed = 4), tolerance = 1e-12)

  # Without a response_seed the stream replays the generator's draws for every
  # saved trial, so the file alone reproduces it.
  replayed <- reference_of(rdata, save_rdata = FALSE, reference_stimuli = rev(ids))
  expect_equal(replayed, direct_reference(rdata, ids, 20), tolerance = 1e-12)

  full <- reference_of(rdata, response_seed = 4, save_rdata = FALSE)
  expect_false(isTRUE(all.equal(seeded, full)))
})

test_that("a subset CI is scored against the reference over its own stimuli", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 12)
  ids <- c(1, 4, 7, 8, 10)
  ci <- ci_of(rdata, ids)

  norms <- direct_reference(rdata, ids, 20, response_seed = 6)
  expected <- (norm(matrix(ci$ci), "f") - median(norms)) / mad(norms)
  matched <- quiet(computeInfoVal2IFC(ci, rdata, iter = 20, response_seed = 6,
                                      reference_stimuli = ids))
  expect_equal(matched, expected, tolerance = 1e-10)

  default <- suppressMessages(quiet(computeInfoVal2IFC(ci, rdata, iter = 20, response_seed = 6)))
  expect_false(isTRUE(all.equal(default, matched)))
})

test_that("a subset reference is cached apart from the default one", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 12)
  reference_of(rdata)
  before <- saved(rdata)

  ids <- c(3, 4, 9)
  first <- reference_of(rdata, reference_stimuli = ids)
  after <- saved(rdata)
  expect_identical(setdiff(names(after), names(before)), "reference_norms_by_stimuli")
  expect_identical(after[names(before)], before)

  ci <- ci_of(rdata, ids)
  out <- utils::capture.output(z <- computeInfoVal2IFC(ci, rdata, reference_stimuli = as.integer(ids)))
  expect_true(any(grepl("Using reference distribution found", out, fixed = TRUE)))
  expect_equal(z, (norm(matrix(ci$ci), "f") - median(first)) / mad(first))
  expect_length(saved(rdata)$reference_norms_by_stimuli, 1)

  reference_of(rdata, reference_stimuli = c(3, 4))
  cache <- saved(rdata)$reference_norms_by_stimuli
  expect_length(cache, 2)
  expect_identical(lapply(cache, `[[`, "reference_stimuli"), list(c(3L, 4L, 9L), c(3L, 4L)))
})

test_that("reference_stimuli must be distinct saved stimulus numbers", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 6)

  expect_error(reference_of(rdata, reference_stimuli = c(1, 2, 2)), "more than once")
  expect_error(reference_of(rdata, reference_stimuli = c(1, 7)), "beyond the 6 trials")
  expect_error(reference_of(rdata, reference_stimuli = c(1, 2.5)), "whole-number")
  expect_error(reference_of(rdata, reference_stimuli = integer(0)), "at least one")
  expect_error(reference_of(rdata, reference_stimuli = factor(1:3)), "not factors")
})

test_that("a file's own reference_stimuli object survives both save paths", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 6)
  planted <- c(9, 9, 9)
  mutate_rdata(rdata, reference_stimuli = planted)

  # The default path re-saves its frame; the argument must neither be
  # overridden by the file's object nor replace it.
  explicit <- reference_of(rdata, response_seed = 2, save_rdata = FALSE, reference_stimuli = 1:6)
  expect_identical(explicit, reference_of(rdata, response_seed = 2, save_rdata = FALSE))
  reference_of(rdata, reference_stimuli = 1:6)
  expect_identical(saved(rdata)$reference_stimuli, planted)

  expect_equal(reference_of(rdata, response_seed = 2, reference_stimuli = c(2, 3)),
               direct_reference(rdata, c(2, 3), 20, response_seed = 2), tolerance = 1e-12)
  expect_identical(saved(rdata)$reference_stimuli, planted)
})

test_that("independent bases keep a subset reference per base", {
  tmp <- withr::local_tempdir()
  rdata <- make_independent_fixture(tmp)

  first <- reference_of(rdata, baseimage = "first", reference_stimuli = c(1, 5, 9))
  second <- reference_of(rdata, baseimage = "second", reference_stimuli = c(1, 5, 9))
  expect_false(isTRUE(all.equal(first, second)))

  cache <- saved(rdata)$reference_norms_by_stimuli
  expect_identical(vapply(cache, `[[`, "", "baseimage"), c("first", "second"))
  expect_null(saved(rdata)$reference_norms_by_base)
})

test_that("a seedless file asks for a response_seed on the subset path too", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 6)
  mutate_rdata(rdata, .remove = "seed")

  expect_error(reference_of(rdata, reference_stimuli = c(1, 2)), "reference_stimuli = <the same stimuli>")
  expect_length(reference_of(rdata, response_seed = 1, reference_stimuli = c(1, 2)), 20)
})

test_that("generateCI records its trial design without adding a field", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 6)

  ci <- ci_of(rdata, c(5, 2, 3))
  expect_named(ci, c("ci", "scaled", "base", "combined"))
  expect_identical(attr(ci, "trial_design"),
                   list(stimuli = c(2L, 3L, 5L), repeated = FALSE, n_participants = 1L))

  expect_true(attr(ci_of(rdata, c(1, 2, 2)), "trial_design")$repeated)
  across <- attr(ci_of(rdata, c(1:6, 1:6), participants = rep(c("a", "b"), each = 6)), "trial_design")
  expect_false(across$repeated)
  expect_identical(across$n_participants, 2L)
  expect_true(attr(ci_of(rdata, c(1, 1, 2, 3), participants = c("a", "a", "b", "b")),
                   "trial_design")$repeated)

  wrapped <- quiet(generateCI2IFC(c(5, 2, 3), c(1, -1, 1), "base", rdata, save_as_png = FALSE))
  expect_identical(attr(wrapped, "trial_design"), attr(ci, "trial_design"))
  trials <- data.frame(g = c("x", "x", "x"), s = c(5, 2, 3), r = c(1, -1, 1))
  batch <- quiet(batchGenerateCI(trials, by = "g", stimuli = "s", responses = "r",
                                 baseimage = "base", rdata = rdata, save_as_png = FALSE))
  expect_identical(attr(batch[[1]], "trial_design"), attr(ci, "trial_design"))
})

test_that("the guard names a subset mismatch and recommends reference_stimuli", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 6)
  score <- function(ci, ...) messages_of(computeInfoVal2IFC(ci, rdata, iter = 5, response_seed = 1, ...))

  subset_ci <- ci_of(rdata, c(1, 3, 4))
  unmatched <- score(subset_ci)
  expect_length(unmatched, 1)
  expect_match(unmatched, "built from 3 of the 6 saved stimuli, but the reference is built over 6")
  expect_match(unmatched, "reference_stimuli", fixed = TRUE)
  expect_match(score(subset_ci, reference_stimuli = c(1, 3)), "built over 2")

  expect_length(score(subset_ci, reference_stimuli = c(4, 1, 3)), 0)
  expect_length(score(ci_of(rdata, 1:6)), 0)
  expect_length(score(list(ci = subset_ci$ci)), 0)
})

test_that("the guard never recommends reference_stimuli for repeats or participant averages", {
  tmp <- withr::local_tempdir()
  rdata <- make_fixture_rdata(tmp, n_trials = 6)
  score <- function(ci, ...) messages_of(computeInfoVal2IFC(ci, rdata, iter = 5, response_seed = 1, ...))

  repeats <- ci_of(rdata, c(1:6, 1:6))
  averaged <- ci_of(rdata, c(1:6, 1:6), participants = rep(c("a", "b"), each = 6))
  for (m in list(score(repeats), score(repeats, reference_stimuli = 1:3),
                 score(averaged), score(averaged, reference_stimuli = 1:3))) {
    expect_length(m, 1)
    expect_match(m, "none is defined for this design")
    expect_no_match(m, "reference_stimuli", fixed = TRUE)
  }
  expect_match(score(averaged), "averages 2 participants")
})

test_that("the guard applies to independent bases", {
  tmp <- withr::local_tempdir()
  rdata <- make_independent_fixture(tmp)
  ci <- quiet(generateCI(c(1, 2, 3), c(1, -1, 1), "second", rdata, save_as_png = FALSE))

  m <- messages_of(computeInfoVal2IFC(ci, rdata, iter = 5, response_seed = 1, baseimage = "second"))
  expect_match(m, "built from 3 of the 12 saved stimuli")
  expect_length(messages_of(computeInfoVal2IFC(ci, rdata, iter = 5, response_seed = 1,
                                               baseimage = "second", reference_stimuli = 1:3)), 0)
})
