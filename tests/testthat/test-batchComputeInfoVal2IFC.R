# batchComputeInfoVal2IFC() (#85) returns what a named loop over
# computeInfoVal2IFC() returns, while resolving each distinct reference once.

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

copy_of <- function(rdata) {
  dest <- file.path(withr::local_tempdir(.local_envir = parent.frame()), basename(rdata))
  file.copy(rdata, dest)
  dest
}

loop_oracle <- function(cis, rdata, reference_stimuli = NULL, ...) {
  per_ci <- if (is.list(reference_stimuli)) reference_stimuli else rep(list(reference_stimuli), length(cis))
  stats::setNames(vapply(seq_along(cis), function(i) {
    quietly(computeInfoVal2IFC(cis[[i]], rdata, reference_stimuli = per_ci[[i]], ...))
  }, numeric(1)), names(cis))
}

participant_cis <- function(rdata, n_trials = 12, per = 4) {
  data <- data.frame(participant = rep(paste0("p", seq_len(n_trials / per)), each = per),
                     stimulus = seq_len(n_trials),
                     response = rep(c(1, -1, -1, 1), length.out = n_trials))
  quietly(batchGenerateCI2IFC(data, "participant", "stimulus", "response", "base", rdata,
                              save_as_png = FALSE, scaling = "none"))
}

own_stimuli <- function(cis) lapply(cis, function(ci) attr(ci, "trial_design")$stimuli)

test_that("the result is identical to a named computeInfoVal2IFC() loop", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  cis <- participant_cis(rdata)
  expect_named(cis, paste0("base_participant_p", 1:3))
  full <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  whole <- list(a = quietly(generateCI(1:12, rep(c(1, -1), 6), "base", full, save_as_png = FALSE)),
                b = quietly(generateCI(1:12, rep(c(-1, 1, 1), 4), "base", full, save_as_png = FALSE)))

  cases <- list(
    list(cis = whole, args = list()),
    list(cis = unname(whole), args = list()),
    list(cis = cis, args = list(reference_stimuli = 1:8)),
    list(cis = cis, args = list(reference_stimuli = own_stimuli(cis))),
    list(cis = cis, args = list(reference_stimuli = list(NULL, c(8, 6, 7, 5), 5:8))),
    list(cis = whole, args = list(response_seed = 4)),
    list(cis = whole, args = list(force_gen_ref_dist = TRUE)),
    list(cis = cis, args = list(reference_stimuli = own_stimuli(cis), reference_method = "images"))
  )
  for (k in seq_along(cases)) {
    source_file <- if (identical(cases[[k]]$cis, cis)) rdata else full
    for_loop <- copy_of(source_file)
    for_batch <- copy_of(source_file)
    expected <- do.call(loop_oracle, c(list(cases[[k]]$cis, for_loop, iter = 30), cases[[k]]$args))
    actual <- quietly(do.call(batchComputeInfoVal2IFC,
                              c(list(cases[[k]]$cis, for_batch, iter = 30), cases[[k]]$args)))
    expect_identical(actual, expected, info = k)
    expect_gt(length(unique(actual)), 1)
  }
})

test_that("an independent-base file is scored per base, as the loop scores it", {
  rdata <- make_independent_fixture(withr::local_tempdir())
  cis <- lapply(list(x = rep(c(1, -1), 6), y = rep(c(1, 1, -1), 4)), function(r) {
    quietly(generateCI(1:12, r, "second", rdata, save_as_png = FALSE))
  })
  # Full-set and subset references for the same base, in one call.
  mixed <- list(NULL, 1:9)
  expect_identical(
    quietly(batchComputeInfoVal2IFC(c(cis, cis), copy_of(rdata), iter = 30, baseimage = "second",
                                    reference_stimuli = c(mixed, mixed))),
    loop_oracle(c(cis, cis), copy_of(rdata), reference_stimuli = c(mixed, mixed), iter = 30,
                baseimage = "second")
  )
  for (stimuli in list(NULL, 1:9)) {
    expected <- loop_oracle(cis, copy_of(rdata), reference_stimuli = stimuli, iter = 30,
                            baseimage = "second")
    actual <- quietly(batchComputeInfoVal2IFC(cis, copy_of(rdata), iter = 30, baseimage = "second",
                                              reference_stimuli = stimuli))
    expect_identical(actual, expected)
  }
})

test_that("a reference that cannot be stored is simulated once per distinct stimulus set", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  cis <- rep(participant_cis(rdata), length.out = 5)
  stimuli <- list(1:6, 7:12, 6:1, 1:6, c(12, 7:11))
  calls <- 0
  original <- generateReferenceDistribution2IFC
  local_mocked_bindings(
    writableFile = function(path) FALSE,
    generateReferenceDistribution2IFC = function(...) {
      calls <<- calls + 1
      original(...)
    }
  )
  batch <- quietly(batchComputeInfoVal2IFC(cis, rdata, iter = 30, reference_stimuli = stimuli))
  expect_identical(calls, 2)
  calls <- 0
  expect_identical(batch, loop_oracle(cis, rdata, reference_stimuli = stimuli, iter = 30))
  expect_identical(calls, 5)
})

test_that("the stimulus file is read the same number of times for any number of CIs", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  ci <- quietly(generateCI(1:12, rep(c(1, -1), 6), "base", rdata, save_as_png = FALSE))
  quietly(computeInfoVal2IFC(ci, rdata, iter = 30))
  loads <- 0
  original <- loadRdata
  local_mocked_bindings(loadRdata = function(...) {
    loads <<- loads + 1
    original(...)
  })
  counted <- function(n) {
    loads <<- 0
    quietly(batchComputeInfoVal2IFC(rep(list(ci), n), rdata, iter = 30))
    loads
  }
  expect_identical(counted(2), counted(20))
  loads <- 0
  loop_oracle(rep(list(ci), 20), rdata, iter = 30)
  expect_gt(loads, counted(20))
})

test_that("design messages are collected into one per kind, naming the CIs", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  cis <- participant_cis(rdata)
  pooled <- quietly(generateCI(rep(1:12, 2), rep(c(1, -1), 12), "base", rdata, save_as_png = FALSE))
  cis <- c(cis, list(pooled = pooled))
  seen <- character()
  utils::capture.output(withCallingHandlers(
    suppressWarnings(batchComputeInfoVal2IFC(cis, rdata, iter = 30,
                                             reference_stimuli = list(1:4, NULL, NULL, NULL))),
    message = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  design <- grep("classification images \\(", seen, value = TRUE)
  expect_length(design, 2)
  expect_match(design, "1 of the 4 classification images (pooled) average", fixed = TRUE, all = FALSE)
  expect_match(design, "2 of the 4 classification images (base_participant_p2, base_participant_p3) were built from different stimuli",
               fixed = TRUE, all = FALSE)

  seen <- character()
  utils::capture.output(withCallingHandlers(
    suppressWarnings(batchComputeInfoVal2IFC(cis[1:3], rdata, iter = 30,
                                             reference_stimuli = own_stimuli(cis[1:3]))),
    message = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  expect_length(grep("classification images \\(", seen), 0)
})

test_that("a single CI, a malformed batch and a mismatched reference_stimuli are refused", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 12)
  ci <- quietly(generateCI(1:12, rep(c(1, -1), 6), "base", rdata, save_as_png = FALSE))
  expect_error(batchComputeInfoVal2IFC(ci, rdata), "is a single classification image")
  expect_error(batchComputeInfoVal2IFC(list(), rdata), "must be a list of classification images")
  expect_error(batchComputeInfoVal2IFC(list(a = ci, b = matrix(1, 32, 32)), rdata),
               "Element b of target_cis is not a classification image")
  expect_error(batchComputeInfoVal2IFC(list(ci, ci), rdata, reference_stimuli = list(1:6)),
               "has 1 elements for 2 classification images")

  # A batch may name an element ci.
  named_ci <- list(ci = ci, control = ci)
  expect_identical(quietly(batchComputeInfoVal2IFC(named_ci, rdata, iter = 30)),
                   loop_oracle(named_ci, rdata, iter = 30))
})
