# Every path through resolveReferenceNorms(), with what it prints, what it
# messages, and what it asks the simulation for, pinned exactly (#390).

run <- function(entry, ..., writable = TRUE, simulated = c(1, 2, 3)) {
  calls <- list()
  testthat::local_mocked_bindings(
    generateReferenceDistribution2IFC = function(rdata, iter, response_seed, save_rdata, ...) {
      calls[[length(calls) + 1]] <<- list(iter = iter, response_seed = response_seed,
                                          save_rdata = save_rdata)
      simulated
    },
    writableFile = function(path) writable
  )
  messages <- character()
  out <- utils::capture.output(norms <- withCallingHandlers(
    rcicr:::resolveReferenceNorms(entry, "stim.Rdata", iter = 50, ...),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  list(norms = norms, out = out, messages = messages, calls = calls)
}

vouched <- function(norms) {
  list(norms = norms, response_seed = NULL, source = "saved_noise",
       fingerprint = rcicr:::referenceSnapshot(norms))
}

test_that("a stored reference is used as stored", {
  r <- run(vouched(4:6), force_gen_ref_dist = FALSE, response_seed = NULL)
  expect_identical(r$norms, 4:6)
  expect_identical(r$out, "Using reference distribution found in rdata file.")
  expect_identical(r$messages, character())
  expect_length(r$calls, 0)
})

test_that("a stored reference in a seedless file says how to store a reproducible one", {
  r <- run(vouched(4:6), force_gen_ref_dist = FALSE, response_seed = NULL, seedless = TRUE,
           baseimage = "b", reference_stimuli = 1:3)
  expect_identical(r$out, "Using reference distribution found in rdata file.")
  expect_identical(r$messages, paste0(
    "stim.Rdata has no stimulus seed, so its stored reference distribution cannot be ",
    "regenerated from it. It is used as stored. Store a reproducible one with ",
    "generateReferenceDistribution2IFC(\"stim.Rdata\", response_seed = <n>, baseimage = \"b\", ",
    "reference_stimuli = <the same stimuli>); later computeInfoVal2IFC() calls reuse it. If the ",
    "file cannot be written, pass response_seed = <n> to computeInfoVal2IFC() instead.\n"
  ))
})

test_that("a missing reference is simulated and saved", {
  r <- run(NULL, force_gen_ref_dist = FALSE, response_seed = NULL)
  expect_identical(r$norms, c(1, 2, 3))
  expect_identical(r$out, "The reference distribution has been saved to the .Rdata file for reuse.")
  expect_identical(r$calls, list(list(iter = 50, response_seed = NULL, save_rdata = TRUE)))
})

test_that("a read-only file is scored without storing, naming the base when there is one", {
  r <- run(NULL, force_gen_ref_dist = FALSE, response_seed = NULL, writable = FALSE)
  expect_identical(r$out, paste0("Built the reference from saved noise, but stim.Rdata is not ",
                                 "writable, so the values were used without being stored. ",
                                 "The next call will build them again."))
  expect_identical(r$calls[[1]]$save_rdata, FALSE)
  r <- run(NULL, force_gen_ref_dist = FALSE, response_seed = NULL, writable = FALSE,
           baseimage = "b")
  expect_match(r$out, "^Built the reference for baseimage b from saved noise", all = FALSE)
})

test_that("a seeded draw is never saved, and says so", {
  r <- run(vouched(4:6), force_gen_ref_dist = FALSE, response_seed = 5, writable = FALSE)
  expect_identical(r$out, paste0("Reference distribution simulated with response_seed = 5. ",
                                 "This independent draw has deliberately not been saved."))
  expect_identical(r$calls, list(list(iter = 50, response_seed = 5, save_rdata = FALSE)))
  expect_identical(r$messages, character())
})

test_that("forcing replaces a stored reference and saves the new one", {
  r <- run(vouched(4:6), force_gen_ref_dist = TRUE, response_seed = NULL)
  expect_identical(r$norms, c(1, 2, 3))
  expect_identical(r$out, "The reference distribution has been saved to the .Rdata file for reuse.")
  expect_identical(r$messages, character())
})

test_that("an unvouched reference is rebuilt at its own size, and a change is reported", {
  stale <- list(norms = c(7, 8, 9, 10), response_seed = NULL)
  r <- run(stale, force_gen_ref_dist = FALSE, response_seed = NULL)
  expect_identical(r$calls, list(list(iter = 4L, response_seed = NULL, save_rdata = TRUE)))
  expect_identical(r$messages, paste0(
    "This stimulus file carried a reference distribution from before rcicr built references ",
    "from the saved noise, and rebuilding it from the saved noise gave different values. The ",
    "InfoVal this returns supersedes any computed from this file before.\n"
  ))
  same <- run(stale, force_gen_ref_dist = FALSE, response_seed = NULL, simulated = c(7, 8, 9, 10))
  expect_identical(same$messages, character())
})

test_that("the seedless advice names each restriction the stored reference has", {
  advice <- function(...) {
    run(vouched(4:6), force_gen_ref_dist = FALSE, response_seed = NULL, seedless = TRUE, ...)$messages
  }
  head <- paste0("stim.Rdata has no stimulus seed, so its stored reference distribution cannot be ",
                 "regenerated from it. It is used as stored. Store a reproducible one with ",
                 "generateReferenceDistribution2IFC(\"stim.Rdata\", response_seed = <n>")
  tail <- paste0("); later computeInfoVal2IFC() calls reuse it. If the file cannot be written, ",
                 "pass response_seed = <n> to computeInfoVal2IFC() instead.\n")
  expect_identical(advice(), paste0(head, tail))
  expect_identical(advice(baseimage = "b"), paste0(head, ", baseimage = \"b\"", tail))
  expect_identical(advice(reference_stimuli = 1:3),
                   paste0(head, ", reference_stimuli = <the same stimuli>", tail))
  expect_identical(advice(masked = c(TRUE, FALSE, FALSE, FALSE)),
                   paste0(head, ", mask = <the same mask>", tail))
  expect_identical(advice(baseimage = "b", reference_stimuli = 1:3, masked = c(TRUE, FALSE, FALSE, FALSE)),
                   paste0(head, ", baseimage = \"b\", reference_stimuli = <the same stimuli>, ",
                          "mask = <the same mask>", tail))
})
