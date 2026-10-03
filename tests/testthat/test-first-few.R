# The six messages that list the first few offending items, pinned in full,
# for a list that fits and one that is cut with "...".

msg <- function(expr) tryCatch({ expr; NA_character_ }, error = conditionMessage)

test_that("trial-input messages list up to five trials, then ...", {
  coerce <- function(...) rcicr:::coerceTrialVectors(...)
  expect_identical(msg(coerce(1:8, rep(1, 8), c(NA, "a", NA, "b", "b", "b", NA, "a"))),
                   paste0("participants has no ID for 3 of 8 trials (trials 1, 3, 7). Give every ",
                          "trial an ID, or remove those trials from stimuli, responses and ",
                          "participants alike."))
  expect_identical(msg(coerce(1:8, rep(1, 8), c(rep(NA, 7), "a"))),
                   paste0("participants has no ID for 7 of 8 trials (trials 1, 2, 3, 4, 5, ...). ",
                          "Give every trial an ID, or remove those trials from stimuli, responses ",
                          "and participants alike."))
  expect_identical(msg(coerce(1:8, c(NA, rep(1, 7)), NA)),
                   paste0("responses has no finite value for 1 of 8 trials (trial 1). Give every ",
                          "trial a response, or remove those trials from stimuli, responses and ",
                          "participants alike."))
  expect_identical(msg(coerce(1:8, c(rep(NA, 6), 1, 1), NA)),
                   paste0("responses has no finite value for 6 of 8 trials (trials 1, 2, 3, 4, 5, ...). ",
                          "Give every trial a response, or remove those trials from stimuli, ",
                          "responses and participants alike."))
})

test_that("batch messages list up to five groups, then ...", {
  data <- data.frame(g = rep(letters[1:7], each = 2), pid = 1, resp = 1)
  some <- data
  some$pid[c(1, 5)] <- NA
  all <- data
  all$pid[] <- NA
  expect_identical(msg(rcicr:::requireParticipantIds(some, "g", "pid")),
                   paste0("2 rows have no participant ID in column pid (g a, c). Give every trial ",
                          "an ID, or remove the trials without one."))
  expect_identical(msg(rcicr:::requireParticipantIds(all, "g", "pid")),
                   paste0("14 rows have no participant ID in column pid (g a, b, c, d, e, ...). ",
                          "Give every trial an ID, or remove the trials without one."))
  some$resp[3] <- NA
  all$resp[] <- NA
  expect_identical(msg(rcicr:::requireBatchResponses(some, "g", "resp")),
                   paste0("1 rows have no finite response in column resp (g b). Give every trial ",
                          "a response, or remove those trials."))
  expect_identical(msg(rcicr:::requireBatchResponses(all, "g", "resp")),
                   paste0("14 rows have no finite response in column resp (g a, b, c, d, e, ...). ",
                          "Give every trial a response, or remove those trials."))
})

test_that("the batch trial-design message names up to five CIs, then ...", {
  averaged <- list(averaged = TRUE)
  three <- capture_messages(rcicr:::reportBatchTrialDesign(list(averaged, NULL, averaged, averaged),
                                                           c("p1", "p2", "p3", "p4")))
  expect_match(three, "^3 of the 4 classification images \\(p1, p3, p4\\) average")
  seven <- capture_messages(rcicr:::reportBatchTrialDesign(rep(list(averaged), 7), paste0("p", 1:7)))
  expect_match(seven, "^7 of the 7 classification images \\(p1, p2, p3, p4, p5, \\.\\.\\.\\) average")
})

test_that("repeated reference_stimuli are listed sorted, up to five, then ...", {
  expect_identical(msg(rcicr:::canonicalReferenceStimuli(c(4, 2, 4, 2, 9), 10)),
                   paste0("reference_stimuli lists stimulus 2, 4 more than once. A reference ",
                          "assumes one response per stimulus from one responder, so it cannot ",
                          "describe repeated presentations or several participants."))
  expect_identical(msg(rcicr:::canonicalReferenceStimuli(rep(c(7, 1, 6, 2, 5, 3), 2), 10)),
                   paste0("reference_stimuli lists stimulus 1, 2, 3, 5, 6, ... more than once. A ",
                          "reference assumes one response per stimulus from one responder, so it ",
                          "cannot describe repeated presentations or several participants."))
})

test_that("the PNG collision message lists up to three files, then ...", {
  dir <- withr::local_tempdir()
  paths <- file.path(dir, paste0("s", 1:5, ".png"))
  file.create(paths[1])
  expect_match(msg(rcicr:::reserveStimulusPngs(paths)),
               paste0("^1 of the 5 PNG files this call would write already exist \\(",
                      paths[1], "\\), and"))
  file.create(paths)
  expect_match(msg(rcicr:::reserveStimulusPngs(paths)),
               paste0("exist \\(", paste(paths[1:3], collapse = ", "), ", \\.\\.\\.\\), and"))
})

test_that("firstFew cuts only past n, not at n", {
  expect_identical(rcicr:::firstFew(1:5), "1, 2, 3, 4, 5")
  expect_identical(rcicr:::firstFew(1:6), "1, 2, 3, 4, 5, ...")
  expect_identical(rcicr:::firstFew(1:3, 3), "1, 2, 3")
  expect_identical(rcicr:::firstFew(1:4, 3), "1, 2, 3, ...")
  expect_identical(rcicr:::firstFew(character()), "")
})
