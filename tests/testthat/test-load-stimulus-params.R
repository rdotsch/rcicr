# Direct tests for loadStimulusParams(), which reads a stimulus .Rdata file in a
# frame of its own.
#
# The validation errors are also reached end to end through generateCI() by
# test-error-paths.R and test-rdata-writer-note.R, and those stay the guard for
# the message a user actually sees. What is only testable here is the frame:
# that the file's contents come back as a return value and land nowhere else.

good_rdata <- function(dir) {
  make_fixture_rdata(dir, img_size = 32, n_trials = 6) # nolint: object_usage_linter.
}

# Rewrites a fixture with `edit` applied to the loaded environment.
rewrite_rdata <- function(rdata, edit) {
  e <- new.env()
  load(rdata, envir = e)
  edit(e)
  save(list = ls(e), file = rdata, envir = e)
  rdata
}

test_that("loadStimulusParams returns the four objects the pipeline uses", {
  skip_if_not_installed("withr")

  out <- rcicr:::loadStimulusParams(good_rdata(withr::local_tempdir()))

  expect_named(out, c("p", "base_faces", "stimuli_params", "img_size"))
  expect_type(out$p, "list")
  expect_equal(out$img_size, 32)
  expect_true(is.matrix(out$stimuli_params$base))
  expect_true(!is.null(out$base_faces$base))
})

test_that("loadStimulusParams leaves nothing behind in its caller's frame", {
  skip_if_not_installed("withr")

  # The whole point of the extraction. The fixture carries a `sigma` field --
  # the one that used to capture generateCI()'s z-map blur argument -- so if the
  # load leaked into this frame, `sigma` would be bound here afterwards.
  rdata <- good_rdata(withr::local_tempdir())
  e <- new.env()
  load(rdata, envir = e)
  expect_equal(e$sigma, 25)

  # Called from a frame of known contents, so "nothing was added" is checkable
  # rather than a diff against this test's own bookkeeping.
  probe <- function(rdata) {
    sigma <- 3
    out <- rcicr:::loadStimulusParams(rdata)
    list(sigma = sigma, bound = sort(ls()))
  }
  res <- probe(rdata)

  expect_equal(res$sigma, 3)
  expect_equal(res$bound, c("out", "rdata", "sigma"))
})

test_that("loadStimulusParams survives a file that carries its own `rdata`", {
  skip_if_not_installed("withr")

  # The one name that could still collide, because it is the helper's own
  # formal: an older generateReferenceDistribution2IFC() saved its `rdata`
  # argument into the file, so files in the wild really do carry this name
  # (test-load-argument-guards.R:57-59 pins the same case for
  # computeInfoVal2IFC). Loading into a dedicated environment keeps the formal
  # out of reach entirely. Loading into the frame would not fail *today* --
  # nothing reads `rdata` after the load -- which is exactly why this is worth
  # pinning: it is a hazard that only becomes a bug when someone adds such a
  # read, and by then the cause is far from the symptom.
  dir <- withr::local_tempdir()
  rdata <- good_rdata(dir)
  mutate_rdata(rdata, rdata = file.path(dir, "does-not-exist.Rdata")) # nolint: object_usage_linter.

  out <- rcicr:::loadStimulusParams(rdata)

  expect_named(out, c("p", "base_faces", "stimuli_params", "img_size"))
  expect_equal(out$img_size, 32)
})

test_that("loadStimulusParams converts a pre-0.3.3 s into p", {
  skip_if_not_installed("withr")

  # Files written before 0.3.3 carry `s` (sinusoids/sinIdx) instead of `p`.
  rdata <- rewrite_rdata(good_rdata(withr::local_tempdir()), function(e) {
    p <- get("p", envir = e)
    assign("s", list(sinusoids = p$patches, sinIdx = p$patchIdx), envir = e)
    rm("p", envir = e)
  })

  out <- rcicr:::loadStimulusParams(rdata)

  expect_equal(out$p$noise_type, "sinusoid")
  expect_false(is.null(out$p$patches))
  expect_false(is.null(out$p$patchIdx))
})

test_that("loadStimulusParams lets s win over p, as loading in-frame did", {
  skip_if_not_installed("withr")

  # Not a preference: the in-frame version ran the conversion whenever `s`
  # existed, overwriting any `p` the file also held. Pinned so the precedence
  # cannot drift now that the two live in one expression.
  rdata <- rewrite_rdata(good_rdata(withr::local_tempdir()), function(e) {
    assign("s", list(sinusoids = "from-s", sinIdx = "idx-from-s"), envir = e)
  })

  out <- rcicr:::loadStimulusParams(rdata)

  expect_equal(out$p$patches, "from-s")
  expect_equal(out$p$noise_type, "sinusoid")
})

test_that("loadStimulusParams names each missing object, with the writer note", {
  skip_if_not_installed("withr")

  drop_field <- function(name) {
    rdata <- rewrite_rdata(good_rdata(withr::local_tempdir()),
                           function(e) rm(list = name, envir = e))
    tryCatch(rcicr:::loadStimulusParams(rdata), error = conditionMessage)
  }

  expect_match(drop_field("base_faces"), "did not contain base_faces", fixed = TRUE)
  expect_match(drop_field("stimuli_params"), "did not contain stimuli_params", fixed = TRUE)
  expect_match(drop_field("img_size"), "did not contain img_size", fixed = TRUE)
  expect_match(drop_field("p"), "did not contain s or p variable", fixed = TRUE)

  # The suffix rdataWriterNote() adds, which is only correct if it was handed
  # the frame the file was loaded into rather than the helper's own.
  expect_match(drop_field("img_size"),
               paste("written by rcicr", utils::packageVersion("rcicr")),
               fixed = TRUE)
})
