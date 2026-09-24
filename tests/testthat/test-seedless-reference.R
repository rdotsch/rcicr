# A stimulus file without a seed (absent or NULL) has no stream for the default
# reference to replay; set.seed(NULL) reseeds from the clock (#334).

make_seedless <- function(dir, kind = c("removed", "null"), independent = FALSE) {
  kind <- match.arg(kind)
  bases <- list(a = file.path(dir, "a.png"), b = file.path(dir, "b.png"))
  make_square_png(bases$a, size = 32, seed = 1) # nolint: object_usage_linter.
  make_square_png(bases$b, size = 32, seed = 2) # nolint: object_usage_linter.
  if (!independent) bases <- bases["a"]
  seed <- if (kind == "null") NULL else 1
  quietly(
    generateStimuli2IFC(bases, n_trials = 6, img_size = 32, stimulus_path = dir, seed = seed,
                        ncores = 1, nscales = 1, save_as_png = FALSE,
                        use_same_parameters = !independent)
  )
  path <- list.files(dir, "\\.Rdata$", full.names = TRUE)
  if (kind == "removed") edit_rdata(path, function(e) rm("seed", envir = e))
  path
}

edit_rdata <- function(path, f) {
  e <- new.env()
  load(path, envir = e)
  f(e)
  save(list = ls(e, all.names = TRUE), file = path, envir = e)
}

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr))
  out
}

messages_of <- function(expr) {
  seen <- character()
  value <- withCallingHandlers(quietly(expr), message = function(m) {
    seen <<- c(seen, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(value = value, messages = seen)
}

target <- function() list(ci = matrix(withr::with_seed(2, rnorm(32 * 32)), 32, 32))

test_that("every default reference path stops on a seedless file before building noise", {
  local_mocked_bindings(referenceNoise = function(...) stop("reached referenceNoise"))
  for (kind in c("removed", "null")) {
    dir <- withr::local_tempdir()
    rd <- make_seedless(dir, kind)
    expect_error(quietly(generateReferenceDistribution2IFC(rd, iter = 3, ncores = 1)),
                 "no stimulus seed to replay")
    expect_error(quietly(computeInfoVal2IFC(target(), rd, iter = 3)), "no stimulus seed to replay")
    expect_error(quietly(computeInfoVal2IFC(target(), rd, iter = 3, force_gen_ref_dist = TRUE)),
                 "no stimulus seed to replay")

    dir <- withr::local_tempdir()
    rd <- make_seedless(dir, kind, independent = TRUE)
    expect_error(quietly(computeInfoVal2IFC(target(), rd, iter = 3, baseimage = "a")),
                 'baseimage = "a"', fixed = TRUE)
    expect_error(quietly(generateReferenceDistribution2IFC(rd, iter = 3, ncores = 1, baseimage = "b")),
                 'baseimage = "b"', fixed = TRUE)
  }
})

test_that("a stale stored reference on a seedless file stops instead of being rebuilt", {
  for (independent in c(FALSE, TRUE)) {
    dir <- withr::local_tempdir()
    rd <- make_seedless(dir, "null", independent = independent)
    stale <- withr::with_seed(3, runif(5))
    edit_rdata(rd, function(e) {
      if (independent) {
        e$reference_norms_by_base <- list(a = list(norms = stale))
      } else {
        e$reference_norms <- stale
      }
    })
    base <- if (independent) "a" else NULL
    expect_error(quietly(computeInfoVal2IFC(target(), rd, iter = 3, baseimage = base)),
                 "no stimulus seed to replay")
  }
})

test_that("the call the error gives runs as written and makes the file reproducible", {
  for (independent in c(FALSE, TRUE)) {
    dir <- withr::local_tempdir()
    rd <- make_seedless(dir, "null", independent = independent)
    base <- if (independent) "a" else NULL
    err <- tryCatch(quietly(computeInfoVal2IFC(target(), rd, iter = 3, baseimage = base)),
                    error = conditionMessage)
    call <- regmatches(err, regexpr("generateReferenceDistribution2IFC\\([^)]*\\)", err))
    expect_length(call, 1)
    call <- sub("<n>", "7", call, fixed = TRUE)
    call <- sub(")$", ", iter = 3, ncores = 1)", call)
    quietly(eval(parse(text = call)))

    first <- messages_of(computeInfoVal2IFC(target(), rd, iter = 3, baseimage = base))
    second <- messages_of(computeInfoVal2IFC(target(), rd, iter = 3, baseimage = base))
    expect_identical(first$value, second$value)
    expect_false(any(grepl("cannot be regenerated", first$messages)))
  }
})

test_that("a read-only seedless file is scored reproducibly with response_seed, unwritten", {
  dir <- withr::local_tempdir()
  rd <- make_seedless(dir, "removed")
  before <- tools::md5sum(rd)
  local_mocked_bindings(writableFile = function(path) FALSE)
  a <- quietly(computeInfoVal2IFC(target(), rd, iter = 3, response_seed = 5))
  b <- quietly(computeInfoVal2IFC(target(), rd, iter = 3, response_seed = 5))
  expect_identical(a, b)
  expect_identical(tools::md5sum(rd), before)
})

test_that("a vouched stored reference on a seedless file is used unchanged, with a message", {
  dir <- withr::local_tempdir()
  rd <- make_fixture_rdata(dir, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  quietly(generateReferenceDistribution2IFC(rd, iter = 5, ncores = 1))
  with_seed <- messages_of(computeInfoVal2IFC(target(), rd, iter = 5))
  expect_false(any(grepl("cannot be regenerated", with_seed$messages)))

  edit_rdata(rd, function(e) e$seed <- NULL)
  without <- messages_of(computeInfoVal2IFC(target(), rd, iter = 5))
  expect_identical(without$value, with_seed$value)
  expect_true(any(grepl("cannot be regenerated", without$messages)))

  # Not quietly(): suppressWarnings() would hide a warning before warn = 2 made it an error.
  withr::local_options(warn = 2)
  utils::capture.output(strict <- suppressMessages(computeInfoVal2IFC(target(), rd, iter = 5)))
  expect_identical(strict, with_seed$value)
})
