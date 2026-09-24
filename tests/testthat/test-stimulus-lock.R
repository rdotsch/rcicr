# generateStimuli2IFC() never overwrites a stimulus .Rdata file, and reserves
# its name with a lock keyed on seed and minute (#338).

fixed_minute <- as.POSIXct("2026-09-24 10:00:00", tz = "UTC")

generate <- function(dir, label = "rcic", seed = 1, ...) {
  base <- file.path(tempdir(), "lock-base.png")
  if (!file.exists(base)) make_square_png(base, size = 16, seed = 1) # nolint: object_usage_linter.
  suppressWarnings(utils::capture.output(
    generateStimuli2IFC(list(face = base), n_trials = 2, img_size = 16, stimulus_path = dir,
                        label = label, seed = seed, ncores = 1, nscales = 1, ...)
  ))
  invisible(NULL)
}

snapshot <- function(dir) {
  files <- list.files(dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  stats::setNames(unname(tools::md5sum(files)), basename(files))
}

test_that("a second call for the same file stops and leaves the first set untouched", {
  local_mocked_bindings(stimulusTime = function() fixed_minute)
  dir <- withr::local_tempdir()
  generate(dir)
  first <- snapshot(dir)
  expect_length(grep("\\.Rdata$", names(first)), 1)

  expect_error(generate(dir), "already exists, and a stimulus file is never overwritten")
  expect_identical(snapshot(dir), first)
})

test_that("a different label, or save_rdata = FALSE, does not stop", {
  local_mocked_bindings(stimulusTime = function() fixed_minute)
  dir <- withr::local_tempdir()
  generate(dir)
  expect_no_error(generate(dir, label = "other"))
  expect_no_error(generate(dir, save_rdata = FALSE))
  expect_length(list.files(dir, "\\.Rdata$"), 2)
})

test_that("a held lock stops any call with that seed and minute, however the target is spelled", {
  local_mocked_bindings(stimulusTime = function() fixed_minute)
  dir <- withr::local_tempdir()
  dir.create(file.path(dir, "sub"))
  lock <- rcicr:::stimulusLockPath(rcicr:::stimulusRdataPath(dir, "rcic", 1, fixed_minute), 1, fixed_minute)
  dir.create(lock)

  spellings <- list(
    list(path = dir, label = "trial"),
    list(path = dir, label = "Trial"),
    list(path = dir, label = "é"),
    list(path = dir, label = "é"),
    list(path = file.path(dir, "sub", ".."), label = "rcic")
  )
  for (s in spellings) {
    expect_error(generate(s$path, label = s$label),
                 "Delete it only after confirming that no generateStimuli2IFC() call is still writing",
                 fixed = TRUE)
  }
  expect_identical(setdiff(list.files(dir, all.files = TRUE, no.. = TRUE), c("sub", basename(lock))),
                   character(0))
})

test_that("a call that fails after taking the lock leaves neither the lock nor a file", {
  local_mocked_bindings(generateNoiseImage = function(...) stop("interrupted"))
  dir <- withr::local_tempdir()
  expect_error(generate(dir), "interrupted")
  expect_identical(list.files(dir, all.files = TRUE, no.. = TRUE), character(0))
})

test_that("while a call runs, no list.files() pattern for the output finds the lock", {
  real <- rcicr:::generateNoiseImage
  dir <- withr::local_tempdir()
  during <- NULL
  local_mocked_bindings(generateNoiseImage = function(...) {
    during <<- list(all = list.files(dir, all.files = TRUE, no.. = TRUE),
                    rdata = list.files(dir, "Rdata", ignore.case = TRUE, all.files = TRUE))
    real(...)
  })
  generate(dir, label = "myRdata", save_as_png = FALSE)
  expect_true(any(startsWith(during$all, ".rcicr-lock-")))
  expect_identical(during$rdata, character(0))
})

test_that("the file is named for the minute the call started", {
  readings <- list(fixed_minute, fixed_minute + 3600)
  local_mocked_bindings(stimulusTime = function() {
    now <- readings[[1]]
    readings <<- readings[-1]
    now
  })
  dir <- withr::local_tempdir()
  generate(dir, save_as_png = FALSE)
  expect_identical(list.files(dir, "\\.Rdata$"),
                   basename(rcicr:::stimulusRdataPath(dir, "rcic", 1, fixed_minute)))
})

test_that("the lock's folder follows the host's path rules, even for text dirname() cannot translate", {
  for (p in c("/a/b/x.Rdata", "/x.Rdata", "a/x.Rdata", "x.Rdata", "./x.Rdata")) {
    expect_identical(rcicr:::targetDir(p), dirname(p))
  }
  expect_identical(rcicr:::targetDir("/a/\u00e9_seed.Rdata"), "/a")
  if (.Platform$OS.type == "windows") {
    expect_identical(rcicr:::targetDir("C:/x.Rdata"), "C:/")
    expect_identical(rcicr:::targetDir("C:\\a\\x.Rdata"), "C:\\a")
  } else {
    # A backslash is part of a file name here, not a separator.
    expect_identical(rcicr:::targetDir("/a/b\\c_seed.Rdata"), "/a")
    expect_identical(rcicr:::targetDir("C:/x.Rdata"), "C:")
  }
})
