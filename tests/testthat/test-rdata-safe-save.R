# saveRdataSafely() and loadRdata(): a save into a stimulus file must never
# leave it unloadable, and a file damaged by a killed save must point to its
# backup from every entry point (#333).

stim_file <- function(dir, value = 1) {
  path <- file.path(dir, "stim.Rdata")
  x <- value
  save(x, file = path)
  path
}

loaded_x <- function(path) {
  e <- new.env()
  load(path, envir = e)
  e$x
}

save_x <- function(path, value) {
  e <- new.env()
  e$x <- value
  saveRdataSafely("x", path, e)
}

sidecars <- function(dir) {
  f <- list.files(dir, all.files = TRUE, no.. = TRUE)
  f[grepl("rcicr-(backup|staging)", f)]
}

test_that("a successful save matches a plain save() and leaves no backup", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir)
  save_x(path, 2)
  expect_identical(loaded_x(path), 2)
  expect_identical(sidecars(dir), character(0))
})

test_that("a save interrupted mid-write leaves the original intact", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  mode_before <- file.mode(path)
  local_mocked_bindings(writeRdata = function(names, file, envir) {
    writeLines("partial", file)
    stop("interrupted")
  })
  expect_error(save_x(path, 2), "interrupted")
  expect_identical(loaded_x(path), 1)
  expect_identical(file.mode(path), mode_before)
  expect_identical(sidecars(dir), character(0))
})

test_that("an interruption while backing up restores nothing and removes the staging copy", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  restored <- FALSE
  local_mocked_bindings(
    copyInto = function(from, to) stop("interrupted"),
    restoreFromBackup = function(file, backup) restored <<- TRUE
  )
  expect_error(save_x(path, 2), "interrupted")
  expect_false(restored)
  expect_identical(loaded_x(path), 1)
  expect_identical(sidecars(dir), character(0))
})

test_that("a failed backup copy stops before the original is touched", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  local_mocked_bindings(copyInto = function(from, to) FALSE)
  expect_error(save_x(path, 2), "unchanged")
  expect_identical(loaded_x(path), 1)
  expect_identical(sidecars(dir), character(0))
})

test_that("a refused publishing rename stops before the original is touched", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  local_mocked_bindings(renameFile = function(from, to) FALSE)
  expect_error(save_x(path, 2), "unchanged")
  expect_identical(loaded_x(path), 1)
  expect_identical(sidecars(dir), character(0))
})

test_that("a staging file that cannot be created stops, unless the directory is read-only", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  local_mocked_bindings(createStaging = function(path) FALSE, writableDir = function(path) TRUE)
  expect_error(save_x(path, 2), "unchanged")
  expect_identical(loaded_x(path), 1)

  skip_on_os("windows")
  local_mocked_bindings(createStaging = function(path) FALSE, writableDir = function(path) FALSE)
  expect_warning(save_x(path, 3), "without a backup")
  expect_identical(loaded_x(path), 3)
})

test_that("a name too long for the staging suffix is saved without a backup", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  path <- file.path(dir, paste0(strrep("a", 229), ".Rdata"))
  x <- 1
  save(x, file = path)
  expect_warning(save_x(path, 2), "too long")
  expect_identical(loaded_x(path), 2)
})

test_that("the backup is owner-only and the original keeps its mode", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  path <- stim_file(dir)
  Sys.chmod(path, "660", use_umask = FALSE)
  seen <- NULL
  real_write <- writeRdata
  local_mocked_bindings(writeRdata = function(names, file, envir) {
    seen <<- file.mode(rdataBackupPath(file))
    real_write(names, file, envir)
  })
  save_x(path, 2)
  expect_identical(seen, as.octmode("600"))
  expect_identical(file.mode(path), as.octmode("660"))
})

test_that("a chmod that does not take stops before any data is copied", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  copied <- FALSE
  local_mocked_bindings(
    makePrivate = function(path) Sys.chmod(path, "644", use_umask = FALSE),
    copyInto = function(from, to) copied <<- TRUE
  )
  expect_error(save_x(path, 2), "owner")
  expect_false(copied)
  expect_identical(sidecars(dir), character(0))
})

test_that("a committed save is kept when its backup cannot be deleted", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  local_mocked_bindings(removeFile = function(path) FALSE)
  expect_warning(save_x(path, 2), "could not be deleted")
  expect_identical(loaded_x(path), 2)
})

test_that("a leftover backup beside a file that loads is removed and the save proceeds", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  file.copy(path, rdataBackupPath(path))
  expect_message(save_x(path, 2), "Removed")
  expect_identical(loaded_x(path), 2)
  expect_identical(sidecars(dir), character(0))
})

test_that("a leftover backup beside a damaged file stops with the restore command", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir, value = 1)
  file.copy(path, rdataBackupPath(path))
  writeLines("partial", path)
  expect_error(save_x(path, 2), "file.copy", fixed = TRUE)
  expect_identical(readLines(path), "partial")
  expect_identical(loaded_x(rdataBackupPath(path)), 1)
})

test_that("leftover staging files are reported and left in place", {
  dir <- withr::local_tempdir()
  path <- stim_file(dir)
  leftover <- file.path(dir, "stim.Rdata.rcicr-staging-abc123")
  writeLines("partial", leftover)
  expect_warning(save_x(path, 2), "stim.Rdata.rcicr-staging-abc123", fixed = TRUE)
  expect_true(file.exists(leftover))
})

test_that("every reader of a damaged stimulus file points to its backup", {
  dir <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(dir, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  file.copy(rdata_path, rdataBackupPath(rdata_path))
  writeLines("partial", rdata_path)

  # load() warns about the damaged header before it errors; the error is the point.
  suppressWarnings({
    expect_error(
      generateCI(stimuli = 1:6, responses = rep(c(1, -1), 3), baseimage = "base",
                 rdata = rdata_path, save_as_png = FALSE),
      "rcicr-backup", fixed = TRUE
    )
    expect_error(
      generateReferenceDistribution2IFC(rdata_path, iter = 3, ncores = 1),
      "rcicr-backup", fixed = TRUE
    )
    expect_error(
      computeInfoVal2IFC(list(ci = matrix(0, 32, 32)), rdata_path),
      "rcicr-backup", fixed = TRUE
    )
  })
})

test_that("a reference distribution is still saved into the stimulus file", {
  dir <- withr::local_tempdir()
  rdata_path <- make_fixture_rdata(dir, img_size = 32, n_trials = 6, nscales = 1, seed = 1)
  norms <- suppressWarnings(generateReferenceDistribution2IFC(rdata_path, iter = 3, ncores = 1))
  e <- new.env()
  load(rdata_path, envir = e)
  expect_identical(e$reference_norms, norms)
  expect_identical(sidecars(dir), character(0))
})
