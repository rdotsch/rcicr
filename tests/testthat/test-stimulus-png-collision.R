# generateStimuli2IFC() never overwrites stimulus PNGs (#350). Their names carry
# no time, so a later-minute call into the same folder used to overwrite them
# while the earlier .Rdata, which described them, survived.

minute <- function(m) as.POSIXct(sprintf("2026-09-24 10:%02d:00", m), tz = "UTC")

base_png <- function() {
  base <- file.path(tempdir(), "png-collision-base.png")
  if (!file.exists(base)) make_square_png(base, size = 16, seed = 1) # nolint: object_usage_linter.
  base
}

generate <- function(dir, n_trials = 6, nscales = 1, label = "rcic", bases = list(face = base_png()),
                     ...) {
  suppressWarnings(utils::capture.output(
    generateStimuli2IFC(bases, n_trials = n_trials, img_size = 16, stimulus_path = dir,
                        label = label, seed = 1, ncores = 1, nscales = nscales, ...)
  ))
  invisible(NULL)
}

# Directories, such as a lock, are recorded by name; md5sum() reads files only.
snapshot <- function(dir) {
  paths <- list.files(dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  sums <- ifelse(dir.exists(paths), "<directory>", "")
  sums[!dir.exists(paths)] <- unname(tools::md5sum(paths[!dir.exists(paths)]))
  stats::setNames(sums, basename(paths))
}

test_that("a later-minute call that would overwrite the PNGs stops and changes nothing", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(stimulusTime = function() minute(0))
  generate(dir, n_trials = 6, nscales = 1)
  first <- snapshot(dir)
  expect_length(grep("\\.png$", names(first)), 12)

  local_mocked_bindings(stimulusTime = function() minute(1))
  expect_error(generate(dir, n_trials = 6, nscales = 2), "12 of the 12 PNG files")
  expect_identical(snapshot(dir), first)

  # It stops before the noise basis, the slow step, is built.
  local_mocked_bindings(generateNoisePattern = function(...) stop("basis built"))
  expect_error(generate(dir, n_trials = 6, nscales = 2), "12 of the 12 PNG files")
  local_mocked_bindings(generateNoisePattern = rcicr::generateNoisePattern)

  # A shorter run collides on the files it shares, and leaves the rest alone too.
  expect_error(generate(dir, n_trials = 4, nscales = 2), "8 of the 8 PNG files")
  expect_identical(snapshot(dir), first)
})

test_that("a different label or folder, or regenerating after deleting the set, still works", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(stimulusTime = function() minute(0))
  generate(dir)
  local_mocked_bindings(stimulusTime = function() minute(1))
  expect_no_error(generate(dir, label = "other"))
  expect_no_error(generate(withr::local_tempdir()))

  unlink(list.files(dir, pattern = "^rcic_", full.names = TRUE))
  expect_no_error(generate(dir, nscales = 2))
  expect_length(list.files(dir, pattern = "^rcic_face_.*\\.png$"), 12)
})

test_that("save_as_png = FALSE neither checks nor touches existing PNGs", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(stimulusTime = function() minute(0))
  generate(dir)
  pngs <- snapshot(dir)[grep("\\.png$", names(snapshot(dir)))]
  local_mocked_bindings(stimulusTime = function() minute(1))
  expect_no_error(generate(dir, save_as_png = FALSE))
  expect_identical(snapshot(dir)[names(pngs)], pngs)
})

test_that("a held PNG lock stops a same-seed call before anything is written", {
  dir <- withr::local_tempdir()
  lock <- rcicr:::acquirePngLock(dir, 1)
  before <- snapshot(dir)
  expect_error(generate(dir), "are reserved by")
  expect_identical(snapshot(dir), before)
  unlink(lock, recursive = TRUE)
  expect_no_error(generate(dir))
})

test_that("base labels the file system treats as one name stop the call, leaving nothing", {
  dir <- withr::local_tempdir()
  # A case-insensitive file system, modelled on any runner.
  local_mocked_bindings(pathTaken = function(path) {
    tolower(basename(path)) %in% tolower(list.files(dirname(path), all.files = TRUE))
  })
  expect_error(generate(dir, n_trials = 2, bases = list(face = base_png(), Face = base_png())),
               "is the same file, on this file system, as another PNG this call writes")
  expect_length(list.files(dir), 0)
})

test_that("on a case-insensitive file system, face and Face collide for real", {
  probe <- withr::local_tempdir()
  file.create(file.path(probe, "a"))
  skip_if_not(file.exists(file.path(probe, "A")), "this file system is case-sensitive")
  dir <- withr::local_tempdir()
  expect_error(generate(dir, n_trials = 2, bases = list(face = base_png(), Face = base_png())),
               "is the same file, on this file system, as another PNG this call writes")
  expect_length(list.files(dir), 0)
})

test_that("a call that fails part-way leaves no placeholder and no PNG behind", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(generateNoiseImage = function(...) stop("generation failed"))
  expect_error(generate(dir), "generation failed")
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
  local_mocked_bindings(generateNoiseImage = rcicr::generateNoiseImage)
})
