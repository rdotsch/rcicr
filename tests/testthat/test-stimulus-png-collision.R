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

test_that("a dangling symlink at a PNG's name is taken, and neither it nor its target is touched", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  outside <- file.path(withr::local_tempdir(), "target.png")
  link <- rcicr:::stimulusPngPath(dir, "rcic", "face", 1, 1, "ori")
  file.symlink(outside, link)
  expect_error(generate(dir), "1 of the 12 PNG files")
  expect_identical(Sys.readlink(link), outside)
  expect_false(file.exists(outside))
})

test_that("a call that fails part-way leaves no placeholder and no PNG behind", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(generateNoiseImage = function(...) stop("generation failed"))
  expect_error(generate(dir), "generation failed")
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
  local_mocked_bindings(generateNoiseImage = rcicr::generateNoiseImage)
})

test_that("a save that fails part-way leaves no .Rdata behind, so the next call can run", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(stimulusTime = function() minute(0))
  local_mocked_bindings(saveStimulusFile = function(names, file, envir) {
    writeLines("half written", file)
    stop("disk full")
  })
  expect_error(generate(dir), "disk full")
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)

  local_mocked_bindings(saveStimulusFile = function(names, file, envir) {
    save(list = names, file = file, envir = envir)
  })
  expect_no_error(generate(dir))
  expect_length(list.files(dir, pattern = "\\.Rdata$"), 1)
})

interrupt <- function() signalCondition(structure(class = c("interrupt", "condition"), list()))

abort <- function(expr) tryCatch(expr, interrupt = function(e) "aborted")

test_that("an abort removes everything the call created, serially", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(generateNoiseImage = function(...) interrupt())
  expect_identical(abort(generate(dir)), "aborted")
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
})

test_that("a parallel abort leaves the folder empty", {
  dir <- withr::local_tempdir()
  local_mocked_bindings(progressOption = function(pb, cl) list(progress = function(n) interrupt()))
  expect_identical(abort(suppressWarnings(utils::capture.output(
    generateStimuli2IFC(list(face = base_png()), n_trials = 40, img_size = 16, stimulus_path = dir,
                        seed = 1, ncores = 2, nscales = 1)
  ))), "aborted")
  Sys.sleep(1)
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
})

test_that("the cleanup kills a worker still writing, before removing its files", {
  # The parallel abort above interrupts as a result arrives, when the other
  # worker is finishing too, so it cannot catch a write that lands after the
  # cleanup. This drives the cleanup directly against a process that writes a
  # reserved path a second later.
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  reserved <- file.path(dir, sprintf("rcic_face_1_%05d_ori.png", 1:2))
  file.create(reserved)
  # One writer per worker, as with ncores = 2.
  writers <- lapply(reserved, function(path) {
    parallel::mcparallel({
      Sys.sleep(1)
      writeLines("late", path)
    })
  })
  pids <- vapply(writers, `[[`, integer(1), "pid")
  rcicr:::releaseStimulusCall(FALSE, NULL, pids, reserved, NULL, character())
  Sys.sleep(2)
  expect_false(any(file.exists(reserved)))
  parallel::mccollect(writers, wait = FALSE)
})

test_that("an abort during the reservation itself leaves no placeholder", {
  dir <- withr::local_tempdir()
  calls <- 0
  local_mocked_bindings(pathTaken = function(path) {
    calls <<- calls + 1
    # The first 12 calls are the up-front check; interrupt partway through the loop.
    if (calls == 12 + 5) interrupt()
    file.exists(path)
  })
  expect_identical(abort(generate(dir)), "aborted")
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
})

test_that("empty leftovers of a killed call are named as such", {
  dir <- withr::local_tempdir()
  file.create(rcicr:::stimulusPngPaths(dir, "rcic", "face", 1, 6))
  expect_error(generate(dir), "as placeholders left by a generateStimuli2IFC\\(\\) call that was killed")
  writeLines("real", rcicr:::stimulusPngPath(dir, "rcic", "face", 1, 1, "ori"))
  expect_error(generate(dir), "12 of the 12 PNG files")
  expect_no_match(tryCatch(generate(dir), error = conditionMessage), "placeholders")
})

test_that("a save that fails after the workers have exited signals no PID", {
  dir <- withr::local_tempdir()
  released <- NULL
  original <- rcicr:::releaseStimulusCall
  local_mocked_bindings(
    saveStimulusFile = function(names, file, envir) stop("disk full"),
    releaseStimulusCall = function(finished, cl, worker_pids, ...) {
      released <<- list(finished = finished, worker_pids = worker_pids)
      original(finished, cl, worker_pids, ...)
    }
  )
  expect_error(suppressWarnings(utils::capture.output(
    generateStimuli2IFC(list(face = base_png()), n_trials = 4, img_size = 16, stimulus_path = dir,
                        seed = 1, ncores = 2, nscales = 1)
  )), "disk full")
  expect_false(released$finished)
  expect_null(released$worker_pids)
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)
})

test_that("a dangling symlink at the .Rdata name is taken, and a failure never removes it", {
  skip_on_os("windows")
  dir <- withr::local_tempdir()
  local_mocked_bindings(stimulusTime = function() minute(0))
  outside <- file.path(withr::local_tempdir(), "target.Rdata")
  link <- rcicr:::stimulusRdataPath(dir, "rcic", 1, minute(0))
  file.symlink(outside, link)
  expect_error(generate(dir), "already exists, and a stimulus file is never overwritten")
  expect_identical(Sys.readlink(link), outside)
  expect_false(file.exists(outside))
})
