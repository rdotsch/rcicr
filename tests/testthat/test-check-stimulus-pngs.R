# checkStimulusPNGs(): a stimulus .Rdata file checked against stimulus PNGs.
# Archives are written with the original base images; candidates are
# regenerated over a grey stand-in, as a researcher with a lost file would.

quiet <- function(expr) {
  utils::capture.output(value <- suppressMessages(expr))
  value
}

write_archive <- function(bases, ..., seed = 7, nscales = 2, n_trials = 4, label = "rcic") {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  files <- lapply(seq_along(bases), function(i) {
    make_square_png(file.path(dir, paste0("base", i, ".png")), size = 32, seed = i) # nolint: object_usage_linter.
  })
  names(files) <- bases
  quiet(generateStimuli2IFC(files, n_trials = n_trials, img_size = 32, stimulus_path = dir,
                            seed = seed, nscales = nscales, ncores = 1, label = label, ...))
  dir
}

write_candidate <- function(bases, ..., seed = 7, nscales = 2, n_trials = 4) {
  dir <- withr::local_tempdir(.local_envir = parent.frame())
  grey <- file.path(dir, "grey.png")
  png::writePNG(matrix(0.5, 32, 32), grey)
  files <- stats::setNames(rep(list(grey), length(bases)), bases)
  quiet(generateStimuli2IFC(files, n_trials = n_trials, img_size = 32, stimulus_path = dir,
                            seed = seed, nscales = nscales, ncores = 1, save_as_png = FALSE,
                            maximize_baseimage_contrast = FALSE, ...))
  list.files(dir, pattern = "\\.Rdata$", full.names = TRUE)
}

rdata_in <- function(dir) list.files(dir, pattern = "\\.Rdata$", full.names = TRUE)

test_that("the file that wrote the PNGs matches every trial", {
  archive <- write_archive("face")
  result <- expect_silent(checkStimulusPNGs(rdata_in(archive), archive))
  expect_equal(result$share, rep(1, 4))
  expect_false(any(result$missing))
  expect_true(all(result$compared > 0))
  expect_identical(attr(result, "unchecked"), character(0))
})

test_that("a candidate with the right settings over a grey base matches every trial", {
  archive <- write_archive("face")
  expect_equal(checkStimulusPNGs(write_candidate("face"), archive)$share, rep(1, 4))
})

test_that("a wrong nscales matches only a small share", {
  archive <- write_archive("face")
  expect_true(all(checkStimulusPNGs(write_candidate("face", nscales = 3), archive)$share < 0.2))
})

test_that("a wrong seed is scored against the archive's seed, and named when it is not passed", {
  archive <- write_archive("face")
  candidate <- write_candidate("face", seed = 8)
  expect_true(all(checkStimulusPNGs(candidate, archive, seed = 7)$share < 0.2))
  expect_error(checkStimulusPNGs(candidate, archive), "label and seed must be the ones")
})

test_that("a 4096-wide pre-0.3.0 parameter matrix is truncated as generateCI() truncates it", {
  archive <- write_archive("face", nscales = 5)
  stored <- new.env()
  load(rdata_in(archive), envir = stored)
  params <- stored$stimuli_params$face
  expect_identical(ncol(params), 4092L)
  stored$stimuli_params$face <- cbind(params, matrix(0, nrow(params), 4))
  wide <- tempfile(fileext = ".Rdata")
  save(list = ls(stored), file = wide, envir = stored)
  expect_equal(checkStimulusPNGs(wide, archive)$share, rep(1, 4))
})

test_that("every base label is checked, so a wrong use_same_parameters shows after the first", {
  archive <- write_archive(c("first", "second"), use_same_parameters = FALSE)
  result <- checkStimulusPNGs(write_candidate(c("first", "second"), use_same_parameters = TRUE),
                              archive)
  expect_equal(result$share[result$base == "first"], rep(1, 4))
  expect_true(all(result$share[result$base == "second"] < 0.2))
})

test_that("PNG names are built as the generator builds them", {
  fractional <- write_archive("face", seed = 1.5)
  expect_false(any(checkStimulusPNGs(rdata_in(fractional), fractional)$missing))
  odd_label <- write_archive("face", label = "x_ori.png_y")
  expect_false(any(checkStimulusPNGs(rdata_in(odd_label), odd_label)$missing))
})

test_that("RGB copies of the PNGs give the same shares as grey ones", {
  archive <- write_archive("face")
  grey_shares <- checkStimulusPNGs(write_candidate("face", nscales = 3), archive)$share
  for (file in list.files(archive, pattern = "_(ori|inv)\\.png$", full.names = TRUE)) {
    img <- png::readPNG(file)
    png::writePNG(array(img, c(dim(img), 3)), file)
  }
  expect_length(dim(png::readPNG(list.files(archive, "_ori\\.png$", full.names = TRUE)[1])), 3L)
  expect_identical(checkStimulusPNGs(write_candidate("face", nscales = 3), archive)$share,
                   grey_shares)
})

test_that("a missing PNG is reported, and no PNG at all is an error", {
  archive <- write_archive("face")
  unlink(file.path(archive, "rcic_face_7_00002_inv.png"))
  expect_warning(result <- checkStimulusPNGs(rdata_in(archive), archive),
                 "1 trial\\(s\\) have no ori or inv PNG: face 2")
  expect_identical(result$missing, c(FALSE, TRUE, FALSE, FALSE))
  expect_true(is.na(result$share[2]))
  unlink(list.files(archive, pattern = "_inv\\.png$", full.names = TRUE))
  expect_warning(result <- checkStimulusPNGs(rdata_in(archive), archive),
                 "4 trial\\(s\\) have no ori or inv PNG")
  expect_true(all(result$missing))
  unlink(list.files(archive, pattern = "\\.png$", full.names = TRUE))
  expect_error(checkStimulusPNGs(rdata_in(archive), archive), "No stimulus PNG named for label")
})

test_that("archive PNGs a shorter candidate does not cover are reported", {
  # Shared parameters: a shorter draw is the first rows of the longer one, so
  # every row checked matches and only the unchecked PNGs show the shortfall.
  shared <- write_archive(c("first", "second"))
  expect_warning(fewer <- checkStimulusPNGs(write_candidate(c("first", "second"), n_trials = 3),
                                            shared),
                 "4 PNG\\(s\\) named for label")
  expect_equal(fewer$share, rep(1, 6))
  expect_setequal(attr(fewer, "unchecked"),
                  c("rcic_first_7_00004_ori.png", "rcic_first_7_00004_inv.png",
                    "rcic_second_7_00004_ori.png", "rcic_second_7_00004_inv.png"))
  own <- write_archive(c("first", "second"), use_same_parameters = FALSE)
  expect_warning(one_base <- checkStimulusPNGs(write_candidate("first",
                                                               use_same_parameters = FALSE),
                                               own),
                 "8 PNG\\(s\\) named for label")
  expect_equal(one_base$share, rep(1, 4))
  expect_true(all(grepl("^rcic_second_", attr(one_base, "unchecked"))))
})

test_that("archive PNGs with a hidden name are scanned too", {
  archive <- write_archive("face", label = ".rcic")
  expect_warning(result <- checkStimulusPNGs(write_candidate("face", n_trials = 3), archive,
                                             label = ".rcic"),
                 "2 PNG\\(s\\) named for label")
  expect_equal(result$share, rep(1, 3))
  expect_setequal(attr(result, "unchecked"),
                  c(".rcic_face_7_00004_ori.png", ".rcic_face_7_00004_inv.png"))
})

test_that("trial numbers of more than five digits are recognised", {
  archive <- write_archive("face")
  file.copy(file.path(archive, "rcic_face_7_00001_ori.png"),
            file.path(archive, "rcic_face_7_100000_ori.png"))
  expect_warning(result <- checkStimulusPNGs(rdata_in(archive), archive), "1 PNG\\(s\\)")
  expect_identical(attr(result, "unchecked"), "rcic_face_7_100000_ori.png")
})

test_that("png_dir is required and must exist", {
  archive <- write_archive("face")
  expect_error(checkStimulusPNGs(rdata_in(archive)), "png_dir must be")
  expect_error(checkStimulusPNGs(rdata_in(archive), file.path(archive, "nope")), "png_dir must be")
})
