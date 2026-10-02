# Every classification image records the scaling its $scaled was made with
# (#9), and the constant recorded is the one that rendered it.

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr))
  out
}

ci_of <- function(rdata, ..., stimuli = 1:8, responses = c(1, -1, -1, 1, 1, 1, -1, 1)) {
  quietly(generateCI(stimuli, responses, "base", rdata, save_as_png = FALSE, n_cores = 1, ...))
}

test_that("each method records itself and a constant that reproduces $scaled", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  independent <- ci_of(rdata, scaling = "independent")
  k <- attr(independent, "scaling")$constant
  expect_identical(attr(independent, "scaling")$method, "independent")
  expect_identical(k, max(abs(range(independent$ci))))
  expect_equal(independent$scaled, (independent$ci + k) / (2 * k))
  expect_false(isTRUE(all.equal(independent$scaled, (independent$ci + 2 * k) / (4 * k))))

  constant <- ci_of(rdata, scaling = "constant", scaling_constant = 0.4)
  expect_identical(attr(constant, "scaling"), list(method = "constant", constant = 0.4))
  expect_equal(constant$scaled, (constant$ci + 0.4) / 0.8)

  for (method in c("none", "matched")) {
    expect_identical(attr(ci_of(rdata, scaling = method), "scaling"),
                     list(method = method, constant = NA_real_))
  }
  expect_identical(attr(ci_of(rdata, scaling = "squash"), "scaling"),
                   list(method = "none", constant = NA_real_))
})

test_that("an all-zero CI under 'independent' records 0", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  zero <- ci_of(rdata, stimuli = c(1, 1), responses = c(1, -1), scaling = "independent")
  expect_true(all(zero$ci == 0))
  expect_identical(attr(zero, "scaling")$constant, 0)
})

test_that("a masked CI records the constant of its unmasked pixels", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  open_ci <- ci_of(rdata)
  peak <- which.max(abs(open_ci$ci))
  mask <- matrix(1, 32, 32)
  mask[peak] <- 0
  masked <- ci_of(rdata, mask = mask)
  k <- attr(masked, "scaling")$constant
  expect_lt(k, attr(open_ci, "scaling")$constant)
  expect_equal(masked$scaled[-peak], (masked$ci[-peak] + k) / (2 * k))
})

test_that("participants record one constant per individual PNG, by ID, masked as the PNG is", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  dir <- withr::local_tempdir()
  # Out of order, and numeric IDs that sort differently as text.
  participants <- c(10, 10, 9, 9, 2, 2, 10, 9)
  open_ci <- ci_of(rdata, participants = participants, save_individual_cis = TRUE, targetpath = dir)
  individual <- attr(open_ci, "scaling")$individual
  expect_identical(individual$method, "independent")
  expect_named(individual$constant, c("2", "9", "10"))

  e <- new.env()
  load(rdata, envir = e)
  base <- e$base_faces$base
  for (id in names(individual$constant)) {
    rows <- participants == as.numeric(id)
    own <- ci_of(rdata, stimuli = (1:8)[rows], responses = c(1, -1, -1, 1, 1, 1, -1, 1)[rows],
                 scaling = "independent")
    k <- individual$constant[[id]]
    expect_identical(k, attr(own, "scaling")$constant, info = id)
    png_pixels <- png::readPNG(file.path(dir, "individual_cis", paste0("ci_", id, ".png")))
    expect_equal(png_pixels, ((own$ci + k) / (2 * k) + base) / 2, tolerance = 1 / 255, info = id)
  }

  # A mask over each participant's largest absolute value lowers every constant.
  peaks <- vapply(names(individual$constant), function(id) {
    rows <- participants == as.numeric(id)
    own <- ci_of(rdata, stimuli = (1:8)[rows], responses = c(1, -1, -1, 1, 1, 1, -1, 1)[rows])
    which.max(abs(own$ci))
  }, integer(1))
  mask <- matrix(1, 32, 32)
  mask[peaks] <- 0
  masked <- attr(ci_of(rdata, participants = participants, mask = mask), "scaling")$individual
  expect_true(all(masked$constant < individual$constant))

  fixed <- attr(ci_of(rdata, participants = participants, individual_scaling = "constant",
                      individual_scaling_constant = 0.2), "scaling")$individual
  expect_identical(fixed, list(method = "constant", constant = 0.2))
  expect_null(attr(ci_of(rdata), "scaling")$individual)
})

test_that("autoscale() records its constant for $scaled and keeps the record behind $combined", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  cis <- list(a = ci_of(rdata, scaling = "independent"),
              b = ci_of(rdata, participants = c(1, 1, 1, 1, 2, 2, 2, 2), scaling = "constant",
                        scaling_constant = 0.3))
  before <- lapply(cis, attr, "scaling")
  scaled <- quietly(autoscale(cis, save_as_pngs = FALSE))
  k <- attr(scaled$a, "scaling")$constant
  expect_identical(attr(scaled$b, "scaling")$constant, k)
  for (name in names(cis)) {
    record <- attr(scaled[[name]], "scaling")
    expect_identical(record$method, "autoscale")
    expect_equal(scaled[[name]]$scaled, (scaled[[name]]$ci + k) / (2 * k))
    expect_identical(record$combined, before[[name]][c("method", "constant")])
    expect_identical(attr(scaled[[name]], "trial_design"), attr(cis[[name]], "trial_design"))
  }
  expect_identical(attr(scaled$b, "scaling")$individual, before$b$individual)
  # $combined is the image the earlier record describes.
  kb <- before$b$constant
  expect_equal(scaled$b$combined, ((scaled$b$ci + kb) / (2 * kb) + scaled$b$base) / 2)

  twice <- quietly(autoscale(scaled, save_as_pngs = FALSE))
  expect_identical(attr(twice$a, "scaling")$combined, before$a[c("method", "constant")])

  bare <- quietly(autoscale(list(x = list(ci = cis$a$ci, base = cis$a$base)), save_as_pngs = FALSE))
  expect_true("combined" %in% names(attr(bare$x, "scaling")))
  expect_null(attr(bare$x, "scaling")$combined)
})

test_that("the batch functions return the autoscale record, and 'none' behind $combined", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  d <- data.frame(g = rep(c("x", "y"), each = 4), s = 1:8, r = c(1, -1, -1, 1, 1, 1, -1, 1))
  for (f in list(batchGenerateCI, batchGenerateCI2IFC)) {
    cis <- quietly(f(d, "g", "s", "r", "base", rdata, save_as_png = FALSE))
    k <- attr(cis[[1]], "scaling")$constant
    for (ci in cis) {
      expect_identical(attr(ci, "scaling")$method, "autoscale")
      expect_identical(attr(ci, "scaling")$constant, k)
      expect_identical(attr(ci, "scaling")$combined, list(method = "none", constant = NA_real_))
      expect_equal(ci$combined, (ci$ci + ci$base) / 2)
    }
  }
})

test_that("the record changes no pixel", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  ci <- ci_of(rdata)
  e <- new.env()
  load(rdata, envir = e)
  params <- rcicr:::selectStimulusParams(e$stimuli_params, "base", 1:8)
  noise <- generateCINoise(params, c(1, -1, -1, 1, 1, 1, -1, 1), e$p)
  expect_identical(ci$ci, noise)
  expect_identical(names(ci), c("ci", "scaled", "base", "combined"))
})
