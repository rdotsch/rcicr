# batchGenerateCI() and batchGenerateCI2IFC() share one loop, and can average
# participants nested in each unit (#87).

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressWarnings(expr))
  out
}

batch_data <- function() {
  data.frame(
    cond = rep(c("a", "b", NA), c(8, 6, 2)),
    pid = c(1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5),
    stim = c(1:8, 1:6, 7:8),
    resp = c(1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1)
  )
}

# What both bodies did before sharing one: generateCI() per unit, named
# <base>[_<label>]_<by>_<unit>, then autoscale() over the list.
per_unit_oracle <- function(data, rdata, scaling = "autoscale", constant = 0.1, label = "",
                            antiCI = FALSE) {
  data <- as.data.frame(data)
  data <- data[!is.na(data$cond), ]
  inner <- if (scaling == "autoscale") "none" else scaling
  cis <- list()
  for (unit in unique(data$cond)) {
    rows <- data[data$cond == unit, ]
    name <- paste0("base_", if (nzchar(label)) paste0(label, "_") else "", "cond_", unit)
    cis[[name]] <- quietly(generateCI(rows$stim, rows$resp, "base", rdata, save_as_png = FALSE,
                                      antiCI = antiCI, scaling = inner, scaling_constant = constant))
  }
  if (scaling == "autoscale") cis <- quietly(autoscale(cis, save_as_pngs = FALSE))
  cis
}

pixels <- function(cis) lapply(cis, function(ci) unclass(ci)[c("ci", "scaled", "base", "combined")])

test_that("both functions give what a generateCI() call per unit gives", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  cases <- list(
    list(),
    list(scaling = "independent"),
    list(scaling = "constant", constant = 0.3),
    list(label = "run1"),
    list(antiCI = TRUE)
  )
  for (args in cases) {
    expected <- do.call(per_unit_oracle, c(list(batch_data(), rdata), args))
    common <- list(data = batch_data(), by = "cond", stimuli = "stim", responses = "resp",
                   baseimage = "base", rdata = rdata, save_as_png = FALSE)
    for (f in list(batchGenerateCI, batchGenerateCI2IFC)) {
      actual <- quietly(do.call(f, c(common, args)))
      expect_identical(names(actual), names(expected))
      expect_identical(pixels(actual), pixels(expected))
    }
  }
  expect_identical(pixels(quietly(batchGenerateCI(tibble::as_tibble(batch_data()), "cond", "stim",
                                                  "resp", "base", rdata, save_as_png = FALSE))),
                   pixels(per_unit_oracle(batch_data(), rdata)))
})

test_that("positional calls bind every argument as they did", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  dir <- withr::local_tempdir()
  d <- batch_data()
  named <- quietly(batchGenerateCI(d, "cond", "stim", "resp", "base", rdata, save_as_png = TRUE,
                                   targetpath = dir, label = "L", antiCI = TRUE,
                                   scaling = "constant", constant = 0.3))
  positional <- quietly(batchGenerateCI(d, "cond", "stim", "resp", "base", rdata, TRUE, dir, "L",
                                        TRUE, "constant", 0.3))
  expect_identical(pixels(positional), pixels(named))
  expect_named(positional, c("base_L_cond_a", "base_L_cond_b"))

  named2 <- quietly(batchGenerateCI2IFC(d, "cond", "stim", "resp", "base", rdata, save_as_png = TRUE,
                                        targetpath = dir, antiCI = TRUE, scaling = "constant",
                                        constant = 0.3, label = "L"))
  positional2 <- quietly(batchGenerateCI2IFC(d, "cond", "stim", "resp", "base", rdata, TRUE, dir,
                                             TRUE, "constant", 0.3, "L"))
  expect_identical(pixels(positional2), pixels(named2))
  expect_identical(pixels(positional2), pixels(named))
})

test_that("participants nested in a unit are averaged, not pooled", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  d <- batch_data()
  d$pid[d$cond %in% "b"] <- c(1, 1, 1, 2, 2, 2) # IDs reused across conditions
  nested <- quietly(batchGenerateCI(d, "cond", "stim", "resp", "base", rdata, save_as_png = FALSE,
                                    scaling = "none", participants = "pid"))
  pooled <- quietly(batchGenerateCI(d, "cond", "stim", "resp", "base", rdata, save_as_png = FALSE,
                                    scaling = "none"))
  for (unit in c("a", "b")) {
    rows <- d[d$cond %in% unit, ]
    expected <- quietly(generateCI(rows$stim, rows$resp, "base", rdata, participants = rows$pid,
                                   save_as_png = FALSE, scaling = "none", n_cores = 1))
    expect_identical(nested[[paste0("base_cond_", unit)]]$ci, expected$ci)
    expect_identical(attr(nested[[paste0("base_cond_", unit)]], "trial_design")$n_participants, 2L)
  }
  # Condition a: participant 1 gave 6 trials and participant 2 gave 2, so
  # averaging their CIs weights the trials differently from pooling them.
  expect_false(isTRUE(all.equal(nested$base_cond_a$ci, pooled$base_cond_a$ci)))
  expect_identical(
    pixels(quietly(batchGenerateCI2IFC(d, "cond", "stim", "resp", "base", rdata, save_as_png = FALSE,
                                       scaling = "none", participants = "pid"))),
    pixels(nested)
  )
})

test_that("a missing participants column or missing IDs stop before any CI", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  local_mocked_bindings(generateCI = function(...) stop("a CI was computed"))
  call <- function(d, participants = "pid") {
    batchGenerateCI(d, "cond", "stim", "resp", "base", rdata, save_as_png = FALSE,
                    participants = participants)
  }
  expect_error(call(batch_data(), "subject"), "participants must name one column of data; 'subject'")
  partly <- batch_data()
  partly$pid[2] <- NA
  expect_error(call(partly), "1 rows have no participant ID in column pid \\(cond a\\)")
  all_of_b <- batch_data()
  all_of_b$pid[all_of_b$cond %in% "b"] <- NA
  expect_error(call(all_of_b), "6 rows have no participant ID in column pid \\(cond b\\)")
  # Rows of the dropped NA group do not count.
  ungrouped <- batch_data()
  ungrouped$pid[is.na(ungrouped$cond)] <- NA
  expect_error(call(ungrouped), "a CI was computed")
})

test_that("nested participants under the default 'autoscale' scale the averaged CIs together", {
  rdata <- make_fixture_rdata(withr::local_tempdir(), n_trials = 8)
  d <- batch_data()
  d$pid[d$cond %in% "b"] <- c(1, 1, 1, 2, 2, 2)
  nested <- quietly(batchGenerateCI(d, "cond", "stim", "resp", "base", rdata, save_as_png = FALSE,
                                    participants = "pid"))
  expected <- list()
  for (unit in c("a", "b")) {
    rows <- d[d$cond %in% unit, ]
    expected[[paste0("base_cond_", unit)]] <- quietly(generateCI(
      rows$stim, rows$resp, "base", rdata, participants = rows$pid, save_as_png = FALSE,
      scaling = "none", n_cores = 1
    ))
  }
  expected <- quietly(autoscale(expected, save_as_pngs = FALSE))
  expect_identical(pixels(nested), pixels(expected))
  expect_identical(attr(nested$base_cond_a, "scaling")$method, "autoscale")
})
