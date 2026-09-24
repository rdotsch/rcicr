# A tibble must be grouped exactly like the data.frame holding the same rows,
# compared CI by CI: a mis-grouped tibble still returns CIs of the right shape
# (#336).

batch_both <- function(fn, df, rdata_path, by = "pid") {
  run <- function(data) {
    suppressWarnings(fn(
      data = data, by = by, stimuli = "stim", responses = "resp",
      baseimage = "base", rdata = rdata_path, save_as_png = FALSE
    ))
  }
  list(df = run(df), tbl = run(tibble::as_tibble(df)))
}

expect_same_cis <- function(res) {
  expect_named(res$tbl, names(res$df))
  expect_identical(lapply(res$tbl, `[[`, "ci"), lapply(res$df, `[[`, "ci"))
}

for (fn_name in c("batchGenerateCI", "batchGenerateCI2IFC")) {
  fn <- get(fn_name)

  test_that(paste(fn_name, "groups a tibble like a data.frame"), {
    tmp <- withr::local_tempdir()
    rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 8, nscales = 1, seed = 1)

    for (k in 2:3) {
      set.seed(k)
      df <- data.frame(pid = rep(seq_len(k), each = 8), stim = rep(1:8, k),
                       resp = sample(c(1, -1), 8 * k, replace = TRUE))
      res <- batch_both(fn, df, rdata_path)
      expect_length(res$tbl, k)
      expect_same_cis(res)
    }
  })

  test_that(paste(fn_name, "drops rows without a group in a tibble too"), {
    tmp <- withr::local_tempdir()
    rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 8, nscales = 1, seed = 1)

    df <- data.frame(pid = c(rep("p1", 4), NA, rep("p2", 3)), stim = 1:8,
                     resp = c(1, -1, 1, 1, -1, -1, 1, -1))
    res <- batch_both(fn, df, rdata_path)
    expect_named(res$df, c("base_pid_p1", "base_pid_p2"))
    expect_same_cis(res)
  })

  test_that(paste(fn_name, "treats a factor group column like a character one"), {
    tmp <- withr::local_tempdir()
    rdata_path <- make_fixture_rdata(tmp, img_size = 32, n_trials = 8, nscales = 1, seed = 1)

    set.seed(4)
    df <- data.frame(pid = rep(c("p1", "p2"), each = 8), stim = rep(1:8, 2),
                     resp = sample(c(1, -1), 16, replace = TRUE))
    chr <- batch_both(fn, df, rdata_path)
    df$pid <- factor(df$pid)
    fct <- batch_both(fn, df, rdata_path)

    expect_identical(lapply(fct$df, `[[`, "ci"), lapply(chr$df, `[[`, "ci"))
    expect_identical(lapply(fct$tbl, `[[`, "ci"), lapply(chr$df, `[[`, "ci"))
  })
}
