# InfoVal of a masked classification image: the CI's norm over its unmasked
# pixels, against a reference over the same pixels (#374).

quietly <- function(expr) {
  out <- NULL
  utils::capture.output(out <- suppressMessages(suppressWarnings(expr)))
  out
}

responses8 <- c(1, -1, -1, 1, 1, 1, -1, 1)

top_mask <- function(rows = 1:5, size = 32) {
  mask <- matrix(1, size, size)
  mask[rows, ] <- 0
  mask
}

masked_fixture <- function() {
  rdata <- make_fixture_rdata(withr::local_tempdir(.local_envir = parent.frame()), n_trials = 8,
                              nscales = 2)
  list(rdata = rdata,
       ci = quietly(generateCI(1:8, responses8, "base", rdata, save_as_png = FALSE,
                               mask = top_mask())))
}

# Every stimulus rendered, masked rows dropped, times the same seeded responses.
brute_norms <- function(rdata, kept, stimuli, iter, response_seed) {
  e <- new.env()
  load(rdata, envir = e)
  noise <- vapply(stimuli, function(s) as.vector(generateNoiseImage(e$stimuli_params$base[s, ], e$p)),
                  numeric(32 * 32))[kept, , drop = FALSE]
  set.seed(response_seed)
  vapply(seq_len(iter), function(i) {
    r <- ((runif(length(stimuli)) > 0.5) * 2) - 1
    norm(noise %*% r / length(stimuli), "f")
  }, numeric(1))
}

test_that("a masked reference equals one computed by brute force, under both methods", {
  f <- masked_fixture()
  kept <- as.vector(top_mask() == 1)
  for (method in c("gram", "images")) {
    for (stimuli in list(NULL, c(2, 3, 5, 7, 8))) {
      norms <- quietly(generateReferenceDistribution2IFC(
        f$rdata, iter = 30, ncores = 1, response_seed = 5, save_rdata = FALSE, mask = top_mask(),
        reference_stimuli = stimuli, reference_method = method
      ))
      used <- if (is.null(stimuli)) 1:8 else stimuli
      expect_equal(norms, brute_norms(f$rdata, kept, used, 30, 5), tolerance = 1e-12,
                   info = paste(method, length(used)))
      unmasked <- quietly(generateReferenceDistribution2IFC(
        f$rdata, iter = 30, ncores = 1, response_seed = 5, save_rdata = FALSE,
        reference_stimuli = stimuli, reference_method = method
      ))
      expect_false(isTRUE(all.equal(norms, unmasked)))
    }
  }
})

test_that("a masked CI gets a finite InfoVal from its own pixels and reference", {
  f <- masked_fixture()
  value <- quietly(computeInfoVal2IFC(f$ci, f$rdata, iter = 40, response_seed = 3))
  expect_true(is.finite(value))
  norms <- quietly(generateReferenceDistribution2IFC(f$rdata, iter = 40, ncores = 1,
                                                     response_seed = 3, save_rdata = FALSE,
                                                     mask = top_mask()))
  ci_norm <- norm(matrix(f$ci$ci[!is.na(f$ci$ci)]), "f")
  expect_identical(value, (ci_norm - median(norms)) / mad(norms))
})

test_that("each mask is stored once and reused, apart from the default reference", {
  f <- masked_fixture()
  other <- quietly(generateCI(1:8, responses8, "base", f$rdata, save_as_png = FALSE,
                              mask = top_mask(28:32)))
  plain <- quietly(generateCI(1:8, responses8, "base", f$rdata, save_as_png = FALSE))
  first <- quietly(computeInfoVal2IFC(f$ci, f$rdata, iter = 20))
  quietly(computeInfoVal2IFC(other, f$rdata, iter = 20))
  quietly(computeInfoVal2IFC(plain, f$rdata, iter = 20))
  e <- new.env()
  load(f$rdata, envir = e)
  expect_length(e$reference_norms_by_stimuli, 2)
  expect_false(is.null(e$reference_norms))
  expect_null(e$reference_norms_by_stimuli[[1]]$reference_stimuli)
  expect_false(identical(e$reference_norms_by_stimuli[[1]]$mask, e$reference_norms_by_stimuli[[2]]$mask))

  local_mocked_bindings(referenceNorms = function(...) stop("simulated again"))
  expect_identical(quietly(computeInfoVal2IFC(f$ci, f$rdata, iter = 20)), first)
})

test_that("a reference stored with generateReferenceDistribution2IFC(mask =) is the one reused", {
  f <- masked_fixture()
  stored <- quietly(generateReferenceDistribution2IFC(f$rdata, iter = 25, ncores = 1,
                                                      response_seed = 9, mask = top_mask()))
  local_mocked_bindings(referenceNorms = function(...) stop("simulated again"))
  ci_norm <- norm(matrix(f$ci$ci[!is.na(f$ci$ci)]), "f")
  expect_identical(quietly(computeInfoVal2IFC(f$ci, f$rdata)),
                   (ci_norm - median(stored)) / mad(stored))
})

test_that("a mask object the file already holds is saved back unchanged", {
  f <- masked_fixture()
  e <- new.env()
  load(f$rdata, envir = e)
  e$mask <- "the file's own"
  save(list = ls(e), file = f$rdata, envir = e)
  quietly(generateReferenceDistribution2IFC(f$rdata, iter = 5, ncores = 1, mask = top_mask()))
  quietly(generateReferenceDistribution2IFC(f$rdata, iter = 5, ncores = 1))
  e <- new.env()
  load(f$rdata, envir = e)
  expect_identical(e$mask, "the file's own")
})

test_that("the batch scores masked and unmasked CIs as the single function does", {
  f <- masked_fixture()
  plain <- quietly(generateCI(1:8, responses8, "base", f$rdata, save_as_png = FALSE))
  cis <- list(masked = f$ci, plain = plain, again = f$ci)
  calls <- 0
  real <- rcicr:::referenceNorms
  local_mocked_bindings(referenceNorms = function(...) {
    calls <<- calls + 1
    real(...)
  })
  batch <- quietly(batchComputeInfoVal2IFC(cis, f$rdata, iter = 20, response_seed = 4))
  expect_identical(calls, 2)
  single <- vapply(cis, function(ci) quietly(computeInfoVal2IFC(ci, f$rdata, iter = 20,
                                                                response_seed = 4)), numeric(1))
  expect_identical(batch, single)
  expect_identical(batch[["masked"]], batch[["again"]])
})

test_that("a mask that covers everything, or the wrong size, stops", {
  f <- masked_fixture()
  expect_error(quietly(generateReferenceDistribution2IFC(f$rdata, iter = 5, ncores = 1,
                                                         mask = matrix(0, 32, 32))),
               "Every pixel is masked")
  all_na <- f$ci
  all_na$ci[] <- NA
  expect_error(quietly(computeInfoVal2IFC(all_na, f$rdata, iter = 5)), "Every pixel is masked")
  expect_error(quietly(generateReferenceDistribution2IFC(f$rdata, iter = 5, ncores = 1,
                                                         mask = matrix(1, 16, 16))),
               "not of the same dimensions")
})
