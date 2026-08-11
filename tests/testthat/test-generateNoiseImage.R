test_that("generateNoiseImage returns a matrix of the expected size", {
  p <- generateNoisePattern(16, nscales = 1)
  img <- generateNoiseImage(rep(1, max(p$patchIdx)), p)
  expect_equal(dim(img), c(16, 16))
})

test_that("mismatched params length errors", {
  p <- generateNoisePattern(16, nscales = 1)
  expect_error(
    generateNoiseImage(rep(1, max(p$patchIdx) - 5), p),
    "Stimulus generation aborted"
  )
})

test_that("legacy 0-indexed patchIdx warns but still returns a valid image", {
  p0 <- generateNoisePattern(16, nscales = 1, pre_0.3.0 = TRUE)
  expect_equal(min(p0$patchIdx), 0)

  params <- rep(1, max(p0$patchIdx) + 1)
  expect_warning(
    result <- generateNoiseImage(params, p0),
    "patch indices start at 0"
  )
  expect_equal(dim(result), c(16, 16))
  expect_false(anyNA(result))
})

# The 0-based path drops every patchIdx == 0 cell (R drops a 0 subscript) and
# recycles a too-short vector over the patch array, which reads like a
# whole-image misalignment but is masked because those cells sit on an all-zero
# patch layer. This checks the masking implication and pins the output against an
# honest "one patch not shown" oracle. Why pre_0.3.0 = TRUE reconstructs the
# genuine legacy layout, and why the effect is one patch rather than the whole
# image, are in DECISIONS.md, "4096 -> 4092 parameters"; see also issue #221.
test_that("0-indexed warning branch drops exactly one patch, not the whole image", {
  for (ns in 1:2) {
    p0 <- generateNoisePattern(16, nscales = ns, pre_0.3.0 = TRUE)
    # masking implication: patchIdx == 0 => patches == 0 (converse is false)
    expect_true(all(p0$patches[p0$patchIdx == 0] == 0))

    set.seed(42)
    params <- rnorm(max(p0$patchIdx) + 1)
    current <- suppressWarnings(generateNoiseImage(params, p0))

    # Oracle: patch idx 0 contributes nothing, every other patch k uses params[k].
    idx <- p0$patchIdx
    idx[idx == 0] <- NA
    w <- params[idx]
    w[is.na(w)] <- 0
    oracle <- rowMeans(p0$patches * array(w, dim(p0$patches)), dims = 2)

    # expect_identical, not expect_equal: the two are bit-identical (the zero
    # cells are contiguous at the end, so the non-zero positions fill exactly as
    # the oracle does and the same rowMeans runs on the same bits), and the
    # NEWS/DECISIONS claim is "exact" - a tolerant compare would pass on drift.
    expect_identical(current, oracle)
  }
})

# NOTE: support for the pre-0.3.3 sinusoids/sinIdx list shape is currently
# broken; test-known-bugs.R holds a failing test for the intended behaviour.

test_that("all-zero params gives an all-zero image", {
  p <- generateNoisePattern(16, nscales = 1)
  result <- generateNoiseImage(rep(0, max(p$patchIdx)), p)
  expect_equal(result, matrix(0, 16, 16))
})

# The patch averaging is vectorised with rowMeans(..., dims = 2). rowMeans() on
# a 3-D array defaults to dims = 1, which collapses dimensions 2 and 3 and
# returns a vector that gets silently recycled back to img_size x img_size. The
# two tests below pin the intended behaviour against an independent oracle so
# that mistake cannot creep back in. See PR #122.

test_that("noise is the per-pixel mean across patch layers (independent oracle)", {
  # A NON-square patch stack: a square one would hide a transposition.
  set.seed(42)
  patches <- array(rnorm(4 * 6 * 3), dim = c(4, 6, 3))
  p <- list(patches = patches, patchIdx = array(1:3, dim = c(4, 6, 3)),
            noise_type = "sinusoid")
  params <- c(0.5, -2, 1.5)

  # The definition, written out by hand: neither apply() nor rowMeans().
  oracle <- matrix(NA_real_, 4, 6)
  for (i in 1:4) {
    for (j in 1:6) {
      total <- 0
      for (k in 1:3) total <- total + patches[i, j, k] * params[p$patchIdx[i, j, k]]
      oracle[i, j] <- total / 3
    }
  }

  result <- generateNoiseImage(params, p)
  expect_equal(dim(result), c(4, 6))
  expect_equal(result, oracle)
})

test_that("vectorised averaging matches the original apply() implementation", {
  # Before 1.1.0 the average was computed as
  #   apply(p$patches * array(params[p$patchIdx], dim(p$patches)), 1:2, mean)
  # This pins the current implementation to that one across noise types, spatial
  # scales and seeds. The two are not bit-identical -- they sum in a different
  # order, so they differ by about 1 ULP -- hence a tolerance rather than
  # identical(). Anything larger than that means the averaging changed.
  configs <- expand.grid(
    noise_type = c("sinusoid", "gabor"),
    nscales = 1:3,
    seed = c(1, 2),
    stringsAsFactors = FALSE
  )

  worst <- 0
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, ]
    p <- generateNoisePattern(16, nscales = cfg$nscales,
                              noise_type = cfg$noise_type, sigma = 5)
    set.seed(cfg$seed)
    params <- rnorm(max(p$patchIdx))

    original <- apply(p$patches * array(params[p$patchIdx], dim(p$patches)),
                      1:2, mean)
    current <- generateNoiseImage(params, p)

    label <- paste0(cfg$noise_type, ", nscales=", cfg$nscales, ", seed=", cfg$seed)
    expect_equal(current, original, tolerance = 1e-12, info = label)
    worst <- max(worst, max(abs(current - original)))
  }

  # Absolute agreement, independent of the relative tolerance above. Patch
  # values are O(1), so 1e-14 is still several orders of magnitude of headroom
  # over the ~1e-17 that differing summation order actually produces.
  expect_lt(worst, 1e-14)
})

test_that("noise columns are not recycled copies of a single column", {
  # The dims = 1 bug produces an image whose every column is identical.
  set.seed(1)
  p <- generateNoisePattern(16, nscales = 1)
  result <- generateNoiseImage(rnorm(max(p$patchIdx)), p)
  expect_false(all(result[, 1] == result[, 2]))
})
