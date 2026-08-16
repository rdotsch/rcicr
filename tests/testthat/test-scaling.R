# Coverage for the four scaling methods of generateCI().
#
# Scaling is the most consequential user-facing choice in this package -- it
# decides what a classification image actually looks like, and generateCI()'s
# roxygen header spends most of its length on it. Each method is tested against
# the mathematical property that defines it rather than against a stored number,
# so these tests say what each method is *for* and stay meaningful if the
# implementation is rewritten.
#
# The single most important property is the last one: scaling must never touch
# `ci`. Published results are computed from `ci`; `scaled` exists for display.

make_ci_fixture <- function(env = parent.frame()) {
  tmp <- withr::local_tempdir(.local_envir = env)
  rdata <- make_fixture_rdata(tmp, img_size = 32, n_trials = 12, nscales = 1, seed = 1) # nolint: object_usage_linter.
  list(tmp = tmp, rdata = rdata, responses = rep(c(1, -1), 6))
}

ci_with_scaling <- function(f, ...) {
  generateCI(
    stimuli = seq_along(f$responses), responses = f$responses,
    baseimage = "base", rdata = f$rdata, save_as_png = FALSE, ...
  )
}

test_that("scaling = 'none' leaves the CI untouched", {
  f <- make_ci_fixture()
  res <- ci_with_scaling(f, scaling = "none")

  expect_identical(res$scaled, res$ci)
})

test_that("scaling = 'constant' applies (ci + k) / 2k and centres zero at 0.5", {
  f <- make_ci_fixture()
  k <- 0.5
  res <- ci_with_scaling(f, scaling = "constant", scaling_constant = k)

  expect_equal(res$scaled, (res$ci + k) / (2 * k))

  # The defining property: a pixel the participants pushed neither way lands at
  # mid grey, so the base image is unchanged there.
  expect_equal(min(abs(res$scaled - 0.5)), min(abs(res$ci)) / (2 * k))
})

test_that("scaling = 'constant' warns when the constant is too small to avoid clipping", {
  f <- make_ci_fixture()

  # A constant well below the CI's own range forces values outside [0, 1], which
  # png::writePNG() would silently clip. The user has to be told.
  expect_warning(
    ci_with_scaling(f, scaling = "constant", scaling_constant = 1e-4),
    "exceed possible intensity range"
  )
})

test_that("scaling = 'matched' matches the intensity range of the base image", {
  f <- make_ci_fixture()
  res <- ci_with_scaling(f, scaling = "matched")

  expect_equal(range(res$scaled), range(res$base))
})

test_that("scaling = 'independent' uses the whole range without clipping", {
  f <- make_ci_fixture()
  res <- ci_with_scaling(f, scaling = "independent")

  # Nothing clips...
  expect_gte(min(res$scaled), 0)
  expect_lte(max(res$scaled), 1)

  # ...and the constant chosen is the smallest one that achieves that, so the
  # extreme of whichever tail is larger in absolute value lands exactly on the
  # boundary. Only one end is touched: the method is symmetric about 0.5, not
  # a range stretch.
  touches_boundary <- isTRUE(all.equal(max(res$scaled), 1)) ||
    isTRUE(all.equal(min(res$scaled), 0))
  expect_true(touches_boundary)

  # Equivalently: it is 'constant' scaling with k = max(abs(ci)).
  k <- max(abs(res$ci))
  expect_equal(res$scaled, (res$ci + k) / (2 * k))
})

test_that("an unrecognised scaling method warns and falls back to none", {
  f <- make_ci_fixture()

  expect_warning(
    res <- ci_with_scaling(f, scaling = "nonesuch"),
    "not found"
  )
  expect_identical(res$scaled, res$ci)
})

test_that("scaling changes only $scaled, never $ci", {
  # This is the reproducibility guarantee. Researchers compute infoVal, z-maps
  # and correlations from $ci; $scaled exists so the CI can be written to a PNG.
  # If a scaling method ever leaked into $ci, published numbers would depend on
  # a display choice.
  f <- make_ci_fixture()

  methods <- list(
    none = list(scaling = "none"),
    constant = list(scaling = "constant", scaling_constant = 0.5),
    matched = list(scaling = "matched"),
    independent = list(scaling = "independent")
  )

  cis <- lapply(methods, function(args) do.call(ci_with_scaling, c(list(f), args))$ci)

  for (method in names(cis)) {
    expect_identical(cis[[method]], cis$none, info = method)
  }
})
