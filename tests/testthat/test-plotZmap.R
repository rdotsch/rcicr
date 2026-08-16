test_that("plotZmap writes a PNG file", {
  tmp <- withr::local_tempdir()
  zmap <- matrix(5, 8, 8)

  plotZmap(
    zmap = zmap, sigma = 3, threshold = 3, decoration = FALSE,
    targetpath = tmp, filename = "zmap", size = 200
  )

  expect_true(file.exists(file.path(tmp, "zmap.png")))
})

test_that("decoration = FALSE works with a background image", {
  # Regression test. This branch used `if (bgimage != '')`, and bgimage is
  # normally an image matrix, so the condition had length img_size^2 -- an error
  # on R >= 4.2 rather than a silent first-element match. Every generateCI() call
  # passes a background image, so the whole undecorated path was dead. Same root
  # cause as the `mask` bug fixed in 1.1.0.
  tmp <- withr::local_tempdir()

  plotZmap(
    zmap = matrix(5, 8, 8), bgimage = matrix(0.5, 8, 8), sigma = 3,
    threshold = 3, decoration = FALSE, targetpath = tmp,
    filename = "withbg", size = 200
  )

  expect_true(file.exists(file.path(tmp, "withbg.png")))
})

test_that("decoration = FALSE works on a device too small for default margins", {
  # Regression test. plot.new() was called before par(mar = c(0,0,0,0)), and
  # plot.new() rejects a device that cannot fit the current margins. Since
  # generateCI() sizes the device to img_size, small stimulus sets could not
  # produce a z-map at all ("figure margins too large").
  tmp <- withr::local_tempdir()

  plotZmap(
    zmap = matrix(5, 8, 8), sigma = 3, threshold = 3, decoration = FALSE,
    targetpath = tmp, filename = "small", size = 32
  )

  expect_true(file.exists(file.path(tmp, "small.png")))
})

# The three tests above assert only that a file appeared, which was the right
# bar for two hard-error regressions but says nothing about what was drawn: a
# plotZmap() that wrote a blank image would pass all of them. These read the
# PNG back and check the thresholding actually took effect. Rendering onto a
# uniform background makes "nothing was painted" exactly equal to "the output
# is the background", which is a comparison rather than an eyeball.
render_zmap <- function(dir, zmap, name, threshold = 3, mask = NULL, bg = 0.5) {
  plotZmap(
    zmap = zmap, bgimage = matrix(bg, 8, 8), sigma = 3, threshold = threshold,
    mask = mask, decoration = FALSE, targetpath = dir, filename = name,
    size = 64
  )
  png::readPNG(file.path(dir, paste0(name, ".png")))
}

# Colour channels only, dropping any alpha plane.
#
# Whether a PNG device writes an alpha channel is a property of the graphics
# backend, not of what was drawn: cairo (Linux, Windows) writes RGB, macOS
# quartz writes RGBA. An opaque alpha plane is a solid block of 1s, so it adds a
# second distinct value to the array without a single pixel having been painted
# -- which is enough on its own to fail a "the whole image is one value"
# assertion, while every image-to-image comparison in this file stays green
# because both sides gain the same plane.
#
# That is exactly what happened when the suite first ran on macOS (R-hub,
# R-devel, 2026-07-28): 220 assertions passed and the single one below failed,
# reporting 2 distinct values where it wanted 1. It was the only assertion in
# the suite that counted distinct values rather than comparing two renders.
colour_channels <- function(img) {
  d <- dim(img)
  if (length(d) == 2) return(img)  # greyscale, no alpha plane
  if (d[3] %in% c(2, 4)) img[, , -d[3], drop = FALSE] else img
}

test_that("sub-threshold z-scores are not painted and supra-threshold ones are", {
  tmp <- withr::local_tempdir()

  below <- render_zmap(tmp, matrix(1, 8, 8), "below")
  above <- render_zmap(tmp, matrix(5, 8, 8), "above")

  # Every cell is under threshold, so the z-map contributes nothing and the
  # result is the bare background: a single grey value across the image. The
  # info= is deliberate -- a bare "expected length 1, got 2" says nothing about
  # which backend wrote what, which is the one thing worth knowing here.
  below_values <- unique(as.vector(colour_channels(below)))
  expect_true(
    length(below_values) == 1,
    info = paste0(
      "expected one flat value across the colour channels; got ",
      length(below_values), " (",
      paste(signif(sort(below_values), 4), collapse = ", "),
      ") from an image of dim ", paste(dim(below), collapse = "x")
    )
  )

  # ...and that flat value is the background being drawn rather than a colour
  # of plotZmap()'s own, which the assertion above does not distinguish. Tested
  # by ordering, not by value: rendering the same call over a darker background
  # must come out darker.
  #
  # The absolute value is deliberately *not* asserted. It is as
  # device-dependent as the alpha channel above -- cairo renders bgimage 0.5 to
  # 0.502, macOS to roughly 0.573, a colour-management difference -- so pinning
  # it repeats exactly the mistake this helper exists to avoid. That version was
  # written and did fail on macOS (test-plotZmap.R:106, 2026-07-28). Ordering
  # survives any monotone transfer function, which is all a sane one can be.
  darker <- render_zmap(tmp, matrix(1, 8, 8), "below_darker", bg = 0.25)
  darker_values <- unique(as.vector(colour_channels(darker)))
  expect_length(darker_values, 1)
  expect_lt(darker_values[1], below_values[1])

  # ...and when the z-map does clear the threshold, the output is no longer the
  # background. Without this the assertion above is satisfied by a function
  # that never draws anything at all.
  expect_false(isTRUE(all.equal(above, below)))
})

test_that("the threshold is applied per pixel, not to the map as a whole", {
  # Half the z-map above threshold, half below. The painted half must match a
  # fully-painted render and the unpainted half must match a bare background,
  # which pins *where* the thresholding lands rather than just that some of it
  # happened.
  tmp <- withr::local_tempdir()

  split <- matrix(1, 8, 8)
  split[, 1:4] <- 5

  mixed <- render_zmap(tmp, split, "mixed")
  above <- render_zmap(tmp, matrix(5, 8, 8), "all_above")
  below <- render_zmap(tmp, matrix(1, 8, 8), "all_below")

  expect_equal(mixed[, 1:32, ], above[, 1:32, ])
  expect_equal(mixed[, 33:64, ], below[, 33:64, ])
})

test_that("z-scores beyond the negative threshold are painted too", {
  # The threshold is applied to abs(zmap), so a strongly negative region is
  # signal and must be drawn. A `zmap < threshold` implementation would drop
  # it, and nothing else in the suite would notice.
  tmp <- withr::local_tempdir()

  negative <- render_zmap(tmp, matrix(-5, 8, 8), "negative")
  below <- render_zmap(tmp, matrix(1, 8, 8), "neutral")

  expect_false(isTRUE(all.equal(negative, below)))
})

test_that("raising the threshold removes regions from the z-map", {
  # Monotonicity: the same z-map drawn at a stricter threshold must lose the
  # region between the two thresholds and fall back to background there.
  tmp <- withr::local_tempdir()
  zmap <- matrix(4, 8, 8)

  lenient <- render_zmap(tmp, zmap, "lenient", threshold = 3)
  strict <- render_zmap(tmp, zmap, "strict", threshold = 10)
  background <- render_zmap(tmp, matrix(0, 8, 8), "background", threshold = 3)

  expect_false(isTRUE(all.equal(lenient, background)))
  expect_equal(strict, background)
})

test_that("a masked region is dropped from the z-map", {
  # Regression test for plotZmap(mask=). Commit 18e07cb (2016) landed the
  # mask import as "todo: applying the mask" and the todo was never picked up,
  # so `mask` was validated and then discarded in every released version. Half
  # a fully-painted z-map is masked here: that half must fall back to bare
  # background while the other half stays painted, which pins *where* the
  # masking lands rather than only that something changed.
  tmp <- withr::local_tempdir()

  # 0 = masked, matching generateCI()'s applyMask() for both matrices and PNGs.
  mask <- matrix(1, 8, 8)
  mask[, 1:4] <- 0

  masked <- render_zmap(tmp, matrix(5, 8, 8), "masked", mask = mask)
  above <- render_zmap(tmp, matrix(5, 8, 8), "unmasked")
  background <- render_zmap(tmp, matrix(1, 8, 8), "bare")

  expect_equal(masked[, 1:32, ], background[, 1:32, ])
  expect_equal(masked[, 33:64, ], above[, 33:64, ])
})

test_that("an all-zero mask blanks the z-map entirely", {
  # The boolean conversion collapsed every cell to FALSE regardless of input
  # (`mask[mask == 0] <- TRUE` coerces TRUE to 1, so the following
  # `mask[mask == 1] <- FALSE` unset it again). A mask that masks everything is
  # the case that distinguishes a working conversion from that one: under the
  # old code it masked nothing.
  tmp <- withr::local_tempdir()

  masked <- render_zmap(tmp, matrix(5, 8, 8), "all", mask = matrix(0, 8, 8))
  background <- render_zmap(tmp, matrix(0, 8, 8), "none")

  expect_equal(masked, background)
})

test_that("a matrix mask and a PNG mask select the same region", {
  # Both input paths use one convention (0 = masked), so the same mask
  # expressed either way must produce the same image. Asserted rather than
  # assumed: the roxygen claimed a matrix used the opposite polarity to a PNG,
  # and implementing that would have made a mask mask complementary halves in
  # plotZmap() and generateCI().
  tmp <- withr::local_tempdir()

  mask_file <- file.path(tmp, "mask.png")
  mask_image <- matrix(1, 8, 8)
  mask_image[, 1:4] <- 0
  png::writePNG(mask_image, mask_file)

  as_matrix <- render_zmap(tmp, matrix(5, 8, 8), "as_matrix", mask = mask_image)
  as_png <- render_zmap(tmp, matrix(5, 8, 8), "as_png", mask = mask_file)
  expect_equal(as_matrix, as_png)

  # ...and the shared region is the black half, not its complement.
  above <- render_zmap(tmp, matrix(5, 8, 8), "png_unmasked")
  background <- render_zmap(tmp, matrix(1, 8, 8), "png_bare")

  expect_equal(as_png[, 1:32, ], background[, 1:32, ])
  expect_equal(as_png[, 33:64, ], above[, 33:64, ])
})

test_that("plotZmap and generateCI mask the same half of one mask", {
  # The two functions have separate mask-import code (issue #89), so a mask
  # built for one has to be checked against the other or they can drift into
  # opposite polarities -- which is what their documentation asserted until
  # 1.1.0, and what an earlier attempt at this fix implemented. Comparing the
  # rendered z-map against the CI directly is not possible (one is a PNG, the
  # other a matrix), so both are reduced to "which half went away".
  mask <- matrix(1, 8, 8)
  mask[, 1:4] <- 0

  # generateCI(): masked cells become NA in the CI itself.
  ci_masked <- rcicr:::applyMask(matrix(1, 8, 8), mask, img_size = 8)
  ci_dropped_left <- all(is.na(ci_masked[, 1:4])) && !anyNA(ci_masked[, 5:8])

  # plotZmap(): masked cells fall back to the background in the rendered PNG.
  tmp <- withr::local_tempdir()
  zmap_masked <- render_zmap(tmp, matrix(5, 8, 8), "cross", mask = mask)
  background <- render_zmap(tmp, matrix(1, 8, 8), "cross_bare")
  painted <- render_zmap(tmp, matrix(5, 8, 8), "cross_full")
  zmap_dropped_left <-
    isTRUE(all.equal(zmap_masked[, 1:32, ], background[, 1:32, ])) &&
    isTRUE(all.equal(zmap_masked[, 33:64, ], painted[, 33:64, ]))

  expect_true(ci_dropped_left)
  expect_true(zmap_dropped_left)
})

test_that("mismatched mask and zmap dimensions error", {
  tmp <- withr::local_tempdir()

  expect_error(
    plotZmap(
      zmap = matrix(0, 8, 8), mask = matrix(0, 4, 4), sigma = 3,
      targetpath = tmp, size = 200
    ),
    "same dimensions"
  )
})

test_that("plotZmap rejects a mask that is neither a string nor a matrix", {
  tmp <- withr::local_tempdir()

  expect_error(
    plotZmap(
      zmap = matrix(0, 8, 8), mask = TRUE, sigma = 3,
      targetpath = tmp, size = 200
    ),
    "neither a string nor a matrix"
  )
})

test_that("the z-map is drawn in the orientation raster::plot() used", {
  # #186 swapped raster::plot() for graphics::image(), which lays matrix rows # nolint: commented_code_linter.
  # along x and counts y upward. Drawn naively the z-map comes out transposed
  # and vertically flipped -- over a base face, drawn separately by
  # rasterImage(), that is a silently wrong figure rather than an error.
  #
  # A suprathreshold blob in one corner of the matrix has to land in the
  # corresponding corner of the PNG. Row 1 is the TOP, column 1 the LEFT.
  tmp <- withr::local_tempdir()
  zmap <- matrix(0, 8, 8)
  zmap[1:2, 1:2] <- 9

  img <- render_zmap(tmp, zmap, "corner")[, , 1]
  n <- nrow(img)
  quadrant <- function(rows, cols) mean(img[rows, cols])

  top_left <- quadrant(1:(n / 2), 1:(n / 2))
  others <- c(
    quadrant(1:(n / 2), (n / 2 + 1):n),
    quadrant((n / 2 + 1):n, 1:(n / 2)),
    quadrant((n / 2 + 1):n, (n / 2 + 1):n)
  )

  # Stated relative to the other quadrants, never against an absolute grey.
  # The background here is a 0.5 matrix, but what a device paints for it is not
  # 0.5 everywhere: quartz renders that mid-grey at roughly 0.573 where cairo
  # gives 0.502, which is documented in ?plotZmap and failed this assertion on
  # the macOS job when it was first written against the literal value.
  expect_false(isTRUE(all.equal(top_left, others[1], tolerance = 0.05)))
  expect_equal(others, rep(others[1], 3), tolerance = 0.05)
})

test_that("raster is no longer a dependency", {
  expect_false("raster" %in% names(getNamespaceImports(asNamespace("rcicr"))))

  desc <- read.dcf(system.file("DESCRIPTION", package = "rcicr"))
  expect_false(grepl("\\braster\\b", paste(desc[1, c("Imports", "Depends")], collapse = " ")))
})

test_that("the rendered z-map matches what raster::plot() drew, pixel for pixel", {
  # The equivalence test for #186, which replaced raster::plot() with
  # graphics::image(). The reference PNGs in fixtures/ were rendered by the
  # raster implementation (1.2.3.9000, commit 437f755) and committed, because
  # once the dependency is gone they cannot be regenerated here.
  #
  # This is the only check that the *data* survived the swap: value-to-colour
  # mapping, cell geometry and orientation all at once. The release gate cannot
  # see it -- that harness captures ci$zmap, the matrix, computed before
  # anything is drawn -- so a rendering regression would otherwise pass green.
  #
  # The comparison comes in two halves, because a PNG's absolute pixel values
  # are a property of the graphics device and not of what was drawn -- the
  # fixtures were rendered on cairo, and quartz paints the same mid-grey at
  # roughly 0.573 against cairo's 0.502 (see ?plotZmap). Asserting exact values
  # everywhere failed the macOS job on correct output.
  #
  #   1. Structural, and device-independent: per-cell intensities must correlate
  #      with the reference. A colour-management difference shifts values
  #      monotonically and leaves this near 1; a transpose, a flip, a half-cell
  #      offset or a different palette moves whole cells and destroys it.
  #   2. Exact, where the fixture's own backend is in use. Tolerance 0.02 is
  #      ~5/255, covering 8-bit rounding only.

  input <- readRDS(test_path("fixtures", "zmap-raster-reference-input.rds"))
  tmp <- withr::local_tempdir()

  # Mean intensity of each z-map cell, the unit the two backends have to agree
  # on. The rendered PNG is 64px for an 8x8 matrix, so each cell is 8x8 pixels.
  cell_means <- function(img) {
    k <- nrow(img) / nrow(input$zmap)
    idx <- seq_len(nrow(input$zmap))
    outer(idx, idx, Vectorize(function(i, j) {
      mean(img[((i - 1) * k + 1):(i * k), ((j - 1) * k + 1):(j * k), 1])
    }))
  }

  for (case in c("nobg", "bg")) {
    plotZmap(
      zmap = input$zmap,
      bgimage = if (case == "bg") input$bg else "",
      sigma = 3, threshold = 3, decoration = FALSE,
      targetpath = tmp, filename = case, size = 64
    )

    got <- png::readPNG(file.path(tmp, paste0(case, ".png")))
    want <- png::readPNG(test_path("fixtures", paste0("zmap-raster-reference-", case, ".png")))

    # Drop any alpha plane: whether one is written is a property of the device,
    # not of what was drawn (see the note on colour channels above).
    got <- got[, , 1:3, drop = FALSE]
    want <- want[, , 1:3, drop = FALSE]

    expect_equal(dim(got), dim(want))
    expect_gt(cor(as.vector(cell_means(got)), as.vector(cell_means(want))), 0.99)

    if (!identical(Sys.info()[["sysname"]], "Darwin")) {
      expect_lt(max(abs(got - want)), 0.02)
    }
  }
})

test_that("a fixed colour scale passed through ... is respected", {
  # Regression test for #186. zlim is a valid argument to both the raster plot
  # method and graphics::image(), and plotZmap(zlim = c(-5, 5)) worked before
  # the swap -- the old code passed no zlim of its own. The replacement computed
  # one from the data and handed image() both, so the call died with
  # `formal argument "zlim" matched by multiple actual arguments`.
  #
  # Defaults in drawZmapLayer() are now merged by name, so this covers the whole
  # class rather than zlim alone: a caller who names any of them replaces it.
  tmp <- withr::local_tempdir()
  zmap <- matrix(c(rep(9, 32), rep(-9, 32)), 8, 8)

  expect_no_error(
    plotZmap(zmap, sigma = 3, threshold = 3, decoration = TRUE,
             targetpath = tmp, filename = "fixed_scale", size = 300,
             zlim = c(-20, 20))
  )

  # Honoured, not merely tolerated: a wider scale puts the same z-scores in a
  # different part of the palette, so the image has to differ.
  plotZmap(zmap, sigma = 3, threshold = 3, decoration = TRUE,
           targetpath = tmp, filename = "data_scale", size = 300)
  expect_false(isTRUE(all.equal(
    png::readPNG(file.path(tmp, "fixed_scale.png")),
    png::readPNG(file.path(tmp, "data_scale.png"))
  )))
})

test_that("custom breaks decide which colour each z-score gets", {
  # graphics::image() gives breaks precedence over zlim when assigning colours,
  # so the colour bar is drawn on them too, or it reports a scale the map was
  # not drawn on.
  #
  # Everything but the data is held fixed -- breaks, palette, and `main`, whose
  # default embeds the differing filename -- so the map is the only thing that
  # can differ between the two renders. An assertion about the pixels of a
  # single render would measure the graphics device instead: absolute colours,
  # and even counts of them, vary across cairo, quartz and Windows.
  tmp <- withr::local_tempdir()

  render <- function(name, value) {
    plotZmap(matrix(value, 8, 8), sigma = 3, threshold = 3, decoration = TRUE,
             targetpath = tmp, filename = name, size = 300, main = "z-map",
             col = c("blue", "red"), breaks = c(-100, 6, 100))
    png::readPNG(file.path(tmp, paste0(name, ".png")))
  }

  # 4 and 9 straddle the break at 6. Ignoring it makes each a constant map over
  # its own degenerate range, which lands every cell in the same colour.
  expect_false(isTRUE(all.equal(render("below_break", 4), render("above_break", 9))))
})

test_that("a palette passed through ... is used, not rejected", {
  # ?plotZmap has always offered `col` as the way to change the palette, and
  # passing it has always failed: the call supplied its own col alongside the
  # caller's, so R stopped with "formal argument 'col' matched by multiple
  # actual arguments" before drawing. Verified against the raster
  # implementation too -- it stacked the arguments the same way, so this is a
  # long-standing bug rather than something #186 introduced.
  tmp <- withr::local_tempdir()
  zmap <- matrix(c(rep(9, 32), rep(-9, 32)), 8, 8)
  palette <- grDevices::heat.colors(20)

  expect_no_error(
    plotZmap(zmap, sigma = 3, threshold = 3, decoration = TRUE,
             targetpath = tmp, filename = "custom", size = 300, col = palette)
  )
  expect_no_error(
    plotZmap(zmap, bgimage = matrix(0.5, 8, 8), sigma = 3, threshold = 3,
             decoration = TRUE, targetpath = tmp, filename = "custom_bg",
             size = 300, col = palette)
  )

  # The palette is honoured rather than merely tolerated: rendering with a
  # different one has to produce a different image.
  other <- file.path(tmp, "default.png")
  plotZmap(zmap, sigma = 3, threshold = 3, decoration = TRUE,
           targetpath = tmp, filename = "default", size = 300)
  expect_false(isTRUE(all.equal(
    png::readPNG(file.path(tmp, "custom.png")),
    png::readPNG(other)
  )))
})

test_that("a decorated z-map too small for its margins is refused, by name", {
  # Base R stops in plot.new() with `figure margins too large`, which names
  # neither this package nor the way out. zmapdecoration = TRUE is the default
  # and generateCI() sizes the device to img_size, so a small stimulus set hit
  # this on an ordinary call.
  tmp <- withr::local_tempdir()

  expect_error(
    plotZmap(matrix(5, 8, 8), sigma = 3, threshold = 3, decoration = TRUE,
             targetpath = tmp, filename = "toosmall", size = 128),
    "decoration = FALSE"
  )

  # png() creates the file when it opens the device, so refusing after that
  # point leaves a stub where a figure should be unless it is cleaned up.
  expect_false(file.exists(file.path(tmp, "toosmall.png")))
})

test_that("a smaller pointsize fits the decoration onto a small device", {
  # 6, not 8: the minimum size is measured in inches, so it depends on the
  # device's resolution. Windows' png() is 96 ppi where cairo is 72, and 128px
  # at pointsize 8 fits the first and not the second.
  tmp <- withr::local_tempdir()

  plotZmap(matrix(5, 8, 8), sigma = 3, threshold = 3, decoration = TRUE,
           targetpath = tmp, filename = "small", size = 128, pointsize = 6)

  expect_true(file.exists(file.path(tmp, "small.png")))
})

test_that("generateCI forwards zmappointsize to the z-map", {
  # The path that matters: generateCI() sizes the z-map device to img_size, so
  # without this argument a small stimulus set cannot produce a decorated z-map
  # at all.
  tmp <- withr::local_tempdir()
  # 128px: too small for the decoration at the default pointsize, large enough
  # to carry it at a smaller one. The default 32px fixture cannot fit it at any
  # sane pointsize, so it could not show the argument working.
  rdata <- make_fixture_rdata(tmp, img_size = 128)
  # The z-map goes to its own directory: the fixture writes a base.png into tmp,
  # and plotZmap() names the z-map after the base image, so a file.exists()
  # check against tmp would pass whether or not a z-map was ever drawn.
  zmaps <- file.path(tmp, "zmaps")
  args <- list(stimuli = 1:6, responses = c(1, -1, 1, -1, 1, -1),
               baseimage = "base", rdata = rdata, save_as_png = FALSE,
               targetpath = tmp, zmap = TRUE, zmapdecoration = TRUE,
               zmaptargetpath = zmaps)

  expect_error(do.call(generateCI, args), "decoration = FALSE")
  expect_false(file.exists(file.path(zmaps, "base.png")))

  expect_no_error(do.call(generateCI, c(args, list(zmappointsize = 6))))
  expect_true(file.exists(file.path(zmaps, "base.png")))
})

test_that("plotZmap treats a scalar NA mask as no mask, like generateCI", {
  # plotZmap() guarded applyMask() with !is.null(mask) while generateCI() used
  # hasMask(), so the same sentinel meant opposite things: NA reached
  # applyMask() here and died on "neither a string nor a matrix", where
  # generateCI() -- whose mask argument *defaults* to NA -- read it as "no
  # mask" (issue #246). NaN takes the same path and was equally rejected.
  #
  # Asserted as a relationship between renders, never absolute pixels: the
  # sentinel renders must equal the no-mask render and differ from a real
  # mask's, so the test cannot pass on a plotZmap() that masks everything.
  tmp <- withr::local_tempdir()

  half <- matrix(1, 8, 8)
  half[, 1:4] <- 0

  none <- render_zmap(tmp, matrix(5, 8, 8), "sentinel_none")
  masked <- render_zmap(tmp, matrix(5, 8, 8), "sentinel_half", mask = half)

  expect_equal(render_zmap(tmp, matrix(5, 8, 8), "sentinel_na", mask = NA), none)
  expect_equal(render_zmap(tmp, matrix(5, 8, 8), "sentinel_nan", mask = NaN), none)
  expect_false(isTRUE(all.equal(masked, none)))
})

test_that("a one-cell NA matrix is validated as a mask, not read as the sentinel", {
  # hasMask() collapsed to length(mask) == 1L && is.na(mask), which matrix(NA,
  # 1, 1) satisfies -- so a malformed one-cell mask was indistinguishable from
  # the scalar NA sentinel and silently discarded. generateCI() returned a
  # fully unmasked CI with no error or warning; plotZmap() only escaped
  # because its own guard was the inconsistent one this change removes.
  # Both must reach applyMask() and be rejected there.
  tmp <- withr::local_tempdir()
  one_cell <- matrix(NA, 1, 1)

  expect_true(rcicr:::hasMask(one_cell))

  # Against an 8x8 target the size check fires first; at matching size the
  # binary check does. Both are the validation this test exists to keep.
  expect_error(
    render_zmap(tmp, matrix(5, 8, 8), "onecell", mask = one_cell),
    "not of the same dimensions"
  )
  expect_error(
    rcicr:::applyMask(matrix(1, 1, 1), one_cell, img_size = 1),
    "values other than 0 or 1"
  )

  # The sentinels it must not swallow, and the valid one-cell mask it must not
  # start rejecting.
  expect_false(rcicr:::hasMask(NA))
  expect_false(rcicr:::hasMask(NaN))
  expect_false(rcicr:::hasMask(NULL))
  expect_false(rcicr:::hasMask(NA_character_))
  expect_true(rcicr:::hasMask(matrix(1, 1, 1)))
  expect_true(rcicr:::hasMask(matrix(NA, 2, 2)))

  # A dim attribute is not the only thing length-1 NA-ness fails to exclude:
  # list(NA) is length 1 and is.na() too, with no dim to catch it. The sentinel
  # is an *atomic* scalar, so malformed containers reach validation as well.
  expect_true(rcicr:::hasMask(list(NA)))
  expect_true(rcicr:::hasMask(array(NA)))
  expect_error(
    render_zmap(tmp, matrix(5, 8, 8), "listna", mask = list(NA)),
    "neither a string nor a matrix"
  )
})
