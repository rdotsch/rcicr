# Synthetic square grayscale PNG base face (never a real photo, avoids licensing concerns).
make_square_png <- function(path, size = 32, seed = 1) {
  withr::with_seed(seed, {
    img <- matrix(runif(size * size), size, size)
  })
  png::writePNG(img, path)
  invisible(path)
}

# Runs generateStimuli2IFC() at minimal size and returns the path to the .Rdata file it writes.
make_fixture_rdata <- function(dir, img_size = 32, n_trials = 6, nscales = 1, seed = 1) {
  png_path <- file.path(dir, "base.png")
  make_square_png(png_path, size = img_size, seed = seed)

  suppressWarnings(
    generateStimuli2IFC(
      base_face_files = list(base = png_path),
      n_trials = n_trials,
      img_size = img_size,
      stimulus_path = dir,
      seed = seed,
      ncores = 1,
      nscales = nscales,
      save_as_png = FALSE,
      save_rdata = TRUE
    )
  )

  rdata_files <- list.files(dir, pattern = "\\.Rdata$", full.names = TRUE)
  rdata_files[1]
}

# Pre-seeds a `reference_norms` vector into an existing .Rdata fixture, so
# computeInfoVal2IFC() takes the "already have a reference distribution" path
# and never reaches generateReferenceDistribution2IFC() or yesno().
seed_reference_norms <- function(rdata_path, n = 50, seed = 1) {
  e <- new.env()
  load(rdata_path, envir = e)
  e$reference_norms <- withr::with_seed(seed, rnorm(n, mean = 10, sd = 2))
  save(list = ls(e), file = rdata_path, envir = e)
  invisible(rdata_path)
}

# Rewrites an existing .Rdata fixture: adds or replaces the objects passed via
# `...`, and drops the names listed in `.remove`. Both directions are needed --
# the load()-collision guards are reached by *adding* a colliding name, and the
# "file did not contain X" guards by *removing* one.
#
# The leading dot on `.path` matters. R partially matches named arguments to
# formals, so a formal called `rdata_path` would swallow a planted
# `rdata = <decoy>` by prefix and make the helper open the decoy instead of the
# fixture. That cost a confusing test failure once already.
mutate_rdata <- function(.path, ..., .remove = character()) {
  e <- new.env()
  load(.path, envir = e)

  objs <- list(...)
  for (nm in names(objs)) assign(nm, objs[[nm]], envir = e)

  present <- intersect(.remove, ls(e))
  if (length(present)) rm(list = present, envir = e)

  save(list = ls(e), file = .path, envir = e)
  invisible(.path)
}
