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

# A set of distinct synthetic faces for the latent module, which needs several
# images rather than one to have a face space at all. Each gets its own seed, so
# the set has variation to decompose; identical images leave nothing to build.
make_face_set <- function(dir, n = 6, size = 16) {
  files <- file.path(dir, sprintf('face%02d.png', seq_len(n)))
  for (i in seq_along(files)) {
    make_square_png(files[i], size = size, seed = i)
  }
  files
}

# A simulated 2IFC observer for the latent module.
#
# The observer has a hidden target latent and, on each trial, picks whichever of
# the two rendered faces is closer to it. In latent space that reduces to the
# sign of the perturbation's projection onto target - base, which is what makes
# the recovered direction predictable rather than merely plausible: under
# Gaussian perturbations of covariance S the response-weighted mean estimates
# S %*% (target - base), not (target - base) itself.
simulate_latent_observer <- function(latent_params, base_latent, target_latent) {
  preference <- target_latent - base_latent
  projection <- as.vector(latent_params %*% preference)

  ifelse(projection > 0, 1, -1)
}

# What the response-weighted mean converges to for that observer: the preference
# direction weighted by the sampling covariance.
expected_latent_direction <- function(generator, base_latent, target_latent, latent_sigma = 1) {
  (latent_sigma * generator$latent_sd)^2 * (target_latent - base_latent)
}

# A whole latent pipeline in one call: a generator, a stimulus set and the
# .Rdata that links them. Lives here rather than in a test file so it is visible
# to every test that needs one, and so it can be extended in one place.
latent_fixture <- function(env = parent.frame(), n_trials = 20, n_components = 3, n_faces = 6) {
  dir <- withr::local_tempdir(.local_envir = env)
  out <- withr::local_tempdir(.local_envir = env)
  generator <- latentGeneratorPCA(make_face_set(dir, n = n_faces),
                                  n_components = n_components, img_size = 16)
  stimuli <- generateStimuliLatent2IFC(generator, n_trials = n_trials,
                                       stimulus_path = out, seed = 1,
                                       save_as_png = FALSE)
  list(generator = generator, stimuli = stimuli, out = out)
}

# A generator with enough dimensions for a recovery target to be interesting.
recovery_generator <- function(env = parent.frame(), n_components = 6, n_faces = 12, size = 16) {
  dir <- withr::local_tempdir(.local_envir = env)
  latentGeneratorPCA(make_face_set(dir, n = n_faces, size = size),
                     n_components = n_components, img_size = size)
}
