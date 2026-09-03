#!/usr/bin/env Rscript
#
# Reproducibility harness: runs one fixed battery of rcicr calls against whatever
# version of the package is on .libPaths(), and saves the numeric results.
#
# This file is always taken from the *current* working tree, so the reference
# version and the version under test execute identical code here -- only the
# package they call differs. Do not make it depend on anything added after the
# oldest version it is ever pointed at (v1.0.1).
#
# Driven by tools/compare-release-output.R; not meant to be called directly.
#
#   Rscript tools/compare-harness.R <outdir> <basedir> <workdir>
#
# RCICR_COMPARE_REF_VERSION tells it which version is on the other side of the
# comparison. Some calls only became comparable once a crash was fixed, so the
# battery grows as the reference moves forward -- but both sides always run the
# same battery, chosen by the *older* of the two.
#
# One <config>.rds per configuration, so the driver can compare a pair at a time
# rather than holding both whole batteries in memory (the 512px noise basis
# alone is ~126 MB per side).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: compare-harness.R <outdir> <basedir> <workdir>")
outdir  <- args[[1]]
basedir <- normalizePath(args[[2]], mustWork = TRUE)
workdir <- normalizePath(args[[3]], mustWork = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outdir <- normalizePath(outdir, mustWork = TRUE)

suppressPackageStartupMessages(library(rcicr))
setwd(workdir)

# --- the battery ------------------------------------------------------------
#
# Coverage rationale: every setting below changes numeric output, so each is
# varied at least once, and every entry point researchers actually call is
# exercised at least once (`extras`).
#
# `defaults-512-sinusoid` is the configuration researchers ran, so it carries
# the headline claim and gets almost everything. The rest trade image size for
# breadth -- a difference in the maths does not care how big the image is, and
# 64/128px keeps the whole battery to minutes rather than hours.
#
# extras (per config, cheapest first):
#   ci2ifc       generateCI2IFC() -- the 2IFC entry point, a wrapper over generateCI()
#   subset       non-contiguous `stimuli` -- the stimulus-number -> parameter-row lookup
#   participants generateCI(participants=) -- the per-participant average path
#   individual_cis save_individual_cis=TRUE, MD5 per filename -- see below
#   batch        batchGenerateCI() + autoscale() -- the whole-batch scaling constant
#   cumulative   computeCumulativeCICorrelation()
#   zmap_quick   generateCI(zmap=TRUE, zmapmethod='quick')  -- includes plotZmap()
#   zmap_ttest   generateCI(zmap=TRUE, zmapmethod='t.test') -- one t-test per pixel
#   mask         generateCI(mask=) -- masked pixels come back NA        [>= 1.1.0]
#   mask_rgba    generateCI(mask=) with a 4-channel PNG mask            [>= 1.1.0]
#   zmap_plain   a z-map with zmapdecoration = FALSE                    [>= 1.1.0]
#
# zmap_plain is also the only way to cover a z-map below 512px here. Decoration
# needs margins a small device does not have: the released versions this gate
# compares against fail with "figure margins too large" below about 160px, and
# have no way to ask for less. The working tree can now draw one, via
# plotZmap(pointsize=) / generateCI(zmappointsize=) from #177, but a battery
# entry the *reference* version cannot execute aborts the comparison, so the
# small configs still leave the decorated z-map out.
#
# The bracketed extras are skipped when the reference version predates them,
# because the reference *crashes* on those calls rather than returning a
# different number: undecorated z-maps die in `if (bgimage != '')` on R >= 4.2,
# z-maps below ~512px in `plot.new()` with "figure margins too large", and
# `mask` in the same length > 1 condition. A fix that turns a crash into a
# number has no old value to compare against -- that is the test suite's job,
# not this gate's. Check the reference version can actually run any extra you
# add, and tag it in SINCE if it needs a floor.
#
# stimulus_pngs = TRUE additionally writes the stimuli to disk and MD5s them:
# the PNG files are what participants actually see, and they go through
# quantisation and clamping that the in-memory matrices do not.
#
# individual_cis is the one extra whose subject is a *filename*, not a number,
# so it is MD5-per-basename rather than a matrix. Its participant IDs are
# deliberately in non-sorted order: the labelling defect it covers (#261) only
# appears when appearance order and sorted order disagree, and `participants`
# above cannot show it -- rep(paste0("p", 1:4), ...) is already sorted. Every
# reference version writes the same *set* of names here, permuted rather than
# different, so the keys line up across runs and only the values move.

ALL_EXTRAS <- c("ci2ifc", "subset", "participants", "individual_cis", "batch",
                "cumulative", "zmap_quick", "zmap_ttest", "mask", "mask_rgba",
                "zmap_plain")

# Oldest reference version that can run each extra without crashing.
# mask_rgba shares mask's floor: 1.1.0 and 1.2.3 were both run against a
# 4-channel mask before it was added here, and each returned a mask rather
# than an error. A 2-channel mask has no floor that works -- every reference
# version indexes channel 3 unconditionally and dies -- so it stays out.
SINCE <- c(mask = "1.1.0", mask_rgba = "1.1.0", zmap_plain = "1.1.0")

REF_VERSION <- package_version(Sys.getenv("RCICR_COMPARE_REF_VERSION", "1.0.1"))
message("battery for reference version ", REF_VERSION)

runnable <- function(extras) {
  keep <- vapply(extras, function(e) is.na(SINCE[e]) || REF_VERSION >= package_version(SINCE[e]),
                 logical(1), USE.NAMES = FALSE)
  extras[keep]
}

cfg_ <- function(name, img_size, nscales = 5, noise_type = "sinusoid", sigma = 25,
                 n_trials = 20, n_base = 1, anti = FALSE, same_params = TRUE,
                 infoval_iter = NA, extras = character(), stimulus_pngs = FALSE,
                 base_type = "plain", maximize_baseimage_contrast = TRUE) {
  stopifnot(all(extras %in% ALL_EXTRAS))
  extras <- runnable(extras)
  as.list(environment())
}

CONFIGS <- list(
  # The default pipeline, at the size and settings the package documents. The
  # only config that can carry the z-map extras (see above).
  cfg_("defaults-512-sinusoid", img_size = 512, stimulus_pngs = TRUE,
       extras = ALL_EXTRAS),

  # Spatial scales: fewer and more than the default 5. nscales also decides how
  # many parameters a trial draws, so it moves the RNG stream too.
  cfg_("sinusoid-128-nscales3", img_size = 128, nscales = 3,
       extras = setdiff(ALL_EXTRAS, c("zmap_quick", "zmap_ttest"))),
  cfg_("sinusoid-128-nscales6", img_size = 128, nscales = 6, n_trials = 12),

  # Alpha changes the stored base face rather than the noise parameters or CI.
  # Cover both affected cases: varying alpha at the default contrast setting,
  # and constant alpha when contrast maximization is disabled.
  cfg_("sinusoid-64-alpha-varying", img_size = 64,
       base_type = "alpha_varying", stimulus_pngs = TRUE),
  cfg_("sinusoid-64-alpha-grey-varying", img_size = 64,
       base_type = "alpha_grey_varying", stimulus_pngs = TRUE),
  cfg_("sinusoid-64-alpha-opaque-default", img_size = 64,
       base_type = "alpha_opaque", stimulus_pngs = TRUE),
  cfg_("sinusoid-64-alpha-opaque-no-max", img_size = 64,
       base_type = "alpha_opaque", maximize_baseimage_contrast = FALSE,
       stimulus_pngs = TRUE),

  # Gabor noise, at two envelope widths.
  cfg_("gabor-128-sigma25", img_size = 128, noise_type = "gabor",
       extras = c("ci2ifc", "batch")),
  cfg_("gabor-128-sigma10", img_size = 128, noise_type = "gabor", sigma = 10),

  # Two base images with antiCI, sharing one parameter set...
  cfg_("sinusoid-128-twobase-anti", img_size = 128, n_base = 2, anti = TRUE,
       stimulus_pngs = TRUE, extras = c("ci2ifc", "batch", "cumulative")),
  # ...and two drawing independent parameters, which is a different RNG stream.
  cfg_("sinusoid-128-twobase-indep", img_size = 128, n_base = 2,
       same_params = FALSE, stimulus_pngs = TRUE),

  # InfoVal, at 64px where the reference distribution is affordable. nscales and
  # sigma are the two settings v1.0.1 failed to save into the .Rdata, so both
  # get a non-default run: that is where its InfoVal was computed against the
  # wrong null.
  cfg_("sinusoid-64-infoval", img_size = 64, infoval_iter = 200),
  cfg_("sinusoid-64-nscales3-infoval", img_size = 64, nscales = 3, infoval_iter = 200),
  cfg_("gabor-64-sigma10-infoval", img_size = 64, noise_type = "gabor", sigma = 10,
       infoval_iter = 200)
)

# --quick (RCICR_COMPARE_QUICK=1) drops the 512px config for local iteration.
# The release gate always runs the full battery.
if (nzchar(Sys.getenv("RCICR_COMPARE_QUICK"))) {
  CONFIGS <- Filter(function(cfg) cfg$img_size < 512, CONFIGS)
  message("quick mode: ", length(CONFIGS), " configs")
}

SCALINGS <- c("none", "constant", "matched", "independent")

run_config <- function(cfg) {
  out <- list()

  # 1. Noise basis, straight from the generator, and 2. one noise image from a
  #    fixed parameter vector, which isolates generateNoiseImage() from the RNG
  #    and from generateCI()'s weighting.
  #
  #    Both are parked on disk and dropped from memory before the pipeline runs:
  #    at 512px the basis is ~250 MB, every pipeline function loads its own copy
  #    from the .Rdata anyway, and a cluster worker gets a third. Holding all of
  #    them at once is what makes this fall over on an 8 GB machine -- the
  #    reference version's worker is killed and the run dies in
  #    `serialize(data, node$con): ignoring SIGPIPE signal`.
  basis_file <- file.path(workdir, paste0("basis-", cfg$name, ".rds"))
  local({
    p <- generateNoisePattern(img_size = cfg$img_size, nscales = cfg$nscales,
                              noise_type = cfg$noise_type, sigma = cfg$sigma)
    set.seed(2)
    params_one <- (runif(max(p$patchIdx)) * 2) - 1
    saveRDS(list(patches = p$patches, patchIdx = p$patchIdx,
                 noise_image = generateNoiseImage(params_one, p)),
            basis_file, compress = FALSE)
  })
  invisible(gc(verbose = FALSE))

  # 3. Stimulus generation. ncores = 1 because the reference version has no
  #    serial fast path, and a one-worker cluster is the cheapest thing both
  #    versions agree on.
  stim_dir <- file.path(getwd(), paste0("stim-", cfg$name))
  ci_dir   <- file.path(getwd(), paste0("ci-", cfg$name))
  zmap_dir <- file.path(getwd(), paste0("zmap-", cfg$name))
  unlink(c(stim_dir, ci_dir, zmap_dir, file.path(getwd(), "zmaps")), recursive = TRUE)

  base_files <- if (identical(cfg$base_type, "plain")) {
    sprintf("base%d_%d.png", seq_len(cfg$n_base), cfg$img_size)
  } else {
    stopifnot(cfg$n_base == 1L)
    sprintf("base_%s_%d.png", cfg$base_type, cfg$img_size)
  }
  bases <- stats::setNames(as.list(file.path(basedir, base_files)),
                           sprintf("base%d", seq_len(cfg$n_base)))

  generateStimuli2IFC(base_face_files = bases, n_trials = cfg$n_trials,
                      img_size = cfg$img_size, stimulus_path = stim_dir,
                      label = "cmp", seed = 1, noise_type = cfg$noise_type,
                      nscales = cfg$nscales, sigma = cfg$sigma, ncores = 1,
                      use_same_parameters = cfg$same_params,
                      maximize_baseimage_contrast = cfg$maximize_baseimage_contrast,
                      save_as_png = cfg$stimulus_pngs, save_rdata = TRUE)

  rdata <- list.files(stim_dir, pattern = "[.]Rdata$", full.names = TRUE)
  if (length(rdata) != 1L) stop("expected exactly one .Rdata in ", stim_dir)
  env <- new.env()
  load(rdata, envir = env)
  for (b in names(bases)) {
    out[[paste0("stimuli_params_", b)]] <- env$stimuli_params[[b]]
    out[[paste0("base_face_", b)]] <- env$base_faces[[b]]
  }

  # The stimuli as they reach a participant: 8-bit, quantised, clamped.
  if (cfg$stimulus_pngs) {
    pngs <- sort(list.files(stim_dir, pattern = "[.]png$", full.names = TRUE))
    if (!length(pngs)) stop("no stimulus PNGs written in ", stim_dir)
    if (identical(cfg$base_type, "plain")) {
      out$stimulus_pngs <- stats::setNames(unname(tools::md5sum(pngs)), basename(pngs))
    } else {
      out$stimulus_pngs <- simplify2array(lapply(pngs, png::readPNG))
      dimnames(out$stimulus_pngs)[[3]] <- basename(pngs)
    }
  }

  # 4. Classification images, one per scaling method.
  set.seed(20260728)
  responses <- sample(c(1, -1), cfg$n_trials, replace = TRUE)
  stimuli <- seq_len(cfg$n_trials)
  for (sc in SCALINGS) {
    ci <- generateCI(stimuli = stimuli, responses = responses,
                     baseimage = "base1", rdata = rdata, scaling = sc,
                     antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir)
    if (identical(sc, SCALINGS[[1]])) {
      out$ci <- ci$ci              # scaling-independent
      out$combined <- ci$combined  # base image + CI, what gets written to disk
    }
    out[[paste0("scaled_", sc)]] <- ci$scaled
  }

  # 5. The other entry points, on the configs that ask for them.
  ex <- cfg$extras

  if ("ci2ifc" %in% ex) {
    ci <- generateCI2IFC(stimuli = stimuli, responses = responses,
                         baseimage = "base1", rdata = rdata,
                         antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir)
    out$ci2ifc <- ci$ci
    out$ci2ifc_scaled <- ci$scaled
  }

  if ("subset" %in% ex) {
    # Non-contiguous, out of order, with a repeat: exercises the stimulus-number
    # to parameter-row lookup rather than a plain 1:n slice.
    idx <- c(seq(2, cfg$n_trials, by = 3), 1, 1)
    ci <- generateCI(stimuli = stimuli[idx], responses = responses[idx],
                     baseimage = "base1", rdata = rdata, scaling = "independent",
                     antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir)
    out$subset_ci <- ci$ci
    out$subset_scaled <- ci$scaled
  }

  if ("participants" %in% ex) {
    pids <- rep(paste0("p", 1:4), length.out = cfg$n_trials)
    ci <- generateCI(stimuli = stimuli, responses = responses, participants = pids,
                     baseimage = "base1", rdata = rdata, scaling = "independent",
                     antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir,
                     n_cores = 1)
    out$participants_ci <- ci$ci
    out$participants_scaled <- ci$scaled
  }

  if ("individual_cis" %in% ex) {
    # Appearance order c("pb", "pa"), sorted order c("pa", "pb"). Kept separate
    # from the `participants` extra above rather than making that one unsorted:
    # changing its IDs would regroup its trials and move participants_ci and
    # participants_scaled too, adding deviations that have nothing to do with
    # the filenames under test.
    ind_dir <- file.path(ci_dir, "individual")
    unlink(ind_dir, recursive = TRUE)
    pids <- rep(c("pb", "pa"), each = ceiling(cfg$n_trials / 2))[seq_len(cfg$n_trials)]
    generateCI(stimuli = stimuli, responses = responses, participants = pids,
               baseimage = "base1", rdata = rdata, scaling = "independent",
               antiCI = cfg$anti, save_as_png = FALSE, save_individual_cis = TRUE,
               targetpath = ind_dir, n_cores = 1)
    pngs <- sort(list.files(file.path(ind_dir, "individual_cis"),
                            pattern = "[.]png$", full.names = TRUE))
    if (!length(pngs)) stop("no individual-CI PNGs written in ", ind_dir)
    out$individual_cis <- stats::setNames(unname(tools::md5sum(pngs)), basename(pngs))
  }

  if ("batch" %in% ex) {
    df <- data.frame(stim = stimuli, resp = responses,
                     grp = rep(c("a", "b"), length.out = cfg$n_trials),
                     stringsAsFactors = FALSE)
    # scaling defaults to 'autoscale', so this covers autoscale()'s shared
    # constant as well as the per-group CIs.
    cis <- batchGenerateCI(data = df, by = "grp", stimuli = "stim",
                           responses = "resp", baseimage = "base1", rdata = rdata,
                           antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir)
    for (nm in sort(names(cis))) {
      out[[paste0("batch_", nm, "_ci")]] <- cis[[nm]]$ci
      out[[paste0("batch_", nm, "_scaled")]] <- cis[[nm]]$scaled
    }
  }

  if ("cumulative" %in% ex) {
    out$cumulative <- computeCumulativeCICorrelation(
      stimuli = stimuli, responses = responses, baseimage = "base1", rdata = rdata)
  }

  if ("mask" %in% ex) {
    ci <- generateCI(stimuli = stimuli, responses = responses,
                     baseimage = "base1", rdata = rdata, scaling = "independent",
                     antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir,
                     mask = file.path(basedir, sprintf("mask_%d.png", cfg$img_size)))
    out$mask_ci <- ci$ci
    out$mask_scaled <- ci$scaled
  }

  if ("mask_rgba" %in% ex) {
    ci <- generateCI(stimuli = stimuli, responses = responses,
                     baseimage = "base1", rdata = rdata, scaling = "independent",
                     antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir,
                     mask = file.path(basedir,
                                      sprintf("mask_rgba_%d.png", cfg$img_size)))
    out$mask_rgba_ci <- ci$ci
    out$mask_rgba_scaled <- ci$scaled
  }

  zmap_calls <- list(
    zmap_quick = list(method = "quick",  decoration = TRUE),
    zmap_ttest = list(method = "t.test", decoration = TRUE),
    zmap_plain = list(method = "quick",  decoration = FALSE)
  )
  for (key in intersect(names(zmap_calls), ex)) {
    call <- zmap_calls[[key]]
    ci <- generateCI(stimuli = stimuli, responses = responses,
                     baseimage = "base1", rdata = rdata, scaling = "independent",
                     antiCI = cfg$anti, save_as_png = FALSE, targetpath = ci_dir,
                     zmap = TRUE, zmapmethod = call$method,
                     zmapdecoration = call$decoration,
                     zmaptargetpath = zmap_dir, n_cores = 1)
    out[[key]] <- ci$zmap
  }

  # 6. InfoVal, where it is cheap enough. Deterministic given the .Rdata: the
  #    reference distribution inherits the stimulus seed (see DECISIONS.md).
  if (!is.na(cfg$infoval_iter)) {
    ci <- generateCI(stimuli = stimuli, responses = responses,
                     baseimage = "base1", rdata = rdata, scaling = "none",
                     save_as_png = FALSE, targetpath = ci_dir)
    out$infoval <- computeInfoVal2IFC(ci, rdata, iter = cfg$infoval_iter)
  }

  unlink(c(stim_dir, ci_dir, zmap_dir, file.path(getwd(), "zmaps")), recursive = TRUE)

  # Basis first, so the report reads generator -> pipeline.
  basis <- readRDS(basis_file)
  unlink(basis_file)
  c(basis, out)
}

saveRDS(list(version = as.character(utils::packageVersion("rcicr")),
             r_version = R.version.string,
             configs = vapply(CONFIGS, function(cfg) cfg$name, character(1))),
        file.path(outdir, "meta.rds"))

for (cfg in CONFIGS) {
  t0 <- Sys.time()
  message("  config: ", cfg$name)
  res <- run_config(cfg)
  # compress = FALSE: these run to hundreds of MB at 512px, and the file lives
  # in tempdir for seconds.
  saveRDS(res, file.path(outdir, paste0(cfg$name, ".rds")), compress = FALSE)
  rm(res)
  invisible(gc(verbose = FALSE))
  message("    done in ", format(round(difftime(Sys.time(), t0, units = "secs"), 1)))
}
message("wrote ", outdir)
