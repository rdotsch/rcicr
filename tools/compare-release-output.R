#!/usr/bin/env Rscript
#
# Release gate: compare this tree's numeric output against the last version of
# rcicr that researchers published results with, across the battery in
# tools/compare-harness.R.
#
# Researchers re-run old analysis scripts and publish what comes out, so the
# question this answers is the only one that matters at release time: does the
# tarball we are about to ship still produce the numbers the reference version
# produced?
#
# `tests/testthat/test-regression-baseline.R` cannot answer it. That golden
# master pins values this repository computed for itself; this runs the actual
# reference tree, installed from its own commit, and compares.
#
#   Rscript tools/compare-release-output.R [--ref=<git-rev>] [--quick]
#                                          [--install-deps] [--keep] [--verbose]
#
# Exit status:  0  every difference is expected and accounted for
#               1  an unexpected difference, or a stale expectation
#               2  the comparison could not be run (setup failure)
#
# A difference is only tolerated if it is listed in EXPECTED below, with a
# reason and the NEWS.md heading that documents it for users. Adding an entry is
# a deliberate, reviewable act -- which is the point. See RELEASING.md for
# where this sits in the release checklist.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default = NULL) {
  hit <- grep(paste0("^--", flag, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", flag, "="), "", hit[[1]]) else default
}
# v1.0.1 (b6ab269) is the last version published before 1.1.0 and the one nearly
# every result in the literature was produced with. It stays pinned there rather
# than advancing to each new release: moving it would let a tree drift away from
# the published numbers one tolerated epsilon at a time, each step "identical to
# the last release". See DECISIONS.md.
REF     <- opt("ref", "v1.0.1")
KEEP    <- "--keep"    %in% args
QUICK   <- "--quick"   %in% args
VERBOSE <- "--verbose" %in% args
# The reference version's DESCRIPTION lists packages this one dropped, and it
# will not install without them. --install-deps fetches the missing ones into a
# throwaway library (what CI does); without it the script says which they are
# and stops, rather than writing to someone's library uninvited.
INSTALL_DEPS <- "--install-deps" %in% args

# --- expected deviations ----------------------------------------------------
#
# key: "<config>/<output>", "<config>/*" for a whole config, or a vector of
#      either when one cause shows up in several outputs. Every key listed must
#      actually deviate, or the run fails as stale -- so a vector is a claim
#      about all of them, not a net cast wide.
# ref: which reference this applies to -- a deviation from v1.0.1 is not a
#      deviation from v1.1.0, and listing it for both would make one of the two
#      runs report a stale expectation. A defect present in both references is
#      the exception, and needs an entry per reference.
# check: optional predicate(ref_value, cur_value). Without one the key excuses
#      *any* deviation in that output, including a later, unrelated regression
#      in the same place. Supply one wherever the expected change has a shape
#      that can be tested.

# The shape of the individual-CI relabelling: the same images under corrected
# names. Filenames must be unchanged and the hashes a permutation of each other,
# which holds only if every pixel is untouched. A regression in individual
# scaling or in PNG writing changes a hash to one that is not in the other
# side's set, and so is reported rather than absorbed by the entries below.
relabelled_only <- function(a, b) {
  identical(names(a), names(b)) && identical(sort(unname(a)), sort(unname(b)))
}

same_files_changed <- function(a, b) {
  identical(names(a), names(b)) && length(a) > 0L && all(a != b)
}

alpha_expectations <- function(ref) {
  list(
    list(ref = ref,
         key = c("sinusoid-64-alpha-varying/base_face_base1",
                 "sinusoid-64-alpha-varying/combined"),
         reason = paste("The reference averaged varying alpha into the colour channels.",
                        "Dropping it changes the stored base face and the CI rendered over it,",
                        "while the noise parameters and raw CI remain identical."),
         news = "Reproducibility impact"),
    list(ref = ref, key = "sinusoid-64-alpha-varying/stimulus_pngs",
         reason = "The corrected varying-alpha base face changes every rendered stimulus.",
         check = same_files_changed,
         news = "Reproducibility impact"),
    list(ref = ref,
         key = c("sinusoid-64-alpha-opaque-no-max/base_face_base1",
                 "sinusoid-64-alpha-opaque-no-max/combined",
                 "sinusoid-64-alpha-opaque-no-max/scaled_matched"),
         reason = paste("With contrast maximization disabled, the reference averaged constant",
                        "alpha into the base. The corrected base changes the rendered CI and",
                        "the matched scaling range, but not the raw CI."),
         news = "Reproducibility impact"),
    list(ref = ref, key = "sinusoid-64-alpha-opaque-no-max/stimulus_pngs",
         reason = "The corrected opaque-alpha base face changes every rendered stimulus.",
         check = same_files_changed,
         news = "Reproducibility impact")
  )
}

EXPECTED <- c(list(
  list(ref = "v1.0.1", key = "sinusoid-64-nscales3-infoval/infoval",
       reason = paste("v1.0.1 did not save nscales/sigma into the .Rdata, so its",
                      "reference distribution was rebuilt at the default nscales = 5",
                      "regardless of how the stimuli were made. At nscales != 5 the old",
                      "InfoVal was computed against the wrong null and is simply wrong."),
       news = "Reproducibility impact"),
  list(ref = "v1.0.1", key = "gabor-64-sigma10-infoval/infoval",
       reason = paste("Same cause as the nscales case: v1.0.1's reference distribution",
                      "ignored noise_type and sigma, so InfoVal for Gabor noise (and for",
                      "any non-default sigma) was measured against a sinusoid null."),
       news = "Reproducibility impact"),

  list(ref = "v1.0.1",
       key = c("defaults-512-sinusoid/individual_cis",
               "sinusoid-128-nscales3/individual_cis"),
       reason = paste("Individual-CI PNGs were named from the participants' order of",
                      "appearance while the loop indexed them by sorted factor level, so",
                      "every file carried another participant's ID whenever the two orders",
                      "disagreed. The pixels are unchanged -- the same set of images is",
                      "written, under corrected names -- so what moves here is which MD5",
                      "sits under which filename. Only the individual-CI PNGs are affected;",
                      "participants_ci and participants_scaled are means across",
                      "participants and must still match."),
       check = relabelled_only,
       news = "Reproducibility impact"),

  # Listed for both references, unlike the entries around it, because the defect
  # is in both: the mislabelling dates to 0.4.0 and every tagged release since
  # carries it, so neither run can report this as stale.
  list(ref = "v1.2.3",
       key = c("defaults-512-sinusoid/individual_cis",
               "sinusoid-128-nscales3/individual_cis"),
       reason = paste("Same defect as the v1.0.1 entry above -- v1.2.3 names individual-CI",
                      "PNGs from order of appearance too. The pixels are unchanged; which",
                      "MD5 sits under which filename is what moves."),
       check = relabelled_only,
       news = "Reproducibility impact"),

  list(ref = "v1.1.0",
       key = c("defaults-512-sinusoid/zmap_quick",
               "defaults-512-sinusoid/zmap_plain",
               "sinusoid-128-nscales3/zmap_plain"),
       reason = paste("The .Rdata's noise sigma was overwriting generateCI()'s z-map blur",
                      "sigma via load(), so 1.1.0 blurred z-maps with 25 instead of 3 and",
                      "ignored the sigma the caller passed. Fixed here,",
                      "which changes every blur-based z-map made from a 1.1.0 stimulus set.",
                      "zmap_ttest is deliberately absent: that method does not blur, so it",
                      "must still match."),
       news = "Reproducibility impact")
), alpha_expectations("v1.0.1"), alpha_expectations("v1.1.0"),
   alpha_expectations("v1.2.3"), alpha_expectations("v1.3.0"))

# Entries are matched by position so a vector `key` can be tracked as one
# expectation: it has fired once any of its keys deviates.
expectation_for <- function(key) {
  cfg <- sub("/.*$", "", key)
  for (i in seq_along(EXPECTED)) {
    e <- EXPECTED[[i]]
    if (key %in% e$key || paste0(cfg, "/*") %in% e$key) return(c(e, list(id = i)))
  }
  NULL
}

# --- tolerances -------------------------------------------------------------
#
# Bit-identity is required wherever a difference could only come from a changed
# random stream, a changed algorithm or a changed file format, rather than from
# summation order. Everything else is allowed to move by a few ULP *of the
# values involved* -- an absolute epsilon would be far too tight for z-maps,
# whose values run to several units, and needlessly loose nowhere.
#
# ULPS = 8 rather than 1 because these are not single passes over the data. A
# z-map is a convolution over 262,144 pixels followed by a standardisation over
# the same 262,144 values, so a 3.5e-18 difference in the classification image
# arrives as ~1.3e-15 -- measured, at 512px. Even at 8 ULP the bar is ~1e-15,
# orders of magnitude below anything that could move a pixel at 8-bit or move a
# cell across the z threshold, and both of those are checked exactly: the NA
# pattern must match cell for cell, which is what caught the sigma bug (1,282
# cells changed side) while every numeric tolerance here would have passed it.

# stimulus_pngs and individual_cis are MD5 digests, not numbers: they belong on
# the bit-identical branch both because a hash has no meaningful tolerance and
# because the numeric branch below calls abs() on them.
EXACT_RE  <- "^(patchIdx|stimuli_params_|base_face_|stimulus_pngs|individual_cis)"
IMAGE_RE  <- "^(patches|noise_image|ci|combined|scaled_|ci2ifc|subset_|participants_|batch_)"
ULP       <- .Machine$double.eps      # 2.22e-16
ULPS      <- 8
INFOVAL_TOL <- 1e-9

say <- function(...) cat(..., "\n", sep = "")
die <- function(...) { say("ERROR: ", ...); quit(status = 2L) }

repo <- normalizePath(file.path(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])), ".."), mustWork = TRUE)
setwd(repo)

# Deliberately *beside* R's session tempdir, not inside it: R removes its own
# tempdir when the driver exits, which silently defeated --keep.
tmp <- file.path(dirname(tempdir()), paste0("rcicr-compare-", Sys.getpid()))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
worktree <- file.path(tmp, "ref-tree")
lib_ref  <- file.path(tmp, "lib-ref")
lib_cur  <- file.path(tmp, "lib-cur")
lib_deps <- file.path(tmp, "lib-deps")
basedir  <- file.path(tmp, "base")
for (d in c(lib_ref, lib_cur, lib_deps, basedir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Snapshot the harness, so both sides provably run the same battery even if the
# working copy is edited while a run is in flight.
harness <- file.path(tmp, "compare-harness.R")
if (!file.copy(file.path(repo, "tools", "compare-harness.R"), harness, overwrite = TRUE))
  die("could not copy tools/compare-harness.R into ", tmp)

cleanup <- function() {
  if (KEEP) { say("\nkept: ", tmp); return(invisible()) }
  system2("git", c("worktree", "remove", "--force", shQuote(worktree)),
          stdout = FALSE, stderr = FALSE)
  unlink(tmp, recursive = TRUE)
}
on.exit(cleanup(), add = TRUE)

# --- set up both versions ---------------------------------------------------

resolve <- function(rev) {
  sha <- suppressWarnings(system2("git", c("rev-parse", "--verify", "--quiet",
                                           paste0(rev, "^{commit}")), stdout = TRUE))
  if (length(sha) == 1L) sha else NA_character_
}
ref_sha <- resolve(REF)
if (is.na(ref_sha)) die("cannot resolve --ref=", REF)
say("reference : ", REF, " (", substr(ref_sha, 1, 7), ")")

for (e in EXPECTED) if (is.na(resolve(e$ref))) die("EXPECTED entry names an unresolvable ref: ", e$ref)
EXPECTED <- Filter(function(e) identical(resolve(e$ref), ref_sha), EXPECTED)
say("expected deviations on file for this reference: ", length(EXPECTED))
say("under test: ", system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE),
    if (length(system2("git", c("status", "--porcelain"), stdout = TRUE))) " + uncommitted changes" else "")
if (QUICK) say("mode      : QUICK -- 512px config skipped, NOT sufficient for a release")

if (system2("git", c("worktree", "add", "--detach", shQuote(worktree), REF),
            stdout = FALSE, stderr = FALSE) != 0L) die("git worktree add failed")

# Anything the reference version imports that this one no longer does.
ref_desc <- read.dcf(file.path(worktree, "DESCRIPTION"))
ref_version <- unname(ref_desc[1, "Version"])
ref_deps <- unlist(strsplit(paste(
  ref_desc[1, intersect(c("Depends", "Imports"), colnames(ref_desc))], collapse = ","), ","))
ref_deps <- trimws(sub("\\(.*", "", ref_deps))
missing_deps <- setdiff(ref_deps[nzchar(ref_deps)],
                        c("R", rownames(utils::installed.packages())))
if (length(missing_deps)) {
  if (!INSTALL_DEPS)
    die(REF, " needs packages this library does not have: ",
        paste(missing_deps, collapse = ", "),
        "\n       install them, or re-run with --install-deps to put them in a temporary library.")
  say("installing ", length(missing_deps), " reference-only dependencies: ",
      paste(missing_deps, collapse = ", "))
  repos <- getOption("repos")
  if (!length(repos) || any(repos == "@CRAN@")) repos <- c(CRAN = "https://cloud.r-project.org")
  utils::install.packages(missing_deps, lib = lib_deps, repos = repos, quiet = TRUE)
  still <- setdiff(missing_deps, rownames(utils::installed.packages(lib.loc = lib_deps)))
  if (length(still)) die("could not install: ", paste(still, collapse = ", "))
}

libpath <- function(...) paste(c(..., lib_deps, .libPaths()), collapse = .Platform$path.sep)

install_into <- function(src, lib) {
  say("installing ", basename(src), " -> ", basename(lib))
  logfile <- file.path(tmp, paste0("install-", basename(lib), ".log"))
  st <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "INSTALL", "--no-docs", "--no-help", "--no-byte-compile",
                  "-l", shQuote(lib), shQuote(src)),
                env = paste0("R_LIBS=", libpath(lib)),
                stdout = logfile, stderr = logfile)
  if (st != 0L) {
    writeLines(tail(readLines(logfile, warn = FALSE), 30))
    die("R CMD INSTALL failed for ", src)
  }
}
install_into(worktree, lib_ref)
install_into(repo, lib_cur)

# Deterministic synthetic base images -- never a real photograph.
make_base <- function(n, path, shift = 0) {
  ax <- seq(-1, 1, length.out = n)
  g <- outer(ax, ax, function(x, y) exp(-((x - shift)^2 + y^2) / 0.5))
  png::writePNG(g, path)
}
make_base_rgba <- function(n, path, varying_alpha) {
  ax <- seq(-1, 1, length.out = n)
  g <- outer(ax, ax, function(x, y) exp(-(x^2 + y^2) / 0.5))
  rgba <- array(0, dim = c(n, n, 4))
  for (i in 1:3) rgba[, , i] <- g
  rgba[, , 4] <- if (varying_alpha) {
    outer(ax, ax, function(x, y) as.numeric(x^2 + y^2 <= 0.6))
  } else {
    1
  }
  png::writePNG(rgba, path)
}
# A mask has to be exactly 0/1, and 0 (black) is the region masked away.
mask_disc <- function(n) {
  ax <- seq(-1, 1, length.out = n)
  outer(ax, ax, function(x, y) as.numeric(x^2 + y^2 <= 0.6))
}
make_mask <- function(n, path) {
  png::writePNG(mask_disc(n), path)
}
# The same disc as an RGBA image. png::readPNG() hands a single-channel PNG
# straight back as a 2D matrix, so the plain mask above never reaches the
# channel-collapse branch of the importer -- this one is the only fixture here
# that does. Alpha is deliberately opaque rather than a copy of the disc: the
# reference versions compare channels 1-3 only, and an alpha plane that differs
# from them must not change the result on either side.
make_mask_rgba <- function(n, path) {
  disc <- mask_disc(n)
  rgba <- array(0, dim = c(n, n, 4))
  for (i in 1:3) rgba[, , i] <- disc
  rgba[, , 4] <- 1
  png::writePNG(rgba, path)
}
for (sz in c(64, 128, 512)) {
  make_base(sz, file.path(basedir, sprintf("base1_%d.png", sz)))
  make_base(sz, file.path(basedir, sprintf("base2_%d.png", sz)), shift = 0.3)
  make_base_rgba(sz, file.path(basedir, sprintf("base_alpha_varying_%d.png", sz)), TRUE)
  make_base_rgba(sz, file.path(basedir, sprintf("base_alpha_opaque_%d.png", sz)), FALSE)
  make_mask(sz, file.path(basedir, sprintf("mask_%d.png", sz)))
  make_mask_rgba(sz, file.path(basedir, sprintf("mask_rgba_%d.png", sz)))
}

run_side <- function(lib, tag) {
  out <- file.path(tmp, paste0("out-", tag))
  wd  <- file.path(tmp, paste0("work-", tag))
  dir.create(wd, showWarnings = FALSE)
  logfile <- file.path(tmp, paste0("harness-", tag, ".log"))
  say("running battery: ", tag, "  (log: ", logfile, ")")
  libs <- libpath(lib)
  st <- system2(file.path(R.home("bin"), "Rscript"),
                c(shQuote(harness), shQuote(out), shQuote(basedir), shQuote(wd)),
                # The battery is chosen by the *reference* version on both sides:
                # some calls only became comparable once a crash was fixed.
                env = c(paste0("R_LIBS=", libs), paste0("R_LIBS_USER=", libs),
                        paste0("RCICR_COMPARE_REF_VERSION=", ref_version),
                        # The 512px battery is memory-bound, not CPU-bound: the
                        # noise basis is ~250 MB and the reference version's
                        # cluster worker gets its own copy, so R's default
                        # eagerness to grow the heap is what tips an 8 GB
                        # machine over. Trading some speed for a smaller
                        # footprint is the right way round here.
                        "R_GC_MEM_GROW=0",
                        if (QUICK) "RCICR_COMPARE_QUICK=1"),
                stdout = if (VERBOSE) "" else logfile,
                stderr = if (VERBOSE) "" else logfile)
  if (st != 0L) {
    if (!VERBOSE) writeLines(tail(readLines(logfile, warn = FALSE), 30))
    die("harness failed for ", tag, " -- the ", tag, " version cannot run the battery")
  }
  out
}

ref_dir <- run_side(lib_ref, "ref")
cur_dir <- run_side(lib_cur, "cur")

ref_meta <- readRDS(file.path(ref_dir, "meta.rds"))
cur_meta <- readRDS(file.path(cur_dir, "meta.rds"))
say("versions  : ", ref_meta$version, " (ref) vs ", cur_meta$version, " (under test)")
if (!identical(ref_meta$configs, cur_meta$configs))
  die("the two sides ran different config lists -- the harness is not shared")

# --- compare ----------------------------------------------------------------

quantise <- function(x) round(pmin(pmax(x, 0), 1) * 255)

compare_one <- function(name, a, b) {
  if (!identical(dim(a), dim(b)) || !identical(length(a), length(b)))
    return(list(status = "DIFF", detail = "different shape"))

  if (grepl(EXACT_RE, name)) {
    if (identical(a, b)) return(list(status = "OK", detail = "bit-identical"))
    n <- if (is.character(a)) sum(a != b) else sum(a != b, na.rm = TRUE)
    return(list(status = "DIFF",
                detail = sprintf("not bit-identical (%d of %d elements differ)", n, length(a))))
  }

  if (identical(name, "infoval")) {
    d <- abs(a - b)
    return(if (d <= INFOVAL_TOL) list(status = "OK", detail = sprintf("|d| = %.3g", d))
           else list(status = "DIFF", detail = sprintf("|d| = %.6g (tol %.0e)", d, INFOVAL_TOL)))
  }

  # NA is meaningful: generateCI(mask=) and the z-map threshold both blank out
  # pixels, so a moved NA is a changed result even where the numbers agree.
  if (!identical(is.na(a), is.na(b)))
    return(list(status = "DIFF",
                detail = sprintf("NA pattern differs (%d vs %d NA)", sum(is.na(a)), sum(is.na(b)))))

  tol <- ULPS * ULP * max(1, max(abs(a), na.rm = TRUE))
  d <- suppressWarnings(max(abs(a - b), na.rm = TRUE))
  if (!is.finite(d)) d <- 0
  ok <- d <= tol
  detail <- sprintf("max|d| = %.3g (tol %.2e)", d, tol)

  if (grepl(IMAGE_RE, name)) {
    px <- sum(quantise(a) != quantise(b), na.rm = TRUE)
    ok <- ok && px == 0L
    detail <- sprintf("%s, %d/%d 8-bit pixels differ", detail, px, length(a))
  }
  list(status = if (ok) "OK" else "DIFF", detail = detail)
}

n_ok <- 0L; unexpected <- list(); accounted <- list()

for (cfg in cur_meta$configs) {
  ref_file <- file.path(ref_dir, paste0(cfg, ".rds"))
  cur_file <- file.path(cur_dir, paste0(cfg, ".rds"))
  if (!file.exists(ref_file) || !file.exists(cur_file)) {
    unexpected[[length(unexpected) + 1L]] <-
      list(key = paste0(cfg, "/*"), detail = "config missing from one side's run")
    next
  }
  ref_res <- readRDS(ref_file); cur_res <- readRDS(cur_file)
  say("\n", cfg)

  only_cur <- setdiff(names(cur_res), names(ref_res))
  only_ref <- setdiff(names(ref_res), names(cur_res))
  for (nm in c(only_cur, only_ref))
    unexpected[[length(unexpected) + 1L]] <-
      list(key = paste0(cfg, "/", nm),
           detail = paste("output present on only one side:",
                          if (nm %in% only_cur) "under test" else REF))

  for (nm in intersect(names(cur_res), names(ref_res))) {
    key <- paste0(cfg, "/", nm)
    res <- compare_one(nm, ref_res[[nm]], cur_res[[nm]])
    exp <- expectation_for(key)
    # An expectation may also carry `check`, a predicate on the two values. A
    # key on its own accepts *any* deviation in that output, so a second,
    # unrelated regression in the same place is waved through as the deviation
    # already on file. Where the expected change has a shape that can be
    # stated, state it: the entry then only excuses the difference it describes.
    if (!is.null(exp) && is.function(exp$check) &&
        !isTRUE(exp$check(ref_res[[nm]], cur_res[[nm]]))) {
      res$detail <- paste0(res$detail, " -- and NOT the change on file for this key")
      exp <- NULL
    }
    if (identical(res$status, "OK")) {
      n_ok <- n_ok + 1L
      say(sprintf("  %-24s OK        %s", nm, res$detail))
    } else if (!is.null(exp)) {
      accounted[[length(accounted) + 1L]] <- list(key = key, detail = res$detail, exp = exp)
      say(sprintf("  %-24s EXPECTED  %s", nm, res$detail))
    } else {
      unexpected[[length(unexpected) + 1L]] <- list(key = key, detail = res$detail)
      say(sprintf("  %-24s DIFF      %s", nm, res$detail))
    }
  }
  rm(ref_res, cur_res)
  invisible(gc(verbose = FALSE))
}

# Staleness is judged per key, not per entry: an entry naming four outputs is a
# claim about all four, so three deviating does not excuse the fourth.
#
# A key for a config this run did not execute is neither fired nor stale, it is
# untested -- otherwise --quick, which skips the 512px config, would report
# every 512px expectation as stale and fail every PR.
fired_keys <- vapply(accounted, function(a) a$key, character(1))
fired_cfgs <- sub("/.*$", "", fired_keys)
all_keys <- unlist(lapply(EXPECTED, function(e) e$key))
key_ran <- vapply(all_keys, function(k) sub("/.*$", "", k) %in% cur_meta$configs, logical(1))
untested <- all_keys[!key_ran]
stale <- Filter(function(k) {
  if (endsWith(k, "/*")) !(sub("/\\*$", "", k) %in% fired_cfgs) else !(k %in% fired_keys)
}, all_keys[key_ran])

# "Deliberate" is not enough on its own: a deviation researchers are not told
# about is still a silent change to their results.
news <- if (file.exists("NEWS.md")) paste(readLines("NEWS.md", warn = FALSE), collapse = "\n") else ""
undocumented <- unique(Filter(function(h) !grepl(h, news, fixed = TRUE),
                              vapply(accounted, function(a) a$exp$news, character(1))))

say("\n", strrep("=", 70))
say(sprintf("%d checks identical within tolerance, %d expected deviations, %d unexpected",
            n_ok, length(accounted), length(unexpected)))

if (length(accounted)) {
  say("\nEXPECTED DEVIATIONS -- these change researchers' numbers, deliberately:")
  for (a in accounted) {
    say("  ", a$key, ": ", a$detail)
    say("    why : ", a$exp$reason)
    say("    news: ", a$exp$news)
  }
}

if (length(undocumented)) {
  say("\nUNDOCUMENTED DEVIATIONS -- an expectation fired, but NEWS.md does not")
  say("mention the heading it points at, so users would never hear about it:")
  for (h in undocumented) say("  missing from NEWS.md: ", h)
}

if (length(untested)) {
  say("\nNOT EXERCISED -- expectations whose configuration this run skipped:")
  for (k in untested) say("  ", k)
  say("  (expected under --quick; a release run covers them)")
}

if (length(stale)) {
  say("\nSTALE EXPECTATIONS -- listed in EXPECTED but no longer deviating:")
  for (k in stale) say("  ", k, " -- remove it from EXPECTED, or find out why it stopped firing")
}

if (length(unexpected)) {
  say("\nUNEXPECTED DIFFERENCES -- this tree does not reproduce ", REF, ":")
  for (u in unexpected) say("  ", u$key, ": ", u$detail)
  say("\nDo not release. Either fix the change, or -- if the change is intended --")
  say("document it in NEWS.md under \"Reproducibility impact\" and add it to EXPECTED")
  say("in this script with a reason. Silence is not an option here.")
}

if (length(unexpected) || length(stale) || length(undocumented)) quit(status = 1L)
if (QUICK) {
  say("\nPASS (quick): this tree reproduces ", REF, " on the configs that were run.")
  say("Run without --quick before releasing -- the 512px default config was skipped.")
} else {
  say("\nPASS: this tree reproduces ", REF, " everywhere it is supposed to.")
}
quit(status = 0L)
