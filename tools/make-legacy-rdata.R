#!/usr/bin/env Rscript

# Generates the legacy `.Rdata` fixtures in tests/testthat/fixtures/, one per row
# of the FIXTURES table below, by installing that version from its own tag and
# letting it write the file. A version may appear in more than one row: v1.0.1
# contributes both a sinusoidal and a gabor fixture. The point is that the fixtures are what old rcicr *actually*
# wrote, not our reconstruction of it -- the same reason
# tools/compare-release-output.R runs the old code rather than trusting the
# golden master. See test-legacy-rdata.R for what they are for.
#
#   Rscript tools/make-legacy-rdata.R [--versions=v1.0.1,v1.1.0] [--keep]
#
# --versions takes tags or row ids ("v1.0.1-gabor"), and rewrites every fixture
# it selects.
#
# Regenerating is a deliberate act, and normally unnecessary: these files are
# frozen artefacts of releases that will not change. Doing it to make a red test
# pass is how a break in the append-only `.Rdata` contract reaches a user.
#
# Old versions need their own dependencies -- v1.0.1 imports raster, dropped in
# #186 -- and install them into the throwaway library rather than the user's.
# Every version here builds a parallel cluster whose workers call
# library(rcicr), so the package has to be installed, which is why this cannot
# just source the old R files.
#
# The fixtures carry `base_face_files` and `stimulus_path` as *absolute paths on
# this machine*, which is what rcicr has always written and cannot be fixed here.
# A reader that re-generates stimuli therefore has to repoint them; the test does
# that, and CI caught it the one time it did not -- the paths existed only on the
# machine that generated the fixtures.

args <- commandArgs(trailingOnly = TRUE)
opt <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}
VERSIONS <- strsplit(opt("versions", "v1.0.1,v1.1.0"), ",")[[1]]
KEEP <- "--keep" %in% args

# One row per fixture. The noise-basis values are load-bearing rather than
# convenient.
#
# v1.0.1 saves neither nscales nor sigma, so a current reader that re-generates
# stimuli falls back to 5 and 25 -- the defaults then and now. Generating those
# fixtures at 5 and 25 makes the fallbacks *correct*, which is the situation a
# returning researcher is actually in. At any other value the null would be
# rebuilt on a different noise basis than the stimuli, giving a wrong infoVal:
# that is what generateReferenceDistribution2IFC() warns about, not something to
# enshrine in a fixture and assert as usable.
#
# v1.1.0 does save nscales, so its value is honoured and no fallback applies.
# Keeping it at 1 -- deliberately unlike the fallback's 5 -- means a fallback
# firing when it should not would change the basis, and be visible.
#
# The gabor row exists because sigma only reaches the basis through
# generateGabor(): with sinusoidal noise the reference norms are identical at
# sigma 25 and 10, so no sinusoidal fixture can exercise the missing-sigma
# fallback at all. Measured at these settings, gabor norms move from
# 0.681/0.689/0.680 at sigma 25 to 0.615/0.620/0.626 at sigma 10, and a 1.0.1
# gabor file is the case that hits both fallbacks at once.
FIXTURES <- list(
  list(tag = "v1.0.1", noise_type = "sinusoid", nscales = 5L, sigma = 25),
  list(tag = "v1.1.0", noise_type = "sinusoid", nscales = 1L, sigma = 25),
  list(tag = "v1.0.1", noise_type = "gabor",    nscales = 5L, sigma = 25)
)

say <- function(...) cat(..., "\n", sep = "")
die <- function(...) { cat("ERROR: ", ..., "\n", sep = ""); quit(status = 2L) }

repo <- normalizePath(".")
if (!file.exists(file.path(repo, "DESCRIPTION"))) die("run this from the package root")
fixtures <- file.path(repo, "tests", "testthat", "fixtures")

tmp <- file.path(dirname(tempdir()), paste0("rcicr-legacy-", Sys.getpid()))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit({
  if (KEEP) say("\nkept: ", tmp) else unlink(tmp, recursive = TRUE)
}, add = TRUE)

# The base image is synthetic and deterministic, never a real photograph, and
# identical to the one the tests would build -- so a fixture differs from a
# freshly generated file only by the version that wrote it.
base_png <- file.path(tmp, "base.png")
set.seed(1)
png::writePNG(matrix(runif(32 * 32), 32, 32), base_png)

# Selectable by tag or by row id ("v1.0.1-gabor"). The row id matters because a
# tag now covers more than one fixture, and every run rewrites what it selects:
# asking for v1.0.1 to add the gabor file would also rewrite the sinusoidal one,
# whose embedded absolute paths differ per run. That is a frozen artefact, so
# name the row when only one of them is meant.
row_tags <- vapply(FIXTURES, `[[`, "", "tag")
row_ids <- paste(row_tags, vapply(FIXTURES, `[[`, "", "noise_type"), sep = "-")
selected <- FIXTURES[row_tags %in% VERSIONS | row_ids %in% VERSIONS]
unknown <- setdiff(VERSIONS, c(row_tags, row_ids))
if (length(unknown))
  die(paste(unknown, collapse = ", "), " has no row in FIXTURES -- decide what its ",
      "noise basis should be, and why, before adding it")

for (fixture in selected) {
  tag <- fixture$tag
  noise_type <- fixture$noise_type
  nscales <- fixture$nscales
  sigma <- fixture$sigma

  sha <- suppressWarnings(system2("git", c("rev-parse", "--verify", "--quiet",
                                           paste0(tag, "^{commit}")), stdout = TRUE))
  if (length(sha) != 1L) die("cannot resolve ", tag)

  # Two rows share v1.0.1, so every path is keyed on the row rather than the tag.
  id <- paste(tag, noise_type, sep = "-")
  worktree <- file.path(tmp, paste0("src-", id))
  lib <- file.path(tmp, paste0("lib-", id))
  outdir <- file.path(tmp, paste0("out-", id))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  say("\n=== ", id, " (", substr(sha, 1, 7), ", nscales ", nscales, ", sigma ", sigma, ") ===")
  if (system2("git", c("worktree", "add", "--detach", shQuote(worktree), tag),
              stdout = FALSE, stderr = FALSE) != 0L) die("git worktree add failed for ", tag)

  desc <- read.dcf(file.path(worktree, "DESCRIPTION"))
  deps <- unlist(strsplit(paste(
    desc[1, intersect(c("Depends", "Imports"), colnames(desc))], collapse = ","), ","))
  deps <- trimws(sub("\\(.*", "", deps))
  missing <- setdiff(deps[nzchar(deps)], c("R", rownames(utils::installed.packages())))
  if (length(missing)) {
    say("installing ", tag, "-only dependencies: ", paste(missing, collapse = ", "))
    repos <- getOption("repos")
    if (!length(repos) || any(repos == "@CRAN@")) repos <- c(CRAN = "https://cloud.r-project.org")
    utils::install.packages(missing, lib = lib, repos = repos, quiet = TRUE)
  }

  logfile <- file.path(tmp, paste0("install-", id, ".log"))
  st <- system2(file.path(R.home("bin"), "R"),
                c("CMD", "INSTALL", "--no-docs", "--no-help", "--no-byte-compile",
                  "-l", shQuote(lib), shQuote(worktree)),
                env = paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep)),
                stdout = logfile, stderr = logfile)
  if (st != 0L) {
    writeLines(utils::tail(readLines(logfile, warn = FALSE), 20))
    die("R CMD INSTALL failed for ", tag)
  }

  # Run the old version out of process, so its namespace never meets this one's.
  # Small on purpose: 32px and 6 trials keep the fixtures well under pre-commit's
  # 1 MB ceiling -- tens of KB, or ~200 KB where nscales is 5 -- and nothing here
  # depends on the stimuli looking like anything.
  script <- file.path(tmp, paste0("gen-", id, ".R"))
  # encodeString(), not shQuote(): this text is R *source*, and on Windows a
  # path's backslashes would be read as escapes -- "C:\\Users\\..." fails to parse
  # on sequences like \\U. Forward slashes work on every platform R supports.
  rsrc <- function(path) encodeString(gsub("\\\\", "/", path), quote = '"')
  writeLines(c(
    "library(rcicr)",
    sprintf("generateStimuli2IFC(base_face_files = list(base = %s),", rsrc(base_png)),
    sprintf("                    n_trials = 6, img_size = 32, nscales = %d, seed = 1,", nscales),
    sprintf("                    noise_type = %s, sigma = %s,",
            encodeString(noise_type, quote = '"'), format(sigma)),
    sprintf("                    stimulus_path = %s, ncores = 1, save_as_png = FALSE)", rsrc(outdir))
  ), script)
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                env = paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep)),
                stdout = file.path(tmp, paste0("gen-", id, ".log")),
                stderr = file.path(tmp, paste0("gen-", id, ".log")))
  if (st != 0L) {
    writeLines(utils::tail(readLines(file.path(tmp, paste0("gen-", id, ".log")), warn = FALSE), 20))
    die("generating stimuli failed for ", id)
  }

  produced <- list.files(outdir, pattern = "\\.Rdata$", full.names = TRUE)
  if (length(produced) != 1L) die("expected one .Rdata from ", id, ", got ", length(produced))

  # Only non-sinusoidal rows are suffixed, so the filenames of the two fixtures
  # that predate this table are unchanged.
  version <- unname(desc[1, "Version"])
  suffix <- if (noise_type == "sinusoid") "" else paste0("-", noise_type)
  dest <- file.path(fixtures, paste0("legacy-rdata-", version, suffix, ".Rdata"))
  if (!file.copy(produced, dest, overwrite = TRUE)) die("could not write ", dest)

  e <- new.env()
  load(dest, envir = e)
  say("wrote ", basename(dest), " (", round(file.size(dest) / 1024, 1), " KB)")
  say("  fields: ", paste(sort(ls(e)), collapse = ", "))
  say("  generator_version: ",
      if (exists("generator_version", e)) paste(format(get("generator_version", e)), collapse = "") else "(absent)",
      "   p$generator_version: ",
      if (!is.null(get("p", e)$generator_version)) paste(format(get("p", e)$generator_version), collapse = "") else "(absent)")

  system2("git", c("worktree", "remove", "--force", shQuote(worktree)),
          stdout = FALSE, stderr = FALSE)
}

say("\ndone. Commit the fixtures if this was a deliberate regeneration.")
