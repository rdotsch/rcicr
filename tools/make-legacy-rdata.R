#!/usr/bin/env Rscript

# Generates the legacy `.Rdata` fixtures in tests/testthat/fixtures/, one per
# released version, by installing that version from its own tag and letting it
# write the file. The point is that the fixtures are what old rcicr *actually*
# wrote, not our reconstruction of it -- the same reason
# tools/compare-release-output.R runs the old code rather than trusting the
# golden master. See test-legacy-rdata.R for what they are for.
#
#   Rscript tools/make-legacy-rdata.R [--versions=v1.0.1,v1.1.0] [--keep]
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

# nscales per version, and these values are load-bearing rather than convenient.
#
# v1.0.1 does not save nscales, so a current reader that re-generates stimuli
# falls back to 5 -- the default then and now. Generating this fixture at 5 makes
# that fallback *correct*, which is the situation a returning researcher is
# actually in. At any other value the null would be rebuilt on a different noise
# basis than the stimuli, giving a wrong infoVal: that is what
# generateReferenceDistribution2IFC() warns about, not something to enshrine in a
# fixture and assert as usable.
#
# v1.1.0 does save it, so its value is honoured and no fallback applies. Keeping
# it at 1 -- deliberately unlike the fallback's 5 -- means a fallback firing when
# it should not would change the basis, and be visible.
NSCALES <- c("v1.0.1" = 5L, "v1.1.0" = 1L)

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

for (tag in VERSIONS) {
  sha <- suppressWarnings(system2("git", c("rev-parse", "--verify", "--quiet",
                                           paste0(tag, "^{commit}")), stdout = TRUE))
  if (length(sha) != 1L) die("cannot resolve ", tag)

  if (!tag %in% names(NSCALES))
    die(tag, " has no entry in NSCALES -- decide what its noise basis should be, ",
        "and why, before adding it")
  nscales <- NSCALES[[tag]]

  worktree <- file.path(tmp, paste0("src-", tag))
  lib <- file.path(tmp, paste0("lib-", tag))
  outdir <- file.path(tmp, paste0("out-", tag))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  say("\n=== ", tag, " (", substr(sha, 1, 7), ") ===")
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

  logfile <- file.path(tmp, paste0("install-", tag, ".log"))
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
  script <- file.path(tmp, paste0("gen-", tag, ".R"))
  # encodeString(), not shQuote(): this text is R *source*, and on Windows a
  # path's backslashes would be read as escapes -- "C:\\Users\\..." fails to parse
  # on sequences like \\U. Forward slashes work on every platform R supports.
  rsrc <- function(path) encodeString(gsub("\\\\", "/", path), quote = '"')
  writeLines(c(
    "library(rcicr)",
    sprintf("generateStimuli2IFC(base_face_files = list(base = %s),", rsrc(base_png)),
    sprintf("                    n_trials = 6, img_size = 32, nscales = %d, seed = 1,", nscales),
    sprintf("                    stimulus_path = %s, ncores = 1, save_as_png = FALSE)", rsrc(outdir))
  ), script)
  st <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                env = paste0("R_LIBS=", paste(c(lib, .libPaths()), collapse = .Platform$path.sep)),
                stdout = file.path(tmp, paste0("gen-", tag, ".log")),
                stderr = file.path(tmp, paste0("gen-", tag, ".log")))
  if (st != 0L) {
    writeLines(utils::tail(readLines(file.path(tmp, paste0("gen-", tag, ".log")), warn = FALSE), 20))
    die("generating stimuli failed for ", tag)
  }

  produced <- list.files(outdir, pattern = "\\.Rdata$", full.names = TRUE)
  if (length(produced) != 1L) die("expected one .Rdata from ", tag, ", got ", length(produced))

  version <- unname(desc[1, "Version"])
  dest <- file.path(fixtures, paste0("legacy-rdata-", version, ".Rdata"))
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
