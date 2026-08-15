#!/usr/bin/env bash
# Rebuilds a full R dev environment for rcicr on a fresh Ubuntu container that
# has only R itself preinstalled -- no compiler, no CRAN binary repo
# configured, so install.packages() fails outright until a toolchain exists,
# then fails to compile from source without libpng/libjpeg/libxml2/libcurl
# headers (and, for devtools' pkgdown/ragg/textshaping/systemfonts/gert
# chain, freetype/tiff/harfbuzz/fribidi/fontconfig/libgit2). Verified
# end-to-end from that starting point: apt packages below plus every
# Imports/Suggests package plus lintr/devtools/roxygen2 all installed
# clean, and `R CMD INSTALL` on the package itself succeeded.
#
# Idempotent -- rerun after a `git pull` that adds a dependency; already-
# installed R packages are skipped.
set -euo pipefail

sudo apt-get update -qq
sudo apt-get install -y -qq \
  build-essential gfortran pandoc \
  libpng-dev libjpeg-dev libcurl4-openssl-dev libxml2-dev \
  libfreetype-dev libtiff-dev libharfbuzz-dev libfribidi-dev libfontconfig-dev \
  libgit2-dev

# R's default site-library is root-owned on a plain apt install, and a
# non-interactive Rscript can neither write there nor prompt to create a
# personal one -- install.packages() aborts. Create and use a personal
# library explicitly rather than relying on ambient permissions.
R_LIBS_USER="$(Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))')"
mkdir -p "$R_LIBS_USER"

# pandoc: both vignettes use rmarkdown::html_vignette, so devtools::check()
# cannot build them without it.
Rscript -e '
lib <- Sys.getenv("R_LIBS_USER")

# Read Imports/Suggests from DESCRIPTION itself rather than hard-coding them
# here -- a hard-coded copy drifts the moment a dependency is added, and
# CONTRIBUTING.md already tells contributors to rerun this script after a
# `git pull` that changes one, which only helps if it actually notices.
desc_field <- function(field) {
  raw <- read.dcf("DESCRIPTION", field)[1, field]
  if (is.na(raw)) return(character(0))
  pkgs <- trimws(strsplit(raw, ",")[[1]])
  sub("\\s*\\(.*\\)$", "", pkgs)
}
pkgs <- c(
  desc_field("Imports"),
  desc_field("Suggests"),
  # dev tooling -- not a package dependency, but needed to work in this repo:
  # lintr/pkgload for #183, devtools/remotes for the CONTRIBUTING.md
  # "Getting set up" workflow
  "lintr", "pkgload", "devtools", "remotes"
)
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org", lib = lib)

# roxygen2 unpinned would drift from DESCRIPTION'"'"'s RoxygenNote, and a
# locally-generated man/NAMESPACE diff from that drift is not a real
# documentation change -- see R-CMD-check.yaml, which pins the same way.
want <- read.dcf("DESCRIPTION", "RoxygenNote")[[1]]
have <- tryCatch(as.character(packageVersion("roxygen2")), error = function(e) NA_character_)
if (!identical(have, want)) {
  remotes::install_version("roxygen2", version = want, repos = "https://cloud.r-project.org", lib = lib)
}

# install.packages() only warns on a failed build, so a script that stops
# only on a nonzero exit status can silently finish "successfully" missing a
# tool -- confirmed: devtools previously failed here (missing system libs,
# now added above) and this step did not catch it. Verify explicitly.
still_missing <- setdiff(c(pkgs, "roxygen2"), rownames(installed.packages()))
if (length(still_missing)) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "))
}
'

# generateStimuli2IFC() spawns parallel workers that each `library(rcicr)`,
# so anything touching it needs the package actually installed -- see
# CONTRIBUTING.md "Getting set up".
R CMD INSTALL --no-multiarch --with-keep.source -l "$R_LIBS_USER" \
  "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
