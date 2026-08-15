#!/usr/bin/env bash
# Rebuilds a full R dev environment for rcicr on a fresh Ubuntu container that
# has only R itself preinstalled -- no compiler, no CRAN binary repo
# configured, so install.packages() fails outright until a toolchain exists,
# then fails to compile from source without libpng/libjpeg/libxml2/libcurl
# headers. Verified end-to-end from that starting point: apt packages below
# plus every Imports/Suggests package plus lintr all installed clean, and
# `R CMD INSTALL` on the package itself succeeded.
#
# Idempotent -- rerun after a `git pull` that adds a dependency; already-
# installed R packages are skipped.
set -euo pipefail

sudo apt-get update -qq
sudo apt-get install -y -qq \
  build-essential gfortran \
  libpng-dev libjpeg-dev libcurl4-openssl-dev libxml2-dev

Rscript -e '
pkgs <- c(
  # DESCRIPTION Imports
  "matlab", "png", "jpeg", "dplyr", "scales", "viridis", "doSNOW", "foreach",
  "spatstat.explore", "spatstat.geom", "tibble", "yesno",
  # DESCRIPTION Suggests
  "testthat", "covr", "withr", "knitr", "rmarkdown",
  # dev tooling -- not a package dependency, but needed to work in this repo
  "lintr", "pkgload"
)
pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs)) install.packages(pkgs, repos = "https://cloud.r-project.org")
'

# generateStimuli2IFC() spawns parallel workers that each `library(rcicr)`,
# so anything touching it needs the package actually installed -- see
# CONTRIBUTING.md "Getting set up". devtools::install() works too if
# devtools is already present; R CMD INSTALL avoids pulling it in just for
# this.
R CMD INSTALL --no-multiarch --with-keep.source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
