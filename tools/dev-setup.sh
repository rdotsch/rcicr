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

# install_deps(".") and read.dcf("DESCRIPTION", ...) below resolve relative to
# the working directory, not this script's location -- so invoking this by
# path from elsewhere would read (or fail to find) the wrong DESCRIPTION.
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# A minimal container commonly runs as root with no sudo installed at all, so
# `sudo apt-get` fails with "sudo: command not found" before anything is
# installed. Only reach for sudo when not already root.
SUDO=""
if [ "$(id -u)" -ne 0 ]; then
  SUDO="sudo"
fi

$SUDO apt-get update -qq
$SUDO apt-get install -y -qq \
  build-essential gfortran pandoc \
  libpng-dev libjpeg-dev libcurl4-openssl-dev libxml2-dev \
  libfreetype-dev libtiff-dev libharfbuzz-dev libfribidi-dev libfontconfig1-dev \
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
repos <- "https://cloud.r-project.org"

# Bootstrap remotes itself -- both install_deps() and the roxygen2 pin below
# need it.
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = repos, lib = lib)
}

# install_deps(dependencies = TRUE) reads Imports/Suggests straight from
# DESCRIPTION -- version constraints, a Remotes: field, all of it -- so
# nothing here needs to duplicate that list by hand and it cannot drift from
# what R CMD INSTALL below actually requires.
remotes::install_deps(".", dependencies = TRUE, repos = repos, upgrade = "never")

# dev tooling -- not a package dependency, but needed to work in this repo:
# lintr/pkgload for #183, devtools for the CONTRIBUTING.md "Getting set up"
# workflow
extra <- c("lintr", "pkgload", "devtools")
missing_extra <- setdiff(extra, rownames(installed.packages()))
if (length(missing_extra)) install.packages(missing_extra, repos = repos, lib = lib)

# roxygen2 unpinned would drift from DESCRIPTION'"'"'s RoxygenNote, and a
# locally-generated man/NAMESPACE diff from that drift is not a real
# documentation change -- see R-CMD-check.yaml, which pins the same way.
want <- read.dcf("DESCRIPTION", "RoxygenNote")[[1]]
have <- tryCatch(as.character(packageVersion("roxygen2")), error = function(e) NA_character_)
if (!identical(have, want)) {
  remotes::install_version("roxygen2", version = want, repos = repos, lib = lib)
}

# dev_package_deps()$diff compares the installed version against CRAN'"'"'s
# latest, not against DESCRIPTION'"'"'s own constraint -- confirmed directly:
# temporarily requiring withr (>= 99.0.0) against an installed 3.0.3 still
# reported diff = 0. So on a rerun after a `git pull` raises a Suggests
# minimum, upgrade = "never" above leaves the old version in place and a
# check based on dev_package_deps() alone would not catch it. Check
# DESCRIPTION'"'"'s declared constraints against installed versions directly,
# and upgrade only the packages that actually fail one.
cmp_ops <- list(">=" = `>=`, ">" = `>`, "<=" = `<=`, "<" = `<`, "==" = `==`, "!=" = `!=`)
constraints <- desc::desc_get_deps()
constraints <- constraints[constraints$package != "R" & constraints$version != "*", ]
unmet <- Map(function(pkg, spec) {
  op <- regmatches(spec, regexpr("^[<>=!]+", spec))
  req <- trimws(sub("^[<>=!]+", "", spec))
  have <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
  if (!is.na(have) && cmp_ops[[op]](package_version(have), package_version(req))) NA_character_ else pkg
}, constraints$package, constraints$version)
unmet <- unique(unlist(unmet, use.names = FALSE))
unmet <- unmet[!is.na(unmet)]
if (length(unmet)) install.packages(unmet, repos = repos, lib = lib)

# install.packages() only warns on a failed build, so a script that stops
# only on a nonzero exit status can silently finish "successfully" missing a
# tool -- confirmed: devtools previously failed here (missing system libs,
# now added above) and this step did not catch it. Verify explicitly, using
# dev_package_deps() and the constraint check above rather than a second
# hard-coded list.
deps <- as.data.frame(remotes::dev_package_deps(".", dependencies = TRUE))
still_unmet <- vapply(unmet, function(pkg) {
  spec <- constraints$version[constraints$package == pkg]
  op <- regmatches(spec, regexpr("^[<>=!]+", spec))
  req <- trimws(sub("^[<>=!]+", "", spec))
  have <- tryCatch(as.character(packageVersion(pkg)), error = function(e) NA_character_)
  is.na(have) || !cmp_ops[[op]](package_version(have), package_version(req))
}, logical(1))
still_missing <- c(
  deps$package[is.na(deps$installed)],
  unmet[still_unmet],
  setdiff(c(extra, "roxygen2"), rownames(installed.packages()))
)
if (length(still_missing)) {
  stop("Failed to install: ", paste(unique(still_missing), collapse = ", "))
}
'

# generateStimuli2IFC() spawns parallel workers that each `library(rcicr)`,
# so anything touching it needs the package actually installed -- see
# CONTRIBUTING.md "Getting set up".
R CMD INSTALL --no-multiarch --with-keep.source -l "$R_LIBS_USER" .
