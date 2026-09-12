#!/usr/bin/env bash
#
# Build an R toolchain able to run the test suite and the release gate inside a
# Claude Code container. Idempotent; safe to re-run.
#
#   bash tools/setup-container-r.sh            # toolchain + rcicr installed
#   bash tools/setup-container-r.sh --no-rcicr # toolchain only
#
# Why this exists rather than the usual `install.packages()`:
#
# The container has no R, and every CRAN mirror is refused by the egress proxy
# (cloud.r-project.org, cran.r-project.org, cran.rstudio.com and
# packagemanager.posit.co all answer 403 to CONNECT). GitHub tarball hosts are
# refused too, so `remotes::install_github()` cannot stand in. Do not try to
# route around that -- it is an organization egress policy, not a fault.
#
# Ubuntu's universe pocket *is* reachable, and it carries Debian builds of
# every dependency this package and its pinned v1.0.1 gate reference need,
# except one. So: apt for everything, and a local stub for the exception.
#
# The exception is `yesno`. It has no Debian package. rcicr reaches it from a
# single `if (interactive())` branch in computeInfoVal2IFC(), guarded by a
# `ref_lookup` table whose data rows have been commented out since 2018 -- so
# every lookup misses and the prompt is unreachable in a batch run. The stub
# below therefore satisfies NAMESPACE's `import(yesno)` without pretending to
# answer: calling it fails loudly, so nothing can quietly depend on a stub.
#
# Two things that are easy to get wrong:
#   - spatstat deb names keep their dots (r-cran-spatstat.explore), unlike the
#     usual lowercase-and-dash convention. Translating the dot gives MISSING.
#   - The gate needs the v1.0.1 reference's *own* imports, which this package
#     has since dropped (raster, sp, ggplot2, plyr and friends). They are in
#     the list below so `--install-deps`, which would reach for CRAN, is never
#     needed. raster also drags in the GDAL/GEOS/PROJ system libraries.
set -euo pipefail

WITH_RCICR=1
[[ "${1:-}" == "--no-rcicr" ]] && WITH_RCICR=0

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Ubuntu's own repos resolve; third-party PPAs in this image do not, and a
# failed PPA must not fail the run.
apt-get update -qq || true

DEBIAN_FRONTEND=noninteractive apt-get install -y -qq \
  r-base-core r-base-dev \
  r-cran-matlab r-cran-png r-cran-jpeg r-cran-dplyr r-cran-scales \
  r-cran-viridis r-cran-tibble \
  r-cran-dosnow r-cran-doparallel r-cran-foreach r-cran-iterators \
  r-cran-spatstat.explore r-cran-spatstat.geom r-cran-spatstat.data \
  r-cran-spatstat.utils r-cran-spatstat.random r-cran-spatstat.sparse \
  r-cran-raster r-cran-sp r-cran-polyclip r-cran-goftest r-cran-deldir \
  r-cran-abind r-cran-tensor r-cran-dbi r-cran-assertthat r-cran-munsell \
  r-cran-plyr r-cran-ggplot2 r-cran-gridextra r-cran-purrr \
  r-cran-testthat r-cran-withr r-cran-knitr r-cran-covr

# --- the yesno stub -------------------------------------------------------
if ! Rscript -e 'quit(status = !requireNamespace("yesno", quietly = TRUE))'; then
  stub="$(mktemp -d)/yesno"
  mkdir -p "$stub/R"
  cat > "$stub/DESCRIPTION" <<'EOF'
Package: yesno
Title: Stub Standing In For The CRAN Package Of The Same Name
Version: 0.0.0.9000
Authors@R: person("rcicr", "container setup", role = c("aut", "cre"),
                  email = "noreply@example.com")
Description: Not the CRAN package. Satisfies rcicr's import(yesno) where no
    Debian build exists and CRAN is unreachable. The prompt it replaces is
    unreachable in a batch run, so this errors rather than answering.
License: MIT + file LICENSE
Encoding: UTF-8
EOF
  echo "YEAR: 2026
COPYRIGHT HOLDER: rcicr container setup" > "$stub/LICENSE"
  cat > "$stub/NAMESPACE" <<'EOF'
export(yesno)
EOF
  cat > "$stub/R/yesno.R" <<'EOF'
#' Refuse to answer an interactive prompt
#'
#' A stub. The real package asks the console a yes/no question; nothing in a
#' batch run may depend on the answer, so reaching this is a bug worth seeing.
#' @param ... Ignored.
#' @export
yesno <- function(...) {
  stop("yesno() is stubbed in this container and has no answer to give; ",
       "the caller reached an interactive prompt it should not have.",
       call. = FALSE)
}
EOF
  R CMD INSTALL --no-docs --no-help "$stub"
  rm -rf "$(dirname "$stub")"
fi

if [[ "$WITH_RCICR" == "1" ]]; then
  # generateStimuli2IFC() spawns workers that library(rcicr) themselves, so the
  # package has to be installed, not just load_all()-ed. See AGENTS.md.
  R CMD INSTALL --no-docs --no-help "$REPO"
fi

Rscript -e 'pkgs <- c("matlab","png","jpeg","dplyr","scales","viridis","tibble",
                      "doSNOW","foreach","spatstat.explore","spatstat.geom",
                      "raster","yesno","testthat","withr")
            missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
            cat("R", as.character(getRversion()), "\n")
            if (length(missing)) {
              cat("MISSING:", paste(missing, collapse = ", "), "\n"); quit(status = 1)
            }
            cat("all dependencies present\n")'
