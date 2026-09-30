# Times generateReferenceDistribution2IFC() itself, at the package's default
# ncores unless told otherwise, in an older build and in the installed one.
# Backs the performance figures in NEWS.md for #354; the ratios, not the
# seconds, are what carry to another machine.
#
# Usage, from the repository root:
#   R CMD INSTALL -l <old-lib> <checkout of the older build>
#   R CMD INSTALL .
#   Rscript analyses/gram-reference-benchmark.R <old-lib>
#
# Each timing runs in a fresh R process, so the two builds never share a
# session. The stimulus sets are generated once, with the installed build.

old_lib <- commandArgs(TRUE)[1]
if (is.na(old_lib) || !dir.exists(old_lib)) {
  stop("Pass the library the older build is installed in.")
}

configs <- list(
  list(img_size = 128, n_trials = 300, ncores = NA),
  list(img_size = 256, n_trials = 300, ncores = 2)
)

scratch <- tempfile("benchmark")
dir.create(scratch)

stimulus_file <- function(img_size, n_trials) {
  path <- file.path(scratch, sprintf("stimuli-%d", img_size))
  dir.create(path)
  base <- file.path(path, "base.png")
  set.seed(1)
  png::writePNG(matrix(runif(img_size^2), img_size, img_size), base)
  invisible(utils::capture.output(rcicr::generateStimuli2IFC(
    list(base = base), n_trials = n_trials, img_size = img_size, stimulus_path = path, seed = 7,
    ncores = 1, save_as_png = FALSE)))
  list.files(path, pattern = "\\.Rdata$", full.names = TRUE)[1]
}

time_reference <- function(lib, rdata, ncores) {
  code <- sprintf(paste0(
    "if (nzchar('%s')) .libPaths(c('%s', .libPaths()));",
    "ncores <- %s; if (is.na(ncores)) ncores <- rcicr:::default_ncores();",
    "t <- system.time(invisible(utils::capture.output(suppressWarnings(",
    "rcicr::generateReferenceDistribution2IFC('%s', iter = 10000, ncores = ncores,",
    "response_seed = 1, save_rdata = FALSE)))))[['elapsed']];",
    "cat(t, ncores)"), lib, lib, ncores, rdata)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(code)), stdout = TRUE,
                 stderr = FALSE)
  as.numeric(strsplit(utils::tail(out, 1), " ")[[1]])
}

for (cfg in configs) {
  rdata <- stimulus_file(cfg$img_size, cfg$n_trials)
  old <- time_reference(old_lib, rdata, cfg$ncores)
  new <- time_reference("", rdata, cfg$ncores)
  cat(sprintf("%dpx, %d trials, ncores = %d: %.1fx faster (%.1f s to %.1f s)\n",
              cfg$img_size, cfg$n_trials, new[2], old[1] / new[1], old[1], new[1]))
}
