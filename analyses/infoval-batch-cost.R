# Per-CI cost of computeInfoVal2IFC() on a 512px stimulus file, which is what
# batchComputeInfoVal2IFC() (#85) saves by loading the file once and resolving
# each reference once. Run from the repository root with rcicr installed:
#   Rscript analyses/infoval-batch-cost.R
suppressMessages(library(rcicr))
d <- tempfile(); dir.create(d)
base <- file.path(d, "b.png"); set.seed(1); png::writePNG(matrix(runif(512^2), 512), base)
invisible(capture.output(generateStimuli2IFC(list(face = base), n_trials = 300, img_size = 512, stimulus_path = d, seed = 1, ncores = 1, save_as_png = FALSE)))
rd <- list.files(d, pattern = "Rdata$", full.names = TRUE)
cat("file MB:", file.size(rd) / 1e6, "\n")
t_ref <- system.time(invisible(capture.output(suppressMessages(generateReferenceDistribution2IFC(rd, iter = 10000, ncores = 1)))))[["elapsed"]]
set.seed(2); ci <- NULL
invisible(capture.output(ci <- generateCI(1:300, sample(c(1, -1), 300, TRUE), "face", rd, save_as_png = FALSE)))
t_hit <- system.time(for (i in 1:5) invisible(capture.output(computeInfoVal2IFC(ci, rd))))[["elapsed"]] / 5
t_load <- system.time(for (i in 1:5) { e <- new.env(); load(rd, envir = e) })[["elapsed"]] / 5
cat(sprintf("reference: %.1fs; cached-hit call: %.2fs; bare load(): %.2fs\n", t_ref, t_hit, t_load))
Sys.chmod(rd, "0444")
t_ro <- system.time(invisible(capture.output(suppressMessages(computeInfoVal2IFC(ci, rd, force_gen_ref_dist = TRUE)))))[["elapsed"]]
cat(sprintf("forced/unstorable reference per call: %.1fs\n", t_ro))

# The same 20 CIs, looped and batched, stored reference and unstorable one.
cis <- lapply(1:20, function(i) {
  set.seed(100 + i)
  out <- NULL
  invisible(capture.output(out <- generateCI(1:300, sample(c(1, -1), 300, TRUE), "face", rd, save_as_png = FALSE)))
  out
})
Sys.chmod(rd, "0644")
timed <- function(expr) system.time(invisible(capture.output(suppressMessages(expr))))[["elapsed"]]
loop_hit <- timed(for (ci in cis) computeInfoVal2IFC(ci, rd))
batch_hit <- timed(batchComputeInfoVal2IFC(cis, rd))
loop_seeded <- timed(for (ci in cis) computeInfoVal2IFC(ci, rd, response_seed = 5))
batch_seeded <- timed(batchComputeInfoVal2IFC(cis, rd, response_seed = 5))
cat(sprintf("20 CIs, stored reference: loop %.1fs, batch %.1fs (%.0fx)\n", loop_hit, batch_hit, loop_hit / batch_hit))
cat(sprintf("20 CIs, response_seed: loop %.1fs, batch %.1fs (%.0fx)\n", loop_seeded, batch_seeded, loop_seeded / batch_seeded))
