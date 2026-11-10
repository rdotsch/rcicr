# Temporary diagnostic for #425: can MKL's run-to-run variation be made to
# appear on this runner? Runs the tests from Ripley's 1.5.0 MKL log against
# the installed package, as CRAN does (NOT_CRAN unset), and reports failures.
library(testthat)
Sys.unsetenv("NOT_CRAN")
files <- c("test-seedless-reference.R", "test-legacy-rdata.R",
           "test-reference-from-saved-noise.R",
           "test-generateReferenceDistribution2IFC.R")
reps <- as.integer(Sys.getenv("PROBE_REPS", "2"))
cfg <- paste(c("MKL_NUM_THREADS", "MKL_CBWR", "MKL_DYNAMIC"),
             Sys.getenv(c("MKL_NUM_THREADS", "MKL_CBWR", "MKL_DYNAMIC"), "-"),
             sep = "=", collapse = " ")
for (r in seq_len(reps)) for (f in files) {
  res <- as.data.frame(test_file(file.path("tests/testthat", f), reporter = "silent",
                                 package = "rcicr", load_package = "installed"))
  failed <- res$test[res$failed > 0]
  cat(sprintf("[%s] rep %d %-42s failed %d of %d%s\n", cfg, r, f, length(failed),
              nrow(res), if (length(failed)) paste0(": ", paste(failed, collapse = " | ")) else ""))
}
