# Plan: tests that hold under MKL, for rcicr 1.6.1

## Verified behaviour

**CRAN's failure.** The additional MKL check of 1.5.0 (Ripley's `tests-MKL`, R-devel r90795, Fedora 44) fails `[ FAIL 12 | WARN 0 | SKIP 13 | PASS 1291 ]`. Every failure is an `expect_identical()` on InfoVal reference norms comparing two computations of the same quantity in one session, and every difference is in the last bits: two identical `computeInfoVal2IFC(..., response_seed = 5)` calls give `62.502882847115842` and `62.502882847115885` (`test-seedless-reference.R:110`). MKL does not promise run-to-run bitwise reproducibility by default; the tests assert what the BLAS does not guarantee. Four of the five failing tests are unchanged in `v1.6.0` (`test-generateReferenceDistribution2IFC.R` "predating noise_type", `test-legacy-rdata.R:29-30`, `test-reference-from-saved-noise.R:33`, `test-seedless-reference.R:8`), so 1.6.0 fails too; the pretest mail is in `notes/cran-review-1.6.0-pretest.md`.

**The `mkl` job on this branch** (x86 `ubuntu-24.04`, Ubuntu's R 4.3.3 and `intel-mkl` 2020.4 selected through update-alternatives, R confirmed on `libmkl_rt.so`):

- Without `MKL_THREADING_LAYER=GNU`, MKL's Intel OpenMP and R's GNU OpenMP in one process give wrong products from base R: `notes/mkl-dgemm-probe.R` found `a %*% b` wrong in almost every column (median ratio of `r'Ar` exactly 2) and `tcrossprod()` wrong, and the tests failed 36 to 40 times with errors up to 45%. That is a broken BLAS setup, not rcicr: the same tests pass under multithreaded OpenBLAS (checked locally, aarch64).
- With `MKL_THREADING_LAYER=GNU`, the probe is clean (0 wrong columns in 20 configurations; `tcrossprod` within 1.1e-15) and `R CMD check` gives `Status: OK` with no test failure (run 42809113898, attempt 2).

So the runner does not reproduce CRAN's failure: on 4 cores MKL returned the same bits on every repetition. A green `mkl` job cannot show that a fix works.

## Change

1. **Fix by reading, not by the job.** The five tests in Ripley's log, then every assertion in `tests/testthat/` that compares two *recomputations* of a BLAS-dependent value with `expect_identical()`: reference norms, InfoVal, classification images, z-maps. A grep for `expect_identical` on such names finds about 80 to sort. Each becomes `expect_equal(..., tolerance = 1e-12)`; the largest difference in Ripley's log is a relative 7e-16.
2. **Keep exact what is exact.** A value stored and read back (a reused reference, `.Rdata` fields, a file's `md5sum`, RNG state, metadata) stays `expect_identical()`; that holds under any BLAS and is what the reuse contract promises. Where one test does both, split the assertion, not the test. No assertion is deleted or loosened beyond 1e-12.
3. **Try to make the variation appear** on the runner before relying on reading alone: rerun the four 1.6.0 tests from Ripley's log in a loop under `MKL_NUM_THREADS` of 1, 2 and 4 and with `MKL_CBWR` unset versus `AUTO`, and compare results across repetitions. If bits move, the job verifies step 1 against the unmodified tests; if not, `cran-comments.md` says the MKL fix is by inspection and could not be reproduced here.
4. **Keep `mkl.yaml`** as a non-required check on pull requests to `main` and on `workflow_dispatch`, with `MKL_THREADING_LAYER=GNU` and a 15-minute limit on the setup step (one attempt hung for 47 minutes on an Ubuntu mirror). It catches outright MKL failures, not run-to-run variation, and says so. `MAINTENANCE.md` lists it. `notes/mkl-dgemm-probe.R` and its step are removed.
5. **Wording.** NEWS for 1.6.1: the test change, that no computed number changes, and that MKL users of R from Debian or Ubuntu need `MKL_THREADING_LAYER=GNU` (base R's products are wrong without it, so rcicr's are too). Where NEWS or help promises a reference "bit for bit", qualify it as holding with the same BLAS.
6. Version 1.6.1, `cran-comments.md` answering the MKL issue, the release gate (expect 0 unexpected against v1.0.1 and v1.6.0), and the release steps as for 1.6.0. Delete this plan before marking the PR ready.

## Inventory and three amendments

**What can move.** In `R/`, BLAS is reached only in `R/reference-base.R` (the reference norms, both methods) and LAPACK through `norm(x, "f")` (`reference-base.R:220`, `reference-stimuli.R:78`): reference norms and every InfoVal. Classification images, z-maps, stimuli and correlations do not touch either.

**The 40 assertions to change**, each an exact comparison of two computations of norms or an InfoVal (an InfoVal always recomputes the CI's norm, so it is never a stored value):

- `test-generateReferenceDistribution2IFC.R:96`
- `test-legacy-rdata.R:145, 146, 191, 192, 220, 221`
- `test-seedless-reference.R:98, 110, 123, 129`
- `test-reference-from-saved-noise.R:342, 442, 454, 815, 820, 912, 990`
- `test-batchComputeInfoVal2IFC.R:58, 70, 81, 101, 170`
- `test-reference-checks.R:40`
- `test-gram-reference.R:102, 119, 132`
- `test-reference-stimuli.R:63, 75, 165`
- `test-rng-kind.R:41, 42, 43, 59, 60`
- `test-masked-infoval.R:68, 87, 96, 128, 129`

Every other `expect_identical()` in the suite compares a stored value, a file, RNG state, a message, a count or a CI, and stays. `expect_false(identical(...))` assertions stay (BLAS noise only widens a difference), except `test-generateReferenceDistribution2IFC.R:297`, which asserts that two routes differ by rounding; it becomes a check of the saved `reference_norms_method`.

**Amendment 1: a helper, not a blanket tolerance.** `expect_same_arithmetic(actual, expected)` in `helper-fixtures.R` is `expect_identical()` when R runs its reference BLAS, and `expect_equal(tolerance = 1e-12)` otherwise, decided from `extSoftVersion()[["BLAS"]]`; anything it does not recognise as the reference BLAS gets the tolerance. CI and every machine with the reference BLAS keep checking the "images reproduces the arithmetic bit for bit" promise exactly; only a BLAS that cannot keep it relaxes.

**Amendment 2: one line of `R/`.** `R/reference-base.R:75` decides whether to tell the user that a rebuilt reference "gave different values" and supersedes earlier InfoVals with `!identical(norms, entry$norms)`. Under MKL a rebuild differing in the last bits gets that message, which is untrue, and `test-reference-from-saved-noise.R:342` ("a refresh that changes nothing is silent") fails on it under any test tolerance. It becomes `!isTRUE(all.equal(norms, entry$norms, tolerance = 1e-12))`. Only the message changes; no returned or stored value does. It goes in NEWS.

**Amendment 3: a product check when rcicr is attached (maintainer's request).** Under the broken threading setup above, rcicr returns reference norms and InfoVals that are wrong by tens of percent, with no warning; NEWS alone reaches few of those users. The probe's numbers (in the `e0a0c96` run) say where a check can see it: a 40 x 40 symmetric matrix times a 40 x 20 matrix of ±1 had 3993 of 4000 columns wrong over 200 repetitions, and with a matrix of 20 x 20 or larger every case was wrong in at least 93.7% of columns; only 8 x 8 by 3 or 5 came out right. A correct BLAS differs from exact arithmetic by about 1e-15.

- **The check.** An internal `blas_product_ok()` multiplies one fixed 40 x 40 symmetric matrix by one fixed 40 x 20 matrix of ±1 with `%*%`, and compares it, at a relative 1e-8, with the same product in plain R arithmetic (`colSums(a * b[, j])`, which no BLAS reaches). It changes no option and no RNG state, and costs microseconds.
- **Where.** `.onAttach()` in `R/zzz.R`, reporting through `packageStartupMessage()`, so `suppressPackageStartupMessages()` silences it and CRAN's rule against messages from `.onLoad()` holds. Only a wrong product prints anything. The message says that this R's matrix product was wrong in a check, that reference norms and InfoVals computed in the session will be wrong, and that with Intel MKL from Debian or Ubuntu, `MKL_THREADING_LAYER=GNU` set before R starts fixed it in our tests.
- **What it does not do.** It does not stop or alter any computation, and it does not run for `rcicr::` calls without `library(rcicr)`. Parallel workers attach rcicr too; a message there goes to the worker's discarded output.
- **Tests.** `blas_product_ok()` is `TRUE` here; with it mocked to `FALSE`, `.onAttach()` gives the message, and with it `TRUE`, none.
- **Measured on the real setup, once.** A temporary step in the `mkl` job attaches rcicr with `MKL_THREADING_LAYER` unset and requires the message, then attaches it with `GNU` and requires none. The step is removed before the PR is ready; the run goes in the squash message.
- NEWS gets it under "Bug fixes", with the threading advice it now carries.

## Main risk

Missing an exact comparison that MKL can move, and CRAN failing 1.6.1 for it. Reading cannot prove completeness. Step 3 is the only way to measure it here; if it does not reproduce the variation, the guard is the sorting rule itself: every `expect_identical()` left on a numeric value must be one that is stored and read back, and the PR lists each that stays, with why.
