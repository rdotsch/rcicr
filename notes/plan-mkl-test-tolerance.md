# Plan: tests that hold under MKL, for rcicr 1.6.1

## Verified behaviour

CRAN's additional MKL check of 1.5.0 (Ripley's `tests-MKL`, R-devel r90795, Fedora 44) fails with `[ FAIL 12 | WARN 0 | SKIP 13 | PASS 1291 ]`. Every failure is an `expect_identical()` on InfoVal reference norms that compares two computations of the same quantity in one session, and every difference is in the last bits of a double. The plainest is `test-seedless-reference.R:110`: two identical `computeInfoVal2IFC(..., response_seed = 5)` calls give `62.502882847115842` and `62.502882847115885`. MKL does not promise run-to-run bitwise reproducibility by default (its results depend on memory alignment and threading), so these tests assert something the BLAS does not guarantee. The package's numbers are not wrong. The 1.6.0 pretest mail, recorded verbatim in `notes/cran-review-1.6.0-pretest.md`, says CRAN will recheck this by hand.

Four of the five failing tests are unchanged in `v1.6.0` and still use `expect_identical()` on recomputed norms: `test-generateReferenceDistribution2IFC.R` ("predating noise_type"), `test-legacy-rdata.R:29-30`, `test-reference-from-saved-noise.R:33`, `test-seedless-reference.R:8`. The Gram method that 1.6.0 makes the default does more of its work through BLAS, so 1.6.0 is expected to fail at least as much.

## Change

1. **Measure first.** A workflow, `mkl.yaml`, on an x86 Ubuntu runner: Ubuntu's R and dependencies from `tools/setup-container-r.sh`, Intel MKL from Ubuntu's `intel-mkl` made the system BLAS and LAPACK, a printed `extSoftVersion()`/`La_library()` proving R uses it, then `R CMD check --no-manual` on the built tarball, as Ripley's check runs it (`skip_on_cran()` fires). Run it on this branch before any test is edited. Its list of failures, not the four above, is the list step 2 fixes.
2. **Loosen only what MKL can move.** A comparison between two *computations* becomes `expect_equal(..., tolerance = 1e-12)`; the largest difference in the log is a relative 7e-16. A comparison of a value *stored and read back* (a reused reference, an unchanged file, its `md5sum`) stays `expect_identical()`: that holds under any BLAS, and it is what the reuse contract promises. Where one test does both, split the assertion, not the test.
3. **Keep `mkl.yaml`** as a non-required check on pull requests to `main` and on `workflow_dispatch`, so a new exact comparison is caught here rather than by CRAN. It joins no required set; `MAINTENANCE.md` lists it.
4. **Wording.** NEWS for 1.6.1 states the test change and that no number changes. Where NEWS or help says a reference is reproduced "bit for bit", qualify it as holding with the same BLAS, which is all MKL allows.
5. Version 1.6.1, `cran-comments.md` answering the MKL issue, the release gate (expect 0 unexpected against v1.0.1 and v1.6.0), and the release steps as for 1.6.0. Delete this plan before marking the PR ready.

## Main risk

That the GitHub runner does not reproduce Ripley's failures: MKL's run-to-run variation depends on CPU, thread count and allocation, and may not appear on a 4-core runner. If step 1 is green on unmodified tests, it proves nothing; force the variation with `MKL_CBWR` settings or more threads, and failing that, fix the tests named in Ripley's log plus every exact comparison of recomputed norms found by reading, and say in `cran-comments.md` that the MKL fix is unverified here.
