# CRAN comments

## Submission type

Update of rcicr 1.5.0, published on CRAN on 2026-09-27. It replaces 1.6.0,
submitted on 2026-11-10 and cancelled at my request the same day because its
tests would also have failed the MKL additional check of 1.5.0. The maintainer
address is unchanged.

## The MKL additional issue

<https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/rcicr.out> shows 12 test
failures. Each was an `expect_identical()` comparing InfoVal reference norms
from two computations in one session; under MKL they differed in the last bits
(largest relative difference 7e-16). The tests asked for bitwise
reproducibility that MKL does not promise. In 1.6.1, comparisons of two
computations of a value that passes through BLAS or LAPACK are exact only
under R's own or the reference BLAS and LAPACK, and to a relative 1e-12 under
any other. Comparisons of stored values, files and RNG state remain exact. One
user-facing message that compared rebuilt norms exactly now uses the same
tolerance. No computed value changes.

We could not reproduce the run-to-run variation itself: `R CMD check` with
Intel MKL as BLAS and LAPACK on Ubuntu 24.04 (x86_64) returned the same bits
on every repetition, on 1 and 4 threads. The fix therefore rests on
reviewing every exact comparison in the tests. Under that MKL setup, with
`MKL_THREADING_LAYER=GNU`, the tests pass.

## What changed

NEWS.md describes the changes. Apart from the tests, 1.6.1 makes
`library(rcicr)` print a startup message if one small matrix product is wrong,
which happens with Ubuntu's `intel-mkl` when `MKL_THREADING_LAYER=GNU` is not
set, and the rebuilt-reference message above.

## Test environments and R CMD check results

The checks below used commit `8f42128`, which the release matches in package
files.

* win-builder, R 4.6.2 (2026-10-29 ucrt) and R-devel (2026-11-07 r90850
  ucrt): Status: OK on both. CRAN incoming feasibility and the PDF and HTML
  manuals passed.
* R-hub, R-devel on Linux (x86_64, 2026-11-09 r90865), Windows (x86_64 ucrt,
  2026-11-09 r90865) and macOS (x86_64, 2026-11-09 r90864)
  ([run 42864865139](https://github.com/rdotsch/rcicr/actions/runs/42864865139)):
  Status: OK on all three; `--no-manual`.
* GitHub Actions ([run 42864817788](https://github.com/rdotsch/rcicr/actions/runs/42864817788)):
  R 4.6.2 on Ubuntu 24.04.5 LTS, Windows Server 2022 and macOS Tahoe
  (aarch64), and R-devel on Ubuntu: Status: OK on all four; `--no-manual`.
* Intel MKL 2020.4 as BLAS and LAPACK, R 4.3.3 on Ubuntu 24.04 (x86_64),
  `MKL_THREADING_LAYER=GNU`
  ([run 42864817618](https://github.com/rdotsch/rcicr/actions/runs/42864817618)):
  `R CMD check` passed with no test failure; `--no-manual`.
* Local: R 4.5.2, Ubuntu 26.04 LTS, aarch64. `R CMD check --as-cran`, PDF and
  HTML manuals included: 0 errors, 0 warnings, 1 NOTE. The NOTE is local: this
  environment could not verify the current time.

The full reproducibility gate also passed
([run 42864817604](https://github.com/rdotsch/rcicr/actions/runs/42864817604)).
Against the published v1.0.1 baseline, 211 checks were identical within
tolerance, 22 documented expected deviations fired, and there were 0
unexpected deviations. Against v1.6.0, 369 checks were identical within
tolerance, with 0 expected and 0 unexpected deviations.

## Downstream dependencies

There are no reverse dependencies on CRAN: `tools::package_dependencies()`
over the current CRAN index, all dependency types, recursive, returns none.
