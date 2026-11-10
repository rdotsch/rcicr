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

To be completed from the checks on the release branch.

## Downstream dependencies

There are no reverse dependencies on CRAN: `tools::package_dependencies()`
over the current CRAN index, all dependency types, recursive, returns none.
