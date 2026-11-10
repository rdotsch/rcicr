# CRAN comments

## Submission type

Update of rcicr 1.5.0, published on CRAN on 2026-09-27. The maintainer address is unchanged.

## What changed

NEWS.md describes the user-facing changes and their reproducibility impact. This release adds references matched to the stimuli and pixels of a classification image, a Gram-matrix method for computing reference norms, batch InfoVal, trial and scaling metadata, stimulus PNG checking, and more explicit validation. A reference stored by 1.4.0 or later, or drawn with a `response_seed`, is reused as stored; other stored references are rebuilt once, with a message. The NEWS.md reproducibility section explains when newly computed values or PNG pixels differ.

## Test environments and R CMD check results

The checks below used commit `5366261493`. The release differs from it in
package files only by the date in the `NEWS.md` heading.

* win-builder, R 4.6.2 (2026-10-29 ucrt) and R-devel (2026-11-07 r90850
  ucrt): Status: OK on both. CRAN incoming feasibility and the PDF and HTML
  manuals passed.
* R-hub, R-devel on Linux (x86_64, 2026-11-08 r90859), Windows (x86_64 ucrt,
  2026-11-09 r90865) and macOS (x86_64, 2026-11-09 r90864)
  ([run 42778203633](https://github.com/rdotsch/rcicr/actions/runs/42778203633)):
  Status: OK on all three; `--no-manual`.
* GitHub Actions ([run 42778195324](https://github.com/rdotsch/rcicr/actions/runs/42778195324)):
  R 4.6.2 on Ubuntu 24.04.5 LTS, Windows Server 2022 and macOS Tahoe
  (aarch64), and R-devel on Ubuntu: Status: OK on all four; `--no-manual`.
* Local: R 4.5.2, Ubuntu 26.04 LTS, aarch64. `R CMD check --as-cran`, PDF and
  HTML manuals included: 0 errors, 0 warnings, 1 NOTE. The NOTE is local: this
  environment could not verify the current time.

The full reproducibility gate also passed
([run 42778195278](https://github.com/rdotsch/rcicr/actions/runs/42778195278)).
Against the published v1.0.1 baseline, 211 checks were identical within
tolerance, 22 documented expected deviations fired, and there were 0
unexpected deviations. Against v1.5.0, 255 checks were identical within
tolerance, with 0 expected and 0 unexpected deviations.

## Downstream dependencies

There are no reverse dependencies on CRAN: `tools::package_dependencies()`
over the current CRAN index, all dependency types, recursive, returns none.
