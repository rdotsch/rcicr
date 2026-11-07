# CRAN comments

## Submission type

Update of rcicr 1.5.0, published on CRAN on 2026-09-27. The maintainer address is unchanged.

## What changed

NEWS.md describes the user-facing changes and their reproducibility impact. This release adds references matched to the stimuli and pixels of a classification image, a Gram-matrix method for computing reference norms, batch InfoVal, trial and scaling metadata, stimulus PNG checking, and more explicit validation. A reference stored by 1.4.0 or later, or drawn with a `response_seed`, is reused as stored; other stored references are rebuilt once, with a message. The NEWS.md reproducibility section explains when newly computed values or PNG pixels differ.

## Test environments and R CMD check results

The checks below used commit `670a3b41e1`.

* Local: R 4.3.3, Ubuntu 24.04.5 LTS, x86_64. `R CMD check --as-cran`: 0 errors,
  0 warnings, 2 NOTEs, PDF and HTML manuals included. Both NOTEs are local:
  this environment could not verify the current time, and R 4.3.3's
  hard-coded list of superseded packages names `doSNOW`, which win-builder did
  not report for 1.5.0.
* GitHub Actions ([run 42373720852](https://github.com/rdotsch/rcicr/actions/runs/42373720852)):
  R release on Ubuntu 24.04.5 LTS, Windows and macOS, and R-devel on Ubuntu;
  all four jobs passed; `--no-manual`.
* win-builder R-release and R-devel: pending.
* R-hub: pending.

The full reproducibility gate also passed. Against the published v1.0.1
baseline (run locally), 211 checks were identical within tolerance, 22
documented expected deviations fired, and there were 0 unexpected deviations.
Against v1.5.0 ([run 42373720808](https://github.com/rdotsch/rcicr/actions/runs/42373720808)),
255 checks were identical within tolerance, with 0 expected and 0 unexpected
deviations.

## Downstream dependencies

Pending: recheck reverse dependencies against the current CRAN index before submission.
