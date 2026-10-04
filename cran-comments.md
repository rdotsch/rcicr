# CRAN comments

## Submission type

Update of rcicr 1.5.0, published on CRAN on 2026-09-27. The maintainer address is unchanged.

## What changed

NEWS.md describes the user-facing changes and their reproducibility impact. This release adds references matched to the stimuli and pixels of a classification image, a Gram-matrix method for computing reference norms, batch InfoVal, trial and scaling metadata, stimulus PNG checking, and more explicit validation. Existing stored references remain reusable. The NEWS.md reproducibility section explains when newly computed values or PNG pixels differ.

## Test environments and R CMD check results

Pending release-branch checks. Record the exact release commit, R versions, platforms, errors, warnings and NOTEs here after the checks complete:

- Local R CMD build and R CMD check --as-cran, including PDF and HTML manuals: pending.
- GitHub Actions R CMD check and the full reproducibility gate against v1.0.1 and v1.5.0: pending.
- win-builder R-release and R-devel: pending.
- R-hub: pending.

## Downstream dependencies

Recheck reverse dependencies against the current CRAN index before submission.

## Notes

The release review is tracking the historical-PNG checker finding in #419. Resolve or document its effect before release. Do not carry the check results from 1.5.0 into this submission.
