# CRAN comments

## Submission type

Submission of rcicr 1.4.0 as an update to 1.3.0, published on 2026-09-02. The maintainer address `rdotsch@gmail.com` is unchanged.

## What changed

`NEWS.md` carries the full list. Its "Reproducibility impact" section describes changes to existing results, who is affected and what to do.

**The informational-value statistic is scored against the right reference distribution.** With
several base images and per-base noise parameters, base images after the first were scored
against the *first* base's reference. Classification images, z-maps, stimuli and responses were
never affected, and the default single-base case is unchanged. Where the saved parameter
matrices differ, `computeInfoVal2IFC()` now requires the `baseimage` label, so a call that
would previously have returned a figure from the wrong reference reports an error instead.

**Reference distributions are built from the noise the file saved** rather than reconstructed
by reopening the base images. Archived experiments whose base images have moved can therefore
be scored at all. Where a reference changes, the returned value is reported as superseding the
earlier one rather than changing silently.

**Read-only archived experiment files can now be scored without modification.** When an InfoVal
reference is missing or must be refreshed, `computeInfoVal2IFC()` computes it in memory if the
stimulus archive is not writable and leaves the file unchanged. Previously these calls could
fail while trying to save an otherwise successfully computed reference. Numerical results are
unchanged relative to the same calculation on a writable copy.

**Malformed trial input is rejected instead of silently changing trial assignments.**
`generateCI()` requires one participant identifier per trial when grouping, and empty, missing,
nonfinite, fractional or out-of-range stimulus identifiers now error. Valid numeric input is
unaffected.

**A classification image with no range renders as a uniform neutral image** instead of `NaN`.
Every value that changes was `NaN` before; no finite number moves.

**Stimuli are written for every base image** when the call also returns a data frame. Only the
first base image's PNGs were written before, though the call succeeded.

Also: an image's alpha channel is ignored when reading a base face; parallel workers start
under any `OutDec` and `scipen` setting, where a comma decimal separator combined with scientific notation previously prevented startup; and reference simulation no longer re-copies its stimulus matrix per
iteration.

## Test environments

Checked release-branch package sources at commit `be704a4b2ebdb5ec85e6d4795728673c9f4665e4` on 2026-09-13. The Linux R-hub artifact contains `rcicr_1.4.0.tar.gz` with SHA-256 `0ab791c88293430044aa429afd840dcdd0d727bdce41387c6efc260bce65a143`.

* GitHub Actions, R 4.6.1 (2026-06-24): Ubuntu 24.04.5 LTS (x86_64), macOS Tahoe 26.6.2 (arm64), and Windows Server 2022 (x86_64, ucrt).
* GitHub Actions, R-devel (2026-09-12 r90533): Ubuntu 24.04.5 LTS (x86_64).
* R-hub, R-devel (2026-09-12 r90533): Ubuntu 24.04.5 LTS (x86_64) and Windows Server 2022 (x86_64, ucrt).
* R-hub, R-devel (2026-09-12 r90532): macOS Sequoia 15.7.9 (x86_64).
* **Pending: win-builder R-release and R-devel for a tarball built from these package sources.**
* **Pending: a complete check including the PDF and HTML manual checks.**

## R CMD check results

All four GitHub Actions jobs report `Status: OK`: 0 errors, 0 warnings, 0 notes. These runs used `--no-manual --as-cran`. The documentation and citation-source checks also pass. [Run 34750613911](https://github.com/rdotsch/rcicr/actions/runs/34750613911).

All three R-hub artifact `00check.log` files report `Status: OK`: 0 errors, 0 warnings, 0 notes. R-hub used `--no-manual --as-cran` and did not run incoming feasibility, so these results do not establish either check. [Run 34750744174](https://github.com/rdotsch/rcicr/actions/runs/34750744174).

The repository's full reproducibility gate passes against both v1.0.1, the published baseline, and v1.3.0, the previous release. [Run 34750613924](https://github.com/rdotsch/rcicr/actions/runs/34750613924).

**Pending: both win-builder logs and their incoming-feasibility results, and a complete manual check. Replace this paragraph with the actual results before submission.** Do not carry forward the reinstatement NOTEs from 1.3.0.

## Downstream dependencies

No reverse dependencies are listed on the [CRAN package page, ETH Zurich mirror](https://stat.ethz.ch/CRAN/web/packages/rcicr/index.html), checked on 2026-09-13. This is a listing check, not a revdepcheck run.

## Notes

* `parallel`/`doSNOW` respect `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2 when it is
  set, so no example, test or vignette uses more than two cores under check.
* Ten test files call `skip_on_cran()`. All are development guards rather than checks that
  validate an installation: a golden-master regression baseline and the cross-platform z-map
  literals, a pipeline smoke test, a signal-recovery test, four serial-versus-parallel
  agreement checks, a progress-reporting check, and the slower path-handling and
  trial-validation cases. They run in CI. Where a file skips only part of itself,
  the single-core and fast cases run everywhere.
