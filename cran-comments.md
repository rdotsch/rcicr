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

Checked release-branch commit `375bbcedd90ac3a5b6af1cd6e62610056056aa08` on 2026-09-12.

**Superseded: the branch now carries the read-only reference fix (#324).** Every result below, the win-builder submissions included, describes a tarball built before it. Rebuild from the current head and redo all of them before submission.

* GitHub Actions, R 4.6.1: Ubuntu 24.04.5 LTS (x86_64), macOS Tahoe 26.6.2 (arm64), and Windows Server 2022 (x86_64, ucrt).
* GitHub Actions, R-devel (2026-09-11 r90528): Ubuntu 24.04.5 LTS (x86_64).
* R-hub, R-devel (2026-09-11 r90528): Ubuntu 24.04.5 LTS (x86_64) and Windows Server 2022 (x86_64, ucrt).
* R-hub, R-devel (2026-09-12 r90532): macOS Sequoia 15.7.9 (x86_64).
* **Pending: win-builder R-release and R-devel.** The 2026-09-13 submissions of `rcicr_1.4.0.tar.gz` (SHA-256 `168aea375d3d67a8d26dd37d159dc66f99b4c4dc9ae0dd8152c383e0089e8b02`) both returned `Status: OK` — R-devel (2026-09-12 r90533 ucrt) and R 4.6.1 (2026-06-24 ucrt) — for the superseded sources above. Resubmit the rebuilt tarball and record its R versions and results.

## R CMD check results

All four GitHub Actions jobs report `Status: OK`: 0 errors, 0 warnings, 0 notes. These runs used `--no-manual --as-cran`. The documentation and citation-source checks also pass. [Run 34700615293](https://github.com/rdotsch/rcicr/actions/runs/34700615293).

All three R-hub artifact `00check.log` files report `Status: OK`: 0 errors, 0 warnings, 0 notes. R-hub used `--no-manual --as-cran` and did not run incoming feasibility, so these results do not establish either check. [Run 34718965345](https://github.com/rdotsch/rcicr/actions/runs/34718965345).

**Pending: both win-builder logs, their incoming-feasibility results, and a complete manual check. Replace this paragraph with the actual results before submission.** Do not carry forward the reinstatement NOTEs from 1.3.0.

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
