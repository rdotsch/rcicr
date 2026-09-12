# CRAN comments

## Submission type

Update to rcicr 1.3.0, which you published on 2026-09-02. No new submission, no archival
involved. The maintainer address `rdotsch@gmail.com` is unchanged and is the address this is
sent from.

## What changed

`NEWS.md` carries the full list, ordered largest-impact first. Five entries change numbers a
researcher may already have, and each says who is affected and what to do; they are grouped
under "Reproducibility impact" so that nobody has to infer it from a diff.

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
under any `OutDec` and `scipen` setting, where a comma decimal separator previously stopped
them starting at all; and reference simulation no longer re-copies its stimulus matrix per
iteration.

## Test environments

*To be filled from the win-builder, R-hub and GitHub Actions runs against the release branch —
see `RELEASING.md` §2. Nothing below this heading has been run for 1.4.0 yet; the figures from
the 1.3.0 submission are deliberately not carried over.*

## R CMD check results

*To be filled from the same runs.*

The incoming-feasibility NOTE from the 1.3.0 submission — `New submission` and
`Package was archived on CRAN` — should not recur for an update, but that is an expectation
rather than a result, and the actual output goes here.

## Downstream dependencies

*To be checked with `revdepcheck` or the CRAN reverse-dependency listing once 1.3.0 has been on
CRAN long enough for any to exist. The 1.3.0 submission reported none, on the grounds that the
package had been off CRAN since 2021; that reasoning no longer applies.*

## Notes

* `parallel`/`doSNOW` respect `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2 when it is
  set, so no example, test or vignette uses more than two cores under check.
* Ten test files call `skip_on_cran()`. All are development guards rather than checks that
  validate an installation: a golden-master regression baseline and the cross-platform z-map
  literals, a pipeline smoke test, a signal-recovery test, four serial-versus-parallel
  agreement checks, a progress-reporting check, and the slower path-handling and
  trial-validation cases. They run on every push in CI. Where a file skips only part of itself,
  the single-core and fast cases run everywhere.
