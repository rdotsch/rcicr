# CRAN comments

## Response to your review of 1.2.1

Thank you for the review, and apologies for the delay in returning it.

**This is the resubmission answering that review**, and nothing has been sent to CRAN in
between. Every point below was addressed shortly afterwards, in versions 1.2.2 and 1.2.3, but
those were released only on GitHub and never submitted, so this is the first time the answers
reach you. The version is 1.3.0 rather than 1.2.2 because development continued in the
meantime; "Also new since 1.2.1" below lists what that added.

1. **"Functions to" at the start of the description** — removed.
2. **References in the description** — added, as
   `Dotsch and Todorov (2012) <doi:10.1177/1948550611430272>` and
   `Brinkman, Todorov and Dotsch (2017) <doi:10.1080/10463283.2017.1381469>`. The package
   help page, regenerated for point 4 below, cites those two and a third, Dotsch, Wigboldus,
   Langner and Van Knippenberg (2008) `<doi:10.1111/j.1467-9280.2008.02186.x>`, which the
   previous page carried. All three DOIs were verified against CrossRef.
3. **`T`/`F` → `TRUE`/`FALSE`** — all 11 occurrences, including the four argument defaults
   visible in the two `.Rd` files you quoted. No value changed. Neither `T` nor `F` is used
   as a variable name anywhere in the package.
4. **Commented-out example lines** — removed. Of the two files you listed,
   `generateNoiseImage.Rd` had them: its whole `\examples` section was three commented
   lines, and it now runs a real example. `generateNoisePattern.Rd` already held the live
   call `generateNoisePattern(256)` and no commented line. Sweeping every `.Rd` rather than
   only those two turned up a third, `man/rcicr-package.Rd`, whose entire `\examples`
   section was `#simple examples will be added soon.` from 2016. That page was the
   package's only hand-maintained `.Rd`, which is why a sweep of the roxygen sources missed
   it. It is now generated and has no `\examples` section, and **every `.Rd` in the package
   is now generated.**
5. **`\dontrun{}` / `\donttest{}`** — all nine wrappers removed. Every example now runs,
   at 32x32 pixels over six trials, and no single example comes close to your five-second
   guideline.
6. **Writing to a default path** — `stimulus_path`, `targetpath` and `zmaptargetpath` have
   lost their defaults (`./stimuli`, `./cis`, `./zmaps`) and are required whenever a call
   writes, so a call that would previously have written now errors naming the argument to
   supply. This is a breaking change and `NEWS.md` leads with the migration. It also
   removed one you did not see: `generateReferenceDistribution2IFC()` created an empty
   `./stimuli` on every call while writing nothing to it. Examples, tests and vignettes
   write only to `tempdir()`.
7. **`par()` in `R/plotZmap.R` and in
   `inst/doc/reverse-correlation-walkthrough.R`** — both reset now. `plotZmap()` saves the
   one parameter it sets and restores it under an immediate `on.exit()`, which also closes
   the PNG device the function opens. The vignette records `par(no.readonly = TRUE)` in its
   setup chunk and restores it in a final chunk.

## Also new since 1.2.1

Development continued on GitHub while your review went unanswered. The full list is in
`NEWS.md`; three items are worth naming here.

**A correctness fix that changes which file an output image is written to.**
`generateCI(participants = ..., save_individual_cis = TRUE)` writes one classification image
per participant. It selected each participant's trials in sorted order but took the output
*filename* from the participants' order of appearance, so whenever those two orders differed
every file was saved under another participant's identifier. The images were computed
correctly; only the names were wrong.

Sorting is lexical, so this is not an exotic case: participants labelled `p1 ... p12` in
collection order sort as `p1, p10, p11, p12, p2, ...`, and in a twelve-participant test
against the last GitHub release, eleven of the twelve files carried the wrong participant's
image. The defect dates to 2017 and is in every GitHub release since, but **it has never been
on CRAN** — the archived version, 0.3.4.1 from 2016, predates the argument entirely, so no
CRAN user has ever received a mislabelled file. We mention it because a researcher installing
from GitHub may have published such a figure; `NEWS.md` documents it under "Reproducibility
impact" with the exact conditions, a check that needs no re-run, and the mapping that
recovers existing output by renaming rather than recomputing.

**Two new arguments.** `plotZmap()` gains `pointsize` and `generateCI()` gains
`zmappointsize`, so a decorated z-map fits on an image too small for the default text size.
Both default to the graphics device's own `12`, so existing calls render exactly as before.

**One fewer dependency.** `plotZmap()` no longer uses `raster`; its `...` arguments now reach
`graphics::image()`.

The version is 1.3.0 rather than 1.2.2 because of those two new arguments and the changed
z-map plotting path, not only because of the fix above.

## Submission type

Resubmission of a package archived on 2021-06-08, "as email to the maintainer was
undeliverable" — administrative, with no policy violation or code defect involved. The
maintainer address `rdotsch@gmail.com` is current, monitored, and the address this is sent
from.

The archived CRAN version is 0.3.4.1. Versions 1.0.1 through 1.2.3 exist only as GitHub
releases made while the package was off CRAN, so nothing is missing between them and this
submission; `NEWS.md` carries their entries.

## Test environments

* win-builder, R-release 4.6.1 (2026-06-24 ucrt) — 1 NOTE, the incoming feasibility one below
* win-builder, R-devel (2026-08-15 r90413 ucrt) — 1 NOTE, the same one
* R-hub, R-devel on Linux (2026-06-21 r90185), macOS (2026-06-24 r90190) and Windows
  (2026-08-15 r90413 ucrt) — `Status: OK` on all three, no notes
* GitHub Actions: ubuntu-latest on R release and R devel, macos-latest and windows-latest
  on R release — all green

## R CMD check results

0 errors | 0 warnings | 1 note, on both win-builder runs.

R-hub reports `Status: OK` with no notes on all three platforms. It does not run the
incoming feasibility check, so that is agreement with win-builder rather than a discrepancy.

**CRAN incoming feasibility** — `New submission` and `Package was archived on CRAN`, both
expected for a reinstatement, alongside the CRAN db override recording the 2021 archival.
The same NOTE lists four possibly misspelled words in `DESCRIPTION`: `Brinkman`, `Dotsch`
and `Todorov` are the author surnames from the references you asked us to add, and
`psychophysical` is the standard adjective for the field the method comes from.

## Downstream dependencies

None. The package has been off CRAN since 2021.

## Notes

* `parallel`/`doSNOW` respect `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2
  when it is set, so no example, test or vignette uses more than two cores under check.
* Five test files call `skip_on_cran()` — development guards (a golden-master regression
  baseline, a pipeline smoke test, a signal-recovery test, a serial/parallel agreement
  check, and a progress-reporting check) that are not needed to validate an installation.
  They run on every push in CI. Only the multi-core cases are skipped in the last of these;
  its single-core cases run everywhere.
