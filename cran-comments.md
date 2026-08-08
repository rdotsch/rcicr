# CRAN comments

## Response to your review of 1.2.1

Thank you for the review. This is 1.2.3, addressing each point.

1. **"Functions to" at the start of the description** — removed.
2. **References in the description** — added, as
   `Dotsch and Todorov (2012) <doi:10.1177/1948550611430272>` and
   `Brinkman, Todorov and Dotsch (2017) <doi:10.1080/10463283.2017.1381469>`. The package
   help page, regenerated for point 4 below, cites those two and a third, Dotsch, Wigboldus,
   Langner and Van Knippenberg (2008) `<doi:10.1111/j.1467-9280.2008.02186.x>`, which the
   previous page carried. All three DOIs were verified against CrossRef.
3. **`T`/`F` → `TRUE`/`FALSE`** — all 11 occurrences, including the four argument defaults
   visible in the two `.Rd` files you quoted. No value changed.
4. **Commented-out example lines** — removed. On `generateNoisePattern.Rd` we found none;
   sweeping every `.Rd` file rather than the one named located it in
   `man/rcicr-package.Rd`, whose `\examples` section was the single line
   `#simple examples will be added soon.` from 2016. That page was the package's only
   hand-maintained `.Rd`, which is why a sweep of the roxygen sources missed it. It is now
   generated and has no `\examples` section, and **every `.Rd` in the package is now
   generated.**
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
7. **`par()` in `R/plotZmap.R` and in the vignette** — both restore now, using the two forms
   your mail shows. `plotZmap()` takes the second: `oldpar <- par(mar = ...)` captures only
   the parameter it sets, restored under an immediate `on.exit()` that also closes the PNG
   device the function opens. We prefer that form in package code, because
   `par(no.readonly = TRUE)` additionally captures derived parameters such as `pin`, which a
   subsequent `plot.window()` invalidates, so restoring them errors. The vignette takes the
   first: it records `par(no.readonly = TRUE)` in its setup chunk and restores it in a final
   chunk, with an inner `par(mar = ...)` and restore around the one figure that needs
   different margins.

## Submission type

Resubmission of a package archived on 2021-06-08, "as email to the maintainer was
undeliverable" — administrative, with no policy violation or code defect involved. The
maintainer address `rdotsch@gmail.com` is current, monitored, and the address this is sent
from.

The archived CRAN version is 0.3.4.1. Versions 1.0.1 through 1.2.2 exist only as GitHub
releases made while the package was off CRAN, so nothing is missing between them and this
submission; `NEWS.md` carries their entries.

## Test environments

* local: Ubuntu 24.04, R 4.3.3 — `R CMD check --as-cran`, with
  `_R_CHECK_CRAN_INCOMING_=TRUE` and `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`
* GitHub Actions: ubuntu-latest on R release and R devel, macos-latest and windows-latest
  on R release — all green
* win-builder, R-devel (2026-08-06 r90366 ucrt) — 1 NOTE, the incoming feasibility one below
* win-builder, R-release (4.6.1) — 1 NOTE, the same one
* R-hub, R-devel on Linux (r90185), Windows (r90366 ucrt) and macOS (r90190) —
  `Status: OK` on all three, no notes

## R CMD check results

0 errors | 0 warnings | 2 notes.

**CRAN incoming feasibility** — `New submission`, `Package was archived on CRAN`, both
expected for a reinstatement.
win-builder flags four spellings in `DESCRIPTION` — `Brinkman`, `Dotsch` and `Todorov` are
the author surnames from the references you asked us to add, and `psychophysical` is the
standard adjective for the field the method comes from.

**Future file timestamps** — our build machine cannot reach the time service the check
uses. Local only; it does not reproduce on win-builder.

## Downstream dependencies

None. The package has been off CRAN since 2021.

## Notes

* `parallel`/`doParallel` respect `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2
  when it is set, so no example, test or vignette uses more than two cores under check.
* Four test files call `skip_on_cran()` — development guards (a golden-master regression
  baseline, a pipeline smoke test, a signal-recovery test, and a serial/parallel agreement
  check) that are not needed to validate an installation. They run on every push in CI.
