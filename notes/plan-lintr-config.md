# Plan: add lintr as a config plus a CI job (#183)

## What was verified, and how

This sandbox initially could not install `lintr` at all — no C compiler, then
missing `libxml2`/`libcurl` headers. Fixed with `apt-get install
build-essential libpng-dev libjpeg-dev libcurl4-openssl-dev gfortran
libxml2-dev`, then `R CMD INSTALL .` after installing every `Imports`/`Suggests`
package from source. lintr 3.4.0 now loads and `lint_package()` runs against
the installed package (needed for `object_usage_linter`, per
`?lintr::executing_linters`).

`lintr::lint_package()` on bare `linters_with_defaults()` returns **1266**
lints. Before writing an exclusion baseline for that number, I checked each
large category against `CONTRIBUTING.md` → "Code conventions" instead of
assuming the defaults are what this project wants:

| linter | default hits | CONTRIBUTING.md says | action |
|---|---|---|---|
| `line_length_linter` | 612 | "Line length is not enforced; 35 lines in R/ already exceed 100 characters" | disable |
| `return_linter` | 18 | "explicit `return()` at the end of exported functions" — lintr's default wants the opposite (implicit return) | disable |
| `quotes_linter` | 263 (double→single direction) | "single quotes... not worth churn to unify, but write new code single-quoted" | disable — the maintainer already declined to enforce this one way; flipping the linter's default direction to single actually makes it worse (**922** hits, because `tests/` idiomatically uses double-quoted `expect_*()` arguments) |
| `object_name_linter` | 51 | "Exported functions are camelCase... frozen" / "Arguments and list fields are snake_case... also frozen" | reconfigure `styles = c("camelCase", "snake_case", "symbols")` — drops to 4 genuine hits |
| *(not in defaults)* `T_and_F_symbol_linter` | — | "`TRUE`/`FALSE`, never `T`/`F`... All occurrences were replaced in 1.2.2" | **enable** — nothing currently guards against reintroducing it, and it's exactly the class of bug the issue names |

Everything else is left at its default. With these five changes, the baseline
drops to **326** lints across **71** (file, linter) pairs — still real, still
worth a to-fix list, but not dominated by fights against decisions already
made. Command used to check, reusable for re-verification:

```r
options(lintr.linter_file = "<candidate .lintr>")
lints <- lintr::lint_package("."); length(lints)
```

One more thing worth surfacing to the reviewer rather than silently fixing:
`assignment_linter` still finds **9** `=`-assignments in `R/`
(`generateGabor.R:19`, `generateNoisePattern.R:{21,24,27,28,59,62,68}`,
`generateSinusoid.R:16`) — real hits, not named-argument false positives (checked
each line). That contradicts CONTRIBUTING.md's "Already universal — there are
zero `=` assignments in `R/`." The claim is stale; out of scope to fix here
(zero-line diff to `R/`), but the PR description will flag it so the sentence
gets corrected separately.

The `vector_logic_linter` hits (8, e.g. `computeCumulativeCICorrelation.R:99`,
`computeInfoVal2IFC.R:{105,183}`, `generateCI.R:{180,492}`,
`generateNoiseImage.R:27`, `plotZmap.R:{132,134}`) are exactly the bug class
the issue cites — `&`/`|` where `&&`/`||` was meant. Not fixed here either
(same zero-diff constraint); called out in the PR description as the
highest-priority item on the resulting to-fix list.

## What ships

1. **`.lintr`** — the tuned `linters:` block above, plus a generated
   `exclusions:` baseline covering the current 326 lints (grouped by file,
   then by linter, listing line numbers — the format `lintr`'s own vignette
   documents for adopting it on an existing codebase).
2. **`tools/regenerate-lintr-baseline.R`** — regenerates the `exclusions:`
   block from a fresh `lint_package()` run (via
   `options(lintr.exclusions = list())` to ignore the current baseline while
   still reading `linters:` from `.lintr`), and rewrites `.lintr` in place.
   This is how a future intentionally-accepted lint gets added to the
   baseline, and how the baseline shrinks as items get fixed — rerun it and
   diff, don't hand-edit the block.
3. **`.github/workflows/lint.yaml`** — one `ubuntu-latest` job, `r-lib/actions`
   pattern matching `R-CMD-check.yaml` (checkout, `setup-r`,
   `setup-r-dependencies` with `extra-packages: any::lintr` and
   `needs: check`), then `Rscript -e 'lintr::lint_package()'` with
   `LINTR_ERROR_ON_LINT: true`. **Not** added to the required-checks ruleset —
   that needs a repo-settings write this environment cannot make (same
   constraint documented in `MAINTENANCE.md` for the stale-`man/` gate), and
   it does not need to be required for "no new lints" to hold: the baseline
   already excludes every existing lint, so any lint the job reports is by
   definition new, and a non-required red check is still visible on the PR.
4. **`.Rbuildignore`**: add `^\.lintr$` (dotfile at top level, not part of the
   built package — same reasoning as the other repo-root config files already
   listed).
5. **`MAINTENANCE.md`**: one row in the workflow table, plus enough of the
   above reasoning that the next person doesn't reconstruct it from a diff.
   Current count is 1546/1800 words — tight, so this stays short and the
   detail above lives in the PR description instead.

No changes to `R/`, `NEWS.md`, or the `.Rdata` contract — this is tooling
only, matching the issue's "zero-line diff to `R/`" framing.

## Step most likely to be wrong

The `.lintr`/`exclusions:` line-number baseline is fragile against reformatting: if
anything in `R/` or `tests/` gets reformatted before this merges, the line
numbers in the generated baseline drift and the CI job would either miss real
new lints (line shifted onto an excluded one) or false-positive (excluded
line shifted onto clean code). Mitigated by generating the baseline as the
last step, immediately before opening the PR for real, and by having
`tools/regenerate-lintr-baseline.R` as the documented fix if it drifts before
merge — rerun it, diff should be line-number-only.
