# Plan: group tibbles correctly in the batch functions (#336)

## Problem

`batchGenerateCI()` and `batchGenerateCI2IFC()` select groups with `data[, by]`. On a tibble that stays a one-column tibble, so the loop runs once with the whole column as `unit`, and `data[, by] == unit` recycles it element by element. With more than one group the call returns one CI, built from alternating trials of different groups and named after the first. Affected: 1.1.0 through 1.4.1, measured by installing v1.0.1 (error) and v1.1.0 (one wrong CI). #336 has the reproduction.

## Change

1. **Both functions read columns as vectors**, with `data[[by]]`, `unitdata[[stimuli]]` and `unitdata[[responses]]`, in all six places that now use `[, col]`: the `NA` filter, the progress bar length, the loop, the row selection, the file name and the two arguments passed on. For a `data.frame` these return exactly what `[, col]` returns today, so its results do not change.
2. **Rows are selected by comparing a vector with a single value**: `data[data[[by]] == unit, , drop = FALSE]`, after the existing `NA` filter.
3. **`NEWS.md`, under "Reproducibility impact"**, first in the section:
   - who is affected: `batchGenerateCI()` or `batchGenerateCI2IFC()` called on a tibble with more than one group, on 1.1.0 to 1.4.1;
   - how to tell: the returned list had one element although the data held several groups;
   - what to do: recompute with this version, or on an older one convert with `as.data.frame()` first. Nothing needs recollecting.
   - It also says that 1.1.0's statement that these functions accept tibble columns was true only for a single group. The released 1.1.0 section itself stays as it is.
4. **README**: a short callout beside the existing 1.3.0 one, since this changes published numbers for the people it affects.

## Verified before planning

- The loop is unchanged since v1.0.1; v1.0.1 errors on a tibble and v1.1.0 returns one wrong CI (both installed from their tags and run).
- For 1, 2, 3 and 5 groups a tibble returns exactly one CI; with one group it equals the `data.frame` result. So "one CI although several groups" identifies every affected call, and nothing else.
- #110 is a real user calling `batchGenerateCI2IFC()` on `readr::read_csv()` output, so the input occurs in practice.
- The release gate runs `batchGenerateCI()` on a `data.frame` only (`tools/compare-harness.R:374`), so it should report no deviation. A tibble configuration cannot be added to the gate: v1.0.1 errors on it, and a reference-side error aborts the gate.

## Tests

For both functions, on the same data as a `data.frame` and as a tibble:

- 2 and 3 groups give the same names and `identical()` CIs;
- a `by` column with `NA` values drops those rows in both;
- a factor `by` column gives the same result as a character one;
- file names, and so the list names, are unchanged for a `data.frame`.

Each new test must fail on the current code; checked with `git stash push -- R/`.

## The step most likely to fail

**That `data.frame` results really do not move.** The reasoning is that `data[[by]]` and `data[, by]` return the same vector for a `data.frame`; the release gate is what measures it, on the first implementation push. A deviation there would mean the change reaches beyond tibbles and must stop the PR.

The `NA` filter is not a risk: `data[!is.na(data[, by]), ]` already keeps the same rows for a tibble as for a `data.frame` (measured: rows 1, 3 and 4 of `c("p1", NA, "p2", "p2", NA)` in both).

## Open question for review

Is a `NEWS.md` entry plus a README callout enough, or should this get an advisory article like the individual-CI one (`vignettes/articles/rcicr-individual-ci-advisory.Rmd`)? The individual-CI case needed one because telling whether you were affected took work. Here the tell is a single glance at the length of the returned list.

## Out of scope

- data.table input. It probably errors today rather than returning a wrong result; not measured, and not claimed as fixed.
- #337 (partly missing participant IDs) and #333: separate changes.
