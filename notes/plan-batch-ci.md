# Plan: one implementation for the batch CI functions, with participants nested in a group (#87)

## Decided with the maintainer

`batchGenerateCI()` and `batchGenerateCI2IFC()` **stay supported and are not deprecated**. They are in published analysis scripts, their output is the natural input to `batchComputeInfoVal2IFC()`, and nothing replaces them yet: `generateCI(participants = ...)` returns only the group CI, never the individual ones. The issue's deprecation and `individual_scaling = 'autoscale'` items are therefore not done. `DECISIONS.md` records that, so the issue's ticked boxes are not taken as a plan later. What remains is the duplication, and the issue's one missing feature: CIs per condition, with participants nested in it.

## What is verified

- The two bodies are the same loop (`diff` of `R/batchGenerateCI.R` and `R/batchGenerateCI2IFC.R` on 7e65245). They differ in argument order (`label` before `antiCI` in one, last in the other), in the `by.levels` local, and in calling `generateCI()` versus `generateCI2IFC()`.
- `generateCI2IFC()` passes `scaling` and `constant` (as `scaling_constant`) and everything else the batch loop gives it straight to `generateCI()`. It leaves `participants` at its default `NA`, as `batchGenerateCI()` passes explicitly. So both batch functions already compute through the same call. The consolidation is the duplicated loop, not a numeric difference.

## The change

1. **One internal function, `batchCIs()`**, holding the loop, the `NA`-group removal, the tibble-safe `[[` column access (#336), the file naming and the autoscale step. Both exports keep their signatures and argument order exactly, and call it. Each keeps its own roxygen page.
2. **A new last argument on both, `participants = NULL`**: the name of a column of participant IDs. It is appended, never inserted: a formal in the middle rebinds every later positional argument in existing scripts. When given, each `by` group's CI is `generateCI(..., participants = <that group's IDs>)`: one CI per participant in the group, then their average, which is what the issue asked for with participants nested in condition. When `NULL`, the call is exactly today's.
   - Individual CIs are not written: there is no `save_individual_cis` passthrough. Scripts that want them can call `generateCI()` per group.
   - `participants` naming a column that does not exist stops before any CI is computed, naming the column.
3. **No other behaviour changes.** The progress bar, names, PNG file names and the autoscale step stay as they are.

## Tests

- **Parity before and after:** for both functions, on the existing fixtures, every pixel field and name identical (`identical()`) to the output of the code on `main`, with `'autoscale'`, `'independent'`, `'constant'`, a `label`, `antiCI = TRUE`, a tibble and a group of `NA`. The `main` values are computed in the test from `generateCI()` calls per group, since that is what both bodies do; a mutation to `batchCIs()` must turn them red.
- **The two exports agree** with each other, given the same named arguments.
- **Positional calls still bind as before:** a call using the full positional order of each signature gives the same result as the named call, so the appended argument cannot have shifted anything.
- **Participants nested in a condition:** each condition's CI equals `generateCI(participants = )` on that condition's rows, and differs from the pooled CI when participants contribute unequal numbers of trials. Two participants with the same ID in different conditions stay separate, because the grouping is per condition.
- **A missing `participants` column** stops before anything is computed.

## Gate, NEWS, DECISIONS

The gate's `batch` extra runs `batchGenerateCI()` with `autoscale()` (`tools/compare-harness.R:368-379`); `batchGenerateCI2IFC()` is covered by the tests above, not by the gate. Expected `0 expected deviations` against `main`, quoted in the PR once run. NEWS: a New features entry for `participants`. `DECISIONS.md`: one short entry on keeping the batch functions instead of deprecating them. It is at its budget, so an equal number of words comes out elsewhere.

## The step most likely to fail

Argument binding. The two exports put `label` in different positions, so a shared internal function must be called by name from both, never by position. The positional-call test is there to catch it.
