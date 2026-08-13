# Plan: fix progress bar under parallel execution (#178)

## The bug

`generateStimuli2IFC()` creates a `txtProgressBar` in the parent process
(`R/generateStimuli2IFC.R:182`) but ticks it inside the `foreach` body (lines 244 and 250).
Under parallel execution — the default — each worker gets its own copy of `pb`, so
`setTxtProgressBar()` updates a bar nobody can see. The user sees 0% until the job finishes.

With `ncores = 1` the loop runs in-process (via `registerDoSEQ()`) and the bar works,
so the defect is parallel-only.

## Spiked

Three approaches tested (`scratchpad/spike-progress.R`):

| approach | result |
|---|---|
| `.combine` wrapper + `.multicombine = TRUE` | Jumps 0% → 100% — `doParallel` batches all results into one combine call |
| `doSNOW` progress callback | Smooth per-task ticks, accurate count |
| `.combine` wrapper + `.multicombine = FALSE` | Smooth but off-by-one (n-1 ticks for n tasks) |

`doSNOW` is the clear winner.

## The fix

1. **Add `doSNOW` to `Imports` in `DESCRIPTION`.** It provides `registerDoSNOW()`, which
   accepts `.options.snow = list(progress = fn)` — a callback that fires in the parent
   process each time a task completes. `snow` (its only dependency) is already a transitive
   dependency of `doParallel`; both are by Revolution Analytics / Microsoft.

2. **In `startBackend()` (`R/zzz.R`)**, replace `doParallel::registerDoParallel(cl)` with
   `doSNOW::registerDoSNOW(cl)`. They do the same thing — register a `foreach` backend on a
   `parallel::makeCluster()` cluster — but `doSNOW` honours `.options.snow`.

3. **In `generateStimuli2IFC()`:**
   - Build `.options.snow = list(progress = function(n) setTxtProgressBar(pb, n))` when
     `cl` is not NULL, empty list otherwise.
   - Pass `.options.snow = opts` to the `foreach()` call.
   - Keep the two `setTxtProgressBar(pb, trial)` calls inside the loop body (lines 244
     and 250), guarded by `is.null(cl)`. `registerDoSEQ()` ignores `.options.snow`
     entirely — measured, see "Serial path" below — so removing them would trade a broken
     parallel bar for a broken serial one.

4. **In `generateCI()`**, two more loops have the same pattern:
   - The **participant-CI loop** (`R/generateCI.R:259-273`): `pb` at 259, ticked at 273.
   - The **z-map loop** (`R/generateCI.R:350-362`): `pb` at 350, ticked at 362.
   Apply the same `.options.snow` fix to both.

5. **Update `NAMESPACE`** via `@importFrom doSNOW registerDoSNOW` in `zzz.R`.

6. **Add a regression test** (`tests/testthat/test-parallel-progress.R`). The numeric tests
   all pass with the bar broken, so nothing today would catch a callback omitted from one of
   the three loops, or the serial ticks removed.

   The parent's bar is observable: `capture.output(type = "output")` sinks it, while worker
   output (`makeCluster(outfile = "")`) writes past the sink to the terminal. So captured
   percentages are exactly the ones the *user* would see. Measured on a bare
   `foreach`/`doSNOW` loop of 6 tasks: `0% 20% 40% 60% 80% 100%` at `ncores = 2`.

   Assert for both `ncores = 1` and `ncores = 2`, on `generateStimuli2IFC()` and on
   `generateCI()`'s participant-CI and z-map paths, that the captured output holds **more
   than one distinct percentage** — i.e. the bar advanced rather than jumping straight to the
   end. Match on distinct values, not an exact sequence: task completion order is not
   deterministic under a cluster, and pinning the exact ticks would make the test flaky.
   `skip_on_cran()` the `ncores = 2` cases, as `test-parallel-equivalence.R` does.

7. **Add a `NEWS.md` entry** under the development-version heading. Progress reporting during
   default (parallel) execution is user-visible behaviour. It goes in the bug-fix section,
   not "Reproducibility impact" — no numeric output changes.

## What does not change

- Numeric output: `doSNOW` dispatches to the same `parallel::makeCluster()` workers. Neither
  parallel loop in this package draws random numbers inside the loop body, so there is no
  RNG stream to diverge.
- Serial execution: `registerDoSEQ()` path is unchanged.
- The `.packages = 'rcicr'` requirement on the `foreach` call is unchanged.

## Serial path

`registerDoSEQ()` does **not** fire `.options.snow`: a 6-iteration `%dopar%` under `doSEQ`
with a counting `progress` callback fired it 0 times. The in-body ticks therefore stay,
guarded by `is.null(cl)`, and the two mechanisms never both run on the same loop.

## Step most likely to fail

The `generateStimuli2IFC()` loop ticks the bar from **two** places (lines 244 and 250) —
line 244 is on the `return_as_dataframe` early-return path, added for issue #82 because that
`return` exits the whole `foreach` body. Both need the same `is.null(cl)` guard, and the
comment at line 244 explaining why it is duplicated has to survive the edit.

## Dependency cost

`doSNOW` is 36 KB installed, has one dependency (`snow`, already transitively present), is
actively maintained, and is by the same team as `doParallel` and `foreach`. CRAN has no
issue with it. Adding it to `Imports` means it is installed automatically; it is not
optional.
