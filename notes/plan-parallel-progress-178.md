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
   - Remove the two `setTxtProgressBar(pb, trial)` calls inside the loop body (lines 244
     and 250) — they are redundant under `doSEQ` too, since the callback fires in-process
     as well. **Verify this**: if `doSEQ` does not fire `.options.snow`, keep the in-body
     ticks guarded by `is.null(cl)`.

4. **In `generateCI()`**, the z-map parallel loop (`R/generateCI.R:449-475`) has the same
   pattern — progress bar created outside, ticked inside. Apply the same fix there.

5. **Update `NAMESPACE`** via `@importFrom doSNOW registerDoSNOW` in `zzz.R`.

## What does not change

- Numeric output: `doSNOW` dispatches to the same `parallel::makeCluster()` workers. Neither
  parallel loop in this package draws random numbers inside the loop body, so there is no
  RNG stream to diverge.
- Serial execution: `registerDoSEQ()` path is unchanged.
- The `.packages = 'rcicr'` requirement on the `foreach` call is unchanged.

## Step most likely to fail

Whether `registerDoSEQ()` fires the `.options.snow` callback. If it does not, the serial
path loses its progress bar unless we keep the in-body ticks for `is.null(cl)`. The spike
did not test this — verify before removing the in-body calls.

## Dependency cost

`doSNOW` is 36 KB installed, has one dependency (`snow`, already transitively present), is
actively maintained, and is by the same team as `doParallel` and `foreach`. CRAN has no
issue with it. Adding it to `Imports` means it is installed automatically; it is not
optional.
