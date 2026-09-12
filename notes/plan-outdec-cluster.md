# Plan: start parallel workers under any `OutDec` / `scipen` (#316)

## The failure, reproduced

Every entry point that spawns workers stalls and then fails for a user who has set **both**
`options(OutDec = ",")` and a `scipen` negative enough to force scientific notation. Measured
here, two workers, R 4.3.3, default `options(timeout)` of 60:

```
baseline                           0.2s  ok
scipen = -10 only                  0.2s  ok
OutDec = ',' only                  0.2s  ok
OutDec = ',' + scipen = -10      240.1s  ERROR: Cluster setup failed. 2 of 2 workers failed to connect.
```

Neither option alone does it. It is a slow failure rather than a hang, and at the default
`ncores = detectCores() - 1` the wait is longer still while the error names the workers rather
than the cause.

## Mechanism

`parallel` renders the worker's connection arguments into a command line as text, and the
parent's formatting options apply:

```
port rendered as: 1,1e+04     # OutDec = ",", scipen = -10
port rendered as: 11000       # restored
```

The child parses that to `NA` and dies in `makeSOCKmaster()` on
`setup_timeout >= 0 is not TRUE`, with `NAs introduced by coercion` ahead of it; the parent
then waits out its full timeout. This is R's `parallel`, not anything this package does — but
`rcicr` is where the user is standing when it happens.

## The change

`startBackend()` in `R/parallel.R` is the only place this package builds a cluster —
`reference-base.R:118`, `zmap-compute.R:36`, `ci-compute.R:28` and `generateStimuli2IFC.R:195`
all route through it — so the fix has one site. Neutralise the two formatting options across
`parallel::makeCluster()` and restore the caller's on exit:

```r
old <- options(scipen = 0, OutDec = ".")
on.exit(options(old), add = TRUE)
cl <- parallel::makeCluster(ncores, outfile = "")
```

**Safe for numeric output.** `OutDec` and `scipen` govern rendering, never arithmetic, and the
window covers cluster construction only — no computation and no random draw happens inside it.
The reproducibility gate and `test-regression-baseline.R` should both be untouched by this;
that is a claim the gate settles, not one to assert here.

## Test

Assert the fix, not the failure: reproducing the stall costs four minutes of CI for a result
already recorded above. Under both options set, a cluster starts promptly and the caller's
`OutDec` and `scipen` survive the call unchanged. Registering the backend and stopping it
immediately keeps it cheap.

`git stash push -- R/` confirms the test fails without the fix — with the four-minute wait,
once, locally rather than in CI.

## Most likely to fail

- **`on.exit()` placement.** `startBackend()` returns the cluster to a caller that stops it
  later; the restore must fire when `startBackend()` returns, not when the cluster dies. Adding
  `add = TRUE` to an existing `on.exit()` in that function would silently change teardown
  order — there is none there today, and the fix must not introduce that coupling.
- **The serial path.** `ncores < 2` returns before `makeCluster()`. The option reset must not
  straddle the `registerDoSEQ()` branch, where it would be dead weight.
- **CI cannot see the bug.** Runners have `OutDec = "."`, so a regression here would be
  invisible to every existing test; the new test is the only thing that will hold the fix in
  place.

## Not in scope

`#314` changed two tests from `scipen = -10` to `-9` to avoid this stall. Those constants are
load-bearing for the reason recorded in #316 and stay as they are; reverting them is a separate
question from fixing the cause.
