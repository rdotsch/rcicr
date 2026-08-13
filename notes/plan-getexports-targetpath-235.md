# Plan: `generateCI()` fails on a dead branch's free variable (#235)

## The bug

```r
generateCI(1:6, responses, "base", rdata,
           save_as_png = FALSE,               # nothing is written
           participants = c(1,1,2,2,3,3),
           n_cores = 2)                       # the default on any multi-core box
#> Error in get(s, env, inherits = FALSE) :
#>   argument "targetpath" is missing, with no default
```

`foreach::getexports()` statically scans the loop body for free variables and `get()`s each
one to ship to the workers. The participant-CI body references `targetpath` at
`R/generateCI.R:292`, inside `if (save_individual_cis)`. `targetpath` is a *required argument
with no default*, so when the caller legitimately omits it — nothing is being written — the
`get()` fails. The branch never runs. The scan does not care.

Found while writing the #178 tests; it reproduces on `main` and predates the
`doParallel` → `doSNOW` swap, since `getexports()` is `foreach`'s and is the same under both.

## Scope, measured not assumed

- **Only `targetpath`.** It is the sole required-and-absent argument reachable from a
  `%dopar%` body (`:292`). `zmaptargetpath` appears only in validation (`:128`) and at
  `:390`, both outside every loop; the z-map loop body references neither.
- **Only the participant path.** The z-map loop body is `generateNoiseImage()` plus the
  progress tick.
- **Only `n_cores > 1`** — which `default_ncores()` makes the default. `registerDoSEQ()`
  exports nothing, so the serial path works.

So the failing combination is the *default* one for anyone passing `participants` who does
not want individual CIs written.

## The fix, and the one that was rejected

Both options named in the issue were spiked (`scratchpad/spike-235.R`):

| | dead branch (`save_individual_cis = FALSE`) | live branch (`TRUE`) |
|---|---|---|
| today | **errors** | works |
| `.noexport = "targetpath"` | fixed | **`object 'targetpath' not found`** |
| resolve to a local | fixed | works |

**`.noexport` is wrong** and must not be used: it fixes the reported call by breaking one
that works today — saving individual CIs under `n_cores > 1`. It suppresses the export
unconditionally, so the live branch then finds nothing.

The fix is to give `targetpath` a bound value before the loop:

```r
targetpath <- if (missing(targetpath)) NULL else targetpath
```

`NULL` is right rather than `""`: the branch that would use it cannot run when it is `NULL`
(`save_individual_cis` is `FALSE`, which is what made `targetpath` optional), and `NULL`
fails loudly rather than writing to a relative path if that reasoning is ever wrong.

## Where the line goes — this is the whole of the risk

It must sit **after the `missing()` validation at `:121`-`:133` and before
`captureArgs()` at `:160`.**

- **After the validation**, because assigning to an argument makes `missing()` unreliable
  for it — measured: `missing(x)` is `TRUE` before `x <- ...` and `FALSE` after. Placing the
  line above `:121` would silently disable the "no targetpath was given" error and let a
  `save_as_png = TRUE` call run on to `paste0(NULL, '/', filename)`.
- **Before `captureArgs()`**, because that is strictly better than after. `captureArgs()`
  deliberately skips required-and-absent arguments (`mget()` raises on them); once
  `targetpath` is bound, it is captured like any other and restored after `load(rdata)`.
  That closes the same `load()`-frame hole for `targetpath` that the guard exists for —
  verified that `mget(names(formals()))` sees `x` once assigned.

## What does not change

- No numeric output. This adds a binding; it changes no computation, and the release gate
  should report `max|d| = 0`.
- Serial runs, which never failed.
- `batchGenerateCI()`, which forwards its own missing `targetpath` — the reason
  `captureArgs()` skips required-and-absent arguments at all. Worth an explicit check that
  it still works, since this is precisely the argument it forwards.

## Step most likely to fail

`batchGenerateCI()`'s forwarding. It relies on `targetpath` staying missing across the call;
binding it inside `generateCI()` is downstream of that, so it should be unaffected — but it
is the one place where "argument is missing" is load-bearing rather than accidental, and the
interaction has not been run.

## Tests

`tests/testthat/test-fixed-bugs.R`, which is where regression tests for this kind of defect
live. Four cases, since only some fail today:

1. `participants` + `n_cores = 2` + no `targetpath` — the reported bug.
2. The same at `n_cores = 1` — passes today; pins that the fix does not break it.
3. `save_individual_cis = TRUE` + `n_cores = 2` + a real `targetpath` — passes today, and is
   what `.noexport` would have broken. **This is the case that must not regress.**
4. `batchGenerateCI()` forwarding a missing `targetpath`.

`skip_on_cran()` the `n_cores = 2` cases, as `test-parallel-equivalence.R` and
`test-parallel-progress.R` do.

Then **remove the workaround in `test-parallel-progress.R`**: its `participant_ci()` helper
passes `targetpath = tmp` purely to dodge this bug, with a comment saying so. Once fixed, it
should stop passing it — otherwise the test no longer covers the call shape users make.

## Also

- `NEWS.md`, bug-fix section. A documented call erroring out under default arguments ranks
  above the message-only fixes.
- `cran-comments.md` if the `skip_on_cran()` file count changes — it will not, since these
  go in an existing file, but the count is a claim made to the reviewer.
