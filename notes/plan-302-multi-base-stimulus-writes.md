# Plan — #302: `return_as_dataframe` skips stimulus writes for all but the first base

## The defect

`generateStimuli2IFC.R:253-262` returns `as.vector(trial_noise)` from inside the loop over base
faces. `return()` exits the whole `%dopar%` body, not the inner loop, so every base after the
first never reaches its `png::writePNG()` calls. The `.Rdata` still records parameters for every
base, and the call reports success, so the stimulus set is silently incomplete.

The roxygen contract at `generateStimuli2IFC.R:27` states the opposite: "Stimuli are still
written to disk for every base image either way."

## What was measured, and how

Locally, R 4.3.3, `rcicr` installed from this working tree at `f110a85`
(`tools/audit-runtime.R`'s `#302` case measured the same on four platforms in #304).

| claim | result |
|---|---|
| Two bases, two trials, `return_as_dataframe = TRUE`, `save_as_png = TRUE` | **4 PNGs, not 8**; all four missing are `base2`'s |
| Same, across `use_same_parameters` TRUE/FALSE × `ncores` 1/2 | 4 PNGs in all four combinations |
| Returned frame's contents | `max|frame - base1 noise| = 0`, `max|frame - base2 noise| = 0.447` |
| `generateNoiseImage()` calls, `save_as_png = FALSE`, 3 bases × 4 trials | **4**, not 12 |

The last row is the one that constrains the fix. `generateReferenceDistribution.R:165` calls this
function with exactly `return_as_dataframe = TRUE, save_as_png = FALSE`, and the early `return`
is what currently stops it computing noise for bases it will never use. Simply continuing the
loop would multiply that work by the number of base images — on the path #309 has just been
optimising.

That path is now reached only for shared parameters: #305 routes the independent-base case
through `generateBaseReference()`, which does not call this function. With
`use_same_parameters = TRUE` every base has the same noise, so returning the first base's is
correct there.

## The change

In the `%dopar%` body:

1. Capture the first base's noise into a local instead of returning from the loop, keeping
   `as.vector()` and its position (before any per-base reassignment of `trial_noise`), so the
   returned frame is unchanged value-for-value.
2. Iterate `names(base_faces)` when `save_as_png` is TRUE and only the first name when it is
   FALSE. Nothing after the first base contributes anything on the no-PNG path, so this holds
   the measured cost at one `generateNoiseImage()` call per trial.
3. End the body with the captured noise when `return_as_dataframe` is TRUE, and `trial`
   otherwise — the value `.combine`/`.final` already expects.
4. Drop the duplicated `setTxtProgressBar()` at `:260` and its comment. With no early return the
   tick at `:267` runs on every path, so the bar advances once per trial exactly as now.

The comment at `:254-259` documents the control flow being removed and goes with it.

Rejected: computing the returned noise unconditionally and picking it up after the loop. With
`use_same_parameters = FALSE`, `trial_noise` holds the *last* base's noise once the loop ends,
which would silently change what the frame returns.

## Expected impact

- **Numeric output: none.** The returned frame keeps the first base's noise; `base1`'s PNGs are
  written from the same values in the same order. No RNG is consumed inside the loop —
  every parameter is drawn before it — so the stream is untouched.
- **New files.** `base2..n` PNGs now appear where a two-base call previously wrote four. That is
  the fix.
- **Release gate: expected clean, no new `EXPECTED` entry.** `tools/compare-harness.R:221` calls
  this function with `return_as_dataframe` at its default FALSE, so the affected combination is
  not in the battery. If the gate does flag something, the change is wrong — it is not a
  deviation to write up.

## Tests

New block in `tests/testthat/test-fixed-bugs.R`, asserting intended behaviour:

- Two named bases × two trials with both flags on writes **eight** PNGs, named for both bases,
  under `use_same_parameters` TRUE and FALSE.
- `base2`'s decoded pixels differ from `base1`'s, so the test cannot pass by writing the same
  image twice.
- `base1`'s decoded pixels equal `(noise + 0.3) / 0.6` combined with its base face, computed
  independently from the saved parameters — pinning that the first base's output did not move.
- The returned frame is 1024 x 2 and equals the first base's noise.
- With `save_as_png = FALSE`, `generateNoiseImage()` is called once per trial regardless of base
  count (via `trace()`, serial only), pinning the property in the fourth row of the table above.

Each will be checked to fail without the change with `git stash push -- R/`, per
`CONTRIBUTING.md` → "Pull requests".

## Most likely to fail

The call-count test. `trace()` on a namespace function only sees calls in this process, so it is
meaningful under `ncores = 1` and silently vacuous in a worker. It is written serial-only for
that reason; if it proves brittle across platforms the property is still covered by the
reproduction in #304, and the test comes out rather than being weakened into a shape assertion.

## Out of scope

`stimulus` and `combined` are computed even when `save_as_png` is FALSE and discarded. With the
first-base shortcut that is one wasted pair of matrix operations per trial — unchanged from
today, and not this fix's to widen into.

## NEWS.md

One entry under a new **Behaviour changes** heading in the development section, above
Performance and below Reproducibility impact. It names who is affected (multi-base callers using
`return_as_dataframe = TRUE` with PNG saving) and what to do: the `.Rdata` already holds every
base's parameters, so re-running the same call with the same `seed` writes the missing stimuli
without invalidating anything already collected.
