# Plan: return the scaling each classification image was rendered with (#9)

## What is missing, verified

The reopened issue's claim holds on `main` at 7e65245:

- `generateCI()` returns `ci`, `scaled`, `base`, `combined` and optionally `zmap`, plus the `trial_design` attribute (`R/generateCI.R:247-255`). The scaling method is not among them. Nor is the constant `applyScaling()` uses: for `'independent'` it is computed from the CI's range and discarded (`R/generateCI.R`, `applyScaling()`).
- `autoscale()` prints its constant (`write(paste0("Using scaling factor constant:", constant), stdout())`) and does not return it; its own documentation says "printed to the console, not returned".
- With `participants` and `save_individual_cis = TRUE`, each participant's PNG is scaled by `individual_scaling`. For `'independent'` that is a different constant per participant, and none of them is kept anywhere.

So a figure made from `$scaled` or a saved PNG cannot be traced back to the constant that produced it, and two `'independent'` images cannot be told apart from two `'constant'` ones.

## The change

A `scaling` **attribute** on every classification image, like `trial_design`, never a new list field: scripts iterate the fields and expect pixel matrices, which is why `trial_design` is an attribute too.

`attr(ci, "scaling")` is `list(method, constant)`:

| method applied | `constant` |
|---|---|
| `'none'` | `NA` |
| `'constant'` | the `scaling_constant` given |
| `'independent'` | the constant computed from this CI (`0` when the CI has no range and was rendered neutral) |
| `'matched'` | `NA`; the mapping follows from `range($ci)` and `range($base)`, both returned |
| `'autoscale'` (set by `autoscale()`) | the constant shared by the whole list |

- **`method` is the method actually applied**, so an unrecognised name, which `applyScaling()` replaces by `'none'` with a warning, records `'none'`.
- **With `participants`**, the attribute also carries `individual = list(method, constant)`. Under `'independent'`, `constant` is a vector named by participant ID, in the same sorted order the individual PNGs are named by. It is recorded whether or not `save_individual_cis` is set, so the individual images can be reproduced later from the same call. It is computed in the parent after the parallel loop, from the per-participant CIs that loop already returns, so the `foreach` body does not change. That stack is unmasked: the loop masks only the copy it renders to PNG (`R/ci-compute.R:53-59`). So the mask is applied to each participant's CI with the same `applyMask()` call before its constant is derived, or a masked-out extreme would set the constant.
- **`autoscale()`** rewrites `$scaled` but deliberately leaves `$combined` as it was passed (`DECISIONS.md`, "`autoscale()` leaves `$combined` untouched"). So the record keeps both apart. The top-level `method` and `constant` always describe `$scaled`, and the PNGs `autoscale()` writes; after `autoscale()` they are `'autoscale'` and the shared constant. A `combined` element holds the record that was there before, which still describes `$combined`. It is present only when the two differ, so a CI straight from `generateCI()` has none. `individual` is carried over unchanged. A CI without a `scaling` attribute (an older version's, or one built by hand) gets `combined = NULL`, meaning unknown. Other attributes, including `trial_design`, are kept. It still prints the constant, and its documentation stops saying the constant is not returned. `batchGenerateCI()` and `batchGenerateCI2IFC()` with their default `'autoscale'` therefore return it too.
- **One source for the constant.** The `'independent'` constant is computed by one helper used by both `applyScaling()` and the record, so the recorded value cannot differ from the one that rendered the image. A test asserts `scaled == (ci + constant) / (2 * constant)` against the recorded constant.

**Rejected: a log file.** Writing one needs a required destination argument (CRAN forbids a default write path; `DECISIONS.md` records why there is none), and it would be a second copy of state that already travels with the result. An attribute survives `saveRDS()` and `save()`. Users who want a file can write the attribute themselves; the documentation shows how.

## Tests

- Each method: the recorded method and constant, and the scaled image reproduced from them. Also assert a wrong constant fails to reproduce it, so the test cannot pass vacuously.
- An unknown method records `'none'`.
- A degenerate all-zero CI under `'independent'` records `0`.
- Masked CI: the constant is computed over the unmasked pixels, as the image was.
- Participants: `individual$constant` named by participant ID in PNG order, with and without a `mask` whose masked-out pixels hold the participant's largest absolute value; with `save_individual_cis = TRUE`, each written PNG equals `(ci_p + k_p) / (2 k_p)` combined with the base, where `ci_p` is participant `p`'s CI and `k_p` their recorded constant. Unsorted participant IDs, the order #267 got wrong.
- `autoscale()` and both batch functions with `'autoscale'`: every element carries the shared constant at the top level, `scaling$combined` reproduces `$combined` from `$ci` and `$base` (with `'none'` for the batch functions, which scale with `'none'` first), `individual` survives on a `participants` CI, and `trial_design` survives. Autoscaling twice keeps the original `combined` record rather than nesting.
- `identical()` of every pixel field before and after the change, on the golden-master inputs: the attribute changes no number.

## Gate and NEWS

No numeric output changes; the gate compares fields, not attributes. Planned check after implementation, with its output quoted in the PR: `Rscript tools/compare-release-output.R --quick --ref="$(git rev-parse origin/main)"`, expecting `0 expected deviations`. NEWS: a New features entry. No Reproducibility impact entry.

## The step most likely to fail

The per-participant constants. They have to be named in the order the individual PNGs are named (`sort(unique(participants))`, the order fixed in #267), not in factor-code order or first-appearance order. A test with participant IDs given out of order, and numeric IDs that sort differently as text, holds them to it.
