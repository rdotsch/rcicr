# Plan: InfoVal for many classification images in one call (#85)

## Why a function, when the reference is already stored for reuse

Measured here (R 4.3.3, one core) by `analyses/infoval-batch-cost.R`, committed with this plan, on a 300-trial 512px stimulus file of 39 MB. Its output:

```
file MB: 39.14279
reference: 6.0s; cached-hit call: 0.84s; bare load(): 0.42s
forced/unstorable reference per call: 4.8s
```

| | per CI today |
|---|---|
| `computeInfoVal2IFC()` with the reference already stored | 0.84 s. A bare `load()` of the file is 0.42 s, and the call loads it twice (`selectReferenceBase()`, then `loadRdata()`) |
| the same when the reference is simulated on every call (a read-only file, `response_seed`, or `force_gen_ref_dist = TRUE`) | 4.8 s |

So a loop over 100 participants spends about 80 s loading one file 200 times, or about 8 minutes simulating one reference 100 times. The batch loads the file once and resolves each distinct reference once. It is also where per-participant `reference_stimuli` (each CI scored over the stimuli it was built from) becomes one argument instead of a hand-written loop.

## The function

`batchComputeInfoVal2IFC(target_cis, rdata, iter = 10000, force_gen_ref_dist = FALSE, response_seed = NULL, baseimage = NULL, reference_stimuli = NULL, reference_method = c("gram", "images"))`

- **Name** (the maintainer's choice). It pairs with `batchGenerateCI2IFC()`, whose result is its natural input, and sits beside it in the reference index.
- `target_cis`: a list of classification images as returned by `generateCI()`. Every element is validated as a CI: a list whose `ci` element is a numeric matrix. A single CI passed by mistake is recognised by its own `ci` element being a numeric matrix, not by the name, since `list(ci = first, control = second)` is a valid batch. It is refused with a pointer to `computeInfoVal2IFC()`.
- Every other argument means what it means in `computeInfoVal2IFC()` and applies to every CI, except `reference_stimuli`, which is one of:
  - `NULL`: every CI against the full set, as `computeInfoVal2IFC()` does;
  - a vector of stimulus numbers: the same subset for every CI;
  - a list with one element per CI (a `NULL` element meaning the full set), for example `lapply(cis, function(ci) attr(ci, "trial_design")$stimuli)`.
- `baseimage` is one label for the whole batch. CIs for different base images go in separate calls, which is how their references are stored anyway.
- **Returns** a numeric vector, one InfoVal per CI, named by `names(target_cis)`.

**The contract:** the result is `identical()` to `stats::setNames(vapply(seq_along(cis), function(i) computeInfoVal2IFC(cis[[i]], rdata, ..., reference_stimuli = <that CI's>), numeric(1)), names(cis))` with the same arguments, names included (and no names when `cis` has none), including under `response_seed` and `force_gen_ref_dist = TRUE`. Both of those give the same reference on every call (one seeded stream, one stimulus stream), so resolving it once and reusing it changes no number. That is what the tests pin.

## Output

- One reference-resolution line per distinct reference ("Using reference distribution found…", or the simulated/saved notes), exactly as `computeInfoVal2IFC()` prints them.
- One line per distinct reference with its median, MAD and iterations, in place of `computeInfoVal2IFC()`'s per-CI line, since the per-CI numbers are the return value.
- The trial-design check (`reportTrialDesign()`) runs per CI but reports **once**: a single message giving how many CIs are scored against a reference over different stimuli, naming up to five, and the same `reference_stimuli` fix. The message about CIs that average participants or repeats is aggregated the same way. With 100 participants, 100 identical messages would bury the one that differs.

## Implementation

`computeInfoVal2IFC()` has three paths (shared, independent base, subset), each ending in "resolve the reference, then score". Each is split at that point into an internal resolver and a scoring step. The resolver returns a summary, `list(median, mad, iter)`, not the norms: a repopulated `ref_lookup` row supplies only those three numbers, and scoring needs nothing else. `computeInfoVal2IFC()` becomes resolver plus score, with its printed lines and messages unchanged. The batch groups CIs by their canonical `reference_stimuli` (`canonicalReferenceStimuli()`, so `c(3,1,2)` and `1:3` share a reference) and calls the resolver once per group.

The shared path's `ref_lookup` block (empty since 2018) stays inside its resolver, unchanged, and a hit on it would serve the whole batch like any other resolved reference.

## Tests

- Parity: `identical()` to the named `computeInfoVal2IFC()` loop above, on a named list straight from `batchGenerateCI2IFC()` and on an unnamed list, for the full set, a shared subset, per-CI subsets (including a `NULL` element and two CIs sharing a subset in a different order), an independent-base file with `baseimage`, `response_seed`, `force_gen_ref_dist = TRUE`, and a read-only file.
- Resolution count: mocking `generateReferenceDistribution2IFC()` to count calls, a read-only file with 5 CIs over 2 distinct subsets simulates exactly 2 references (the loop simulates 5). Stored-reference case: the file is loaded a fixed number of times, independent of the number of CIs.
- Messages: 3 mismatched CIs out of 5 give one message naming those 3; matched CIs give none.
- Refusals: a single CI; `reference_stimuli` a list of the wrong length; an element that is not a CI. Accepted: a batch with an element named `ci`.
- The existing `computeInfoVal2IFC()` tests pass unchanged; the refactor adds no expectation to them.

## Gate

Planned check, run after the implementation, with its output quoted in the PR: `Rscript tools/compare-release-output.R --quick --ref="$(git rev-parse origin/main)"`. It must report `0 expected deviations` and no unexpected one, since the refactor has to leave `computeInfoVal2IFC()` numerically inert.

## NEWS and docs

New features entry; an example in the new function's help, and a pointer from `computeInfoVal2IFC()`'s "Matching the reference" section. No Reproducibility impact entry: no existing call returns a different number.

## The step most likely to fail

The refactor of `computeInfoVal2IFC()`. Its shared path loads the stimulus file into its own frame and depends on `captureArgs()` to protect its arguments from saved names. Moving that code into a resolver has to keep the protection, or a saved field such as `iter` in an old file could silently replace the caller's value. The existing collision tests in `test-computeInfoVal2IFC.R` and `test-reference-from-saved-noise.R` cover this and must stay green without edits.
