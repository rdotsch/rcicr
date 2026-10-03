# Plan: one set of checks for every reference path (#377)

## The defect

The per-base and subset paths validate `iter` and the reference; the default (shared-parameter) path does not:

- `generateReferenceDistribution2IFC(iter = 2.5)` runs 2 iterations; `iter = 0` and `NA` fail inside `txtProgressBar()` or an `if`.
- `subsetReference()` stops on a MAD of 0; `sharedReference()` and `baseReference()` divide by it (infinite or `NaN` InfoVal). `sharedReference()` also does not reject a stored reference with empty or non-finite norms.

## The change

1. `validateIter(iter)`, called at the top of `generateReferenceDistribution2IFC()` before any dispatch, replacing the two copies in `generateBaseReference()` and `generateSubsetReference()`. Not in `computeInfoVal2IFC()`: there `iter` is unused when a stored reference is found, and a call that works today with a cache must keep working.
2. `referenceSummary(norms, note, what)` returning `list(median, mad, iter, note)`, used by all three resolvers: it rejects empty or non-finite norms (with the existing `force_gen_ref_dist = TRUE` advice) and a MAD of 0.

## Behaviour change

Only calls that errored opaquely, ran a truncated `iter`, or returned `Inf`/`NaN` change; every finite InfoVal is identical.

## Tests

`iter` of 2.5, 0, -1, `NA`, `c(1, 2)` on the shared path stop with "iter must be a positive integer" before simulating. A stored shared reference of constant norms, or holding `NA`, stops in `computeInfoVal2IFC()`. Existing per-base and subset messages unchanged.
