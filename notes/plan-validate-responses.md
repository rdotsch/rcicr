# Plan: reject missing and non-numeric responses (#372)

## The defect

One `NA` response makes every pixel of the CI `NA` (black PNG), with only R's own `min()`/`max()` warnings. Factor and character responses do the same. Through the batch functions the call stops in `autoscale()` with a message about masks. Logical responses are weighted 1/0, not ±1, so every `FALSE` trial contributes nothing.

## The change

1. `coerceTrialVectors()` (`R/ci-inputs.R`, used by `generateCI()` and `computeCumulativeCICorrelation()`) stops when `responses`:
   - is not numeric (factor, character, logical), saying so and, for logicals, how to recode (`ifelse(x, 1, -1)` for a 2IFC choice, `as.numeric(x)` if 1/0 weights were intended);
   - contains `NA`/`NaN` or infinite values, naming the count and up to five trial positions, as the #337 participant-ID error does.
   Dropping such trials would decide for the user what a missing response means (DECISIONS.md, "Trial alignment errors stop computation").
2. `batchGenerateCI()` / `batchGenerateCI2IFC()` check the response column of the whole table before computing any CI, naming the groups affected, like `requireParticipantIds()`.

## Behaviour change

A call that used to return an all-`NA` CI, or a CI with `FALSE` trials weighted 0, now stops. No call that returned a meaningful number before is affected, except logical responses, whose result was wrong for a 2IFC task. NEWS.md: behaviour changes section.

## Tests

One per rejection (NA, NaN, Inf, factor, character, logical) in `generateCI()`, pooled and with participants; `computeCumulativeCICorrelation()`; the batch pre-check naming the group; numeric ±1 and rating-scale responses still accepted with identical output. Each fails without the fix (`git stash push -- R/`).

## Step most likely to fail

Rejecting logicals was decided on the issue. If a caller relied on 1/0 weighting through logicals, the error message tells them the one-line conversion that reproduces it exactly.
