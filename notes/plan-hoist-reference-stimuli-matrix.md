# Plan: hoist reference stimulus matrix conversion

Issue: https://github.com/rdotsch/rcicr/issues/306

## Verified state

The work starts from main at `abbb7086459cde5f8df2d3309784020a2d853c11`, the squash merge of #308:

- Repository guidance: [AGENTS.md at the base commit](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/AGENTS.md), [CONTRIBUTING.md at the base commit](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/CONTRIBUTING.md), and [DECISIONS.md at the base commit](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/DECISIONS.md) require preserving published numeric output, planning behaviour changes before implementation, and checking the release gate.
- The shared path [converts `stimuli` inside the simulation loop](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/R/generateReferenceDistribution.R#L193-L209). `stimuli` is the unchanged data frame returned by stimulus regeneration.
- The independent path [already receives a matrix from `independentReferenceNoise()`](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/R/reference-base.R#L1-L41), so its loop [does not need a second conversion](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/R/reference-base.R#L43-L66). It will remain unchanged.
- Existing independent-base tests construct an independent pixel oracle and cover seeded norms, base isolation, and non-default `nscales` response-stream offsets ([test-independent-reference.R](https://github.com/rdotsch/rcicr/blob/abbb7086459cde5f8df2d3309784020a2d853c11/tests/testthat/test-independent-reference.R#L18-L75)). Those assertions use `expect_equal(..., tolerance = 1e-12)`; this plan does not describe them as bit-identical.
- The issue reports the measured repeated conversion cost and requires bit-identical norms for a fixed `response_seed`, including a non-default `nscales` or `noise_type`, with green regression and release gates and no new `EXPECTED` entry.

## Change

In `R/generateReferenceDistribution.R`, reassign `stimuli <- as.matrix(stimuli)` once immediately after `generateStimuli2IFC(..., return_as_dataframe = TRUE, ...)` returns and before response simulation. Reassignment avoids adding a scratch object to the function frame, which is later persisted by the existing `save(list = setdiff(ls(...), internals), ...)` path. The loop retains the existing response draw, multiplication, division, and norm statements, with `as.matrix(responses)` still per-iteration.

Do not modify `R/reference-base.R`: independent reference noise is already a matrix, and the independent path's existing seeded oracle covers it.

## Validation

Add deterministic tests alongside the reference-distribution tests:

1. Run the shared path with a fixed `response_seed` and a non-default `nscales`, reconstruct the same stimulus matrix and response stream independently, and assert `expect_identical()` for the complete norm vector and final `.Random.seed`. Run the same check with the legacy loop expression (matrix coercion repeated per iteration) to make the before/after equivalence explicit.
2. Where practical, use a scoped custom S3 class on a mocked stimulus data frame and count `as.matrix` dispatches. Assert the conversion count is one for the hoisted path, with no wall-clock threshold and no global namespace mutation. If the package's test harness cannot scope this method safely, omit the counter rather than add a brittle structural or AST test; the direct diff and exact numerical/RNG assertions remain the checks.
3. Keep the existing independent-base oracle tests, including `nscales = 3`, unchanged to verify that the untouched path and saved-field compatibility remain intact. Explicitly load the `.Rdata` after a saved reference run and compare the established saved fields, including the existing cache fields.

Run the targeted reference tests, `tests/testthat/test-regression-baseline.R`, and the repository's quick reproducibility gate. Confirm no new `EXPECTED` entry is needed. The full release battery is reserved by the workflow for release versions/tags; report it as not run if no release trigger is available. Do not change `NEWS.md` or `DECISIONS.md` because the operation changes neither arithmetic nor the RNG stream.

## Most likely failure

The main risk is accidentally changing the object saved back into the stimulus file or advancing the random stream in a different order. The reassignment preserves the existing `stimuli` name, and the tests compare both every norm and the final `.Random.seed`; saved-field compatibility is checked by loading a file after reference generation and comparing the established fields.
