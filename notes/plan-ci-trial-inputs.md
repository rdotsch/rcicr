# Validate CI trial alignment (#294, #300)

## Verified evidence

At main 2b1313a, R/ci-inputs.R coerces trial vectors and checks only stimulus/response lengths. selectStimulusParams() indexes the saved matrix without validating IDs. R/generateCINoise.R multiplies stimuli by responses without checking row alignment. R/computeCumulativeCICorrelation.R separately indexes stimuli and has no trial-length check. These statements were checked against the current files through GitHub.

The runtime evidence in #294 and #300 records reproductions on Linux release/devel, macOS and Windows: two participant IDs recycle across twenty trials, fractional IDs truncate, and zero IDs drop rows while responses recycle. This plan adds no new runtime measurements.

## Proposed contract

- After the existing data-frame/tibble coercion, require one participant ID per stimulus/response whenever grouping is enabled. Preserve the existing all-NA no-grouping convention, including the default scalar NA; changing missing-ID semantics is outside this fix.
- Accept numeric integer-valued stimulus IDs from 1 through the selected base's saved trial count, including repeats, nonconsecutive IDs, and single trials. Reject empty, zero, negative, fractional, missing, nonfinite, and out-of-range IDs before aggregation or indexing can discard or reinterpret them.
- Reject factor, character and logical stimulus IDs explicitly rather than interpreting factor codes, row names, or logical masks as trial numbers. Document converting imported labels to verified numeric stimulus numbers, never as.numeric(factor).
- Apply the same ID and trial-length contract to cumulative correlation, while preserving its trial ordering and deliberate nonaggregation.
- At direct generateCINoise() entry, require a matrix row per response, or exactly one response for a single-trial parameter vector. Preserve the arithmetic and legacy parameter truncation.

Errors are preferable to warnings because an ignored warning would still produce a plausible CI with the wrong trial assignments. This intentionally tightens previously accepted malformed calls and must be documented as a compatibility change.

## Implementation and compatibility boundary

Keep these two issues together as one trial-alignment contract. Introduce a shared ID validator, validate before generateCI() aggregates responses, and retain validation at parameter selection so direct helper use cannot bypass it. Check participant lengths before any participant splitting or output writes. Reuse validation in cumulative correlation without refactoring unrelated computation.

Inspect every generateCINoise() caller before adding its guard: the highest-risk step is confusing a dropped single-row vector with multiple trials, or rejecting an existing valid wrapper call. Inspect factor handling before unlist() so coercion cannot erase the information required to reject ambiguous IDs.

Do not change valid numeric results, response weighting, participant ordering, image scaling, RNG behavior, or saved fields. Update roxygen/man documentation, NEWS.md Reproducibility impact, and the relevant DECISIONS.md entry within its word budget. Remove this plan on the same branch during implementation.

## Acceptance evidence required

1. Regression tests reject short (dividing and nondividing) and long participant vectors, including one-column input containers. Keep all-NA no-grouping controls and valid unequal participant groups.
2. Invalid-ID cases cover both CI and cumulative APIs, including the case where aggregate() would otherwise discard an NA ID. Factors, characters, logicals, and empty vectors have explicit errors.
3. Valid repeated/nonconsecutive IDs, one trial, imported one-column numeric data, and legacy saved parameter layouts produce unchanged numerical output.
4. Direct generateCINoise() rejects mismatched dimensions but retains matrix and single-trial-vector results.
5. End-to-end grouping tests compare individual CIs to independently selected participant trials and show the incorrect grouping differs. Exercise serial and parallel execution.
6. Prove the new regression tests fail against unchanged production code. Install the package before worker tests; regenerate docs with the pinned roxygen2.
7. Run the golden master, full suite, R CMD check and the quick release comparisons against v1.0.1 and v1.3.0. Valid output must not need a new EXPECTED entry. Document invalid-call errors through targeted compatibility tests rather than changing the numerical gate to expect crashes.
