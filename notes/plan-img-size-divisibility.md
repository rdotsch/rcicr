# Plan: explain an img_size that nscales cannot tile (#339)

## Problem

`generateNoisePattern()` tiles scale `s` with `2^(s - 1)` patches per side, so each patch is `img_size / 2^(s - 1)` pixels. When `img_size` is not divisible by `2^(nscales - 1)`, the patch size is fractional and the call fails deep inside the tiling with "number of items to replace is not a multiple of replacement length", which does not say what is wrong. `generateStimuli2IFC()` reaches it after reading the base images.

## Change

1. **`generateNoisePattern()` checks first** that `img_size %% 2^(nscales - 1) == 0`, and otherwise stops with a message naming the constraint and two ways out: the nearest valid sizes (the multiples of `2^(nscales - 1)` below and above, the lower one only if it is positive), and the largest `nscales` that tiles this `img_size`. For example, `img_size = 100, nscales = 5`: "use img_size 96 or 112, or nscales = 3 or fewer".
2. **Every call that works today still works:** the check rejects exactly the calls that already fail (below). `generateStimuli2IFC()` needs no change; it calls `generateNoisePattern()` before creating the folder or taking the `.Rdata` lock, so nothing is written.
3. **`NEWS.md`, under "Bug fixes"**, last (a message-only fix).
4. **`?generateStimuli2IFC`** (`img_size`, `nscales`) states the constraint.

## Verified before planning

- On `main`, 110 combinations of `img_size` (16, 24, 30, 32, 48, 50, 64, 96, 100, 120, 128), `nscales` 1 to 5 and both noise types: all 20 non-divisible ones fail with that error, and all 90 divisible ones return a pattern of the right size with finite values. So divisibility is exactly the failure condition, and the check moves no working call.
- Every release-gate configuration in `tools/compare-harness.R` uses a divisible size (15 configurations, including `nscales = 6` at 128), so the gate should report no deviation.
- `generateNoisePattern()` has two callers: `generateStimuli2IFC()`, and `simulateNoiseIntensities()`, which uses the default `nscales = 5` with its `img_size` argument.

## Tests

- `img_size = 100, nscales = 5` stops with the constraint, `96 or 112` and `nscales = 3`;
- a size below `2^(nscales - 1)` (for example 10 with `nscales = 5`) offers only the upper size;
- the grid above: every divisible combination still returns a pattern, every other one stops with the new message;
- `generateStimuli2IFC()` with such a size stops with the message and leaves `stimulus_path` absent.

The failure tests must fail on the current code; checked with `git stash push -- R/`.

## The step most likely to fail

**The grid being wide enough.** The claim that no working call breaks rests on 110 combinations. The reasoning backs it: `matlab::repmat(p, scale)` produces `scale * nrow(p)` rows, which equals `img_size` only when the patch size is a whole number.

## Out of scope

- Other invalid `img_size` or `nscales` values (non-integer, zero, negative): not measured here.
