# Plan: split `generateCI()` (#184)

`R/generateCI.R` is 589 lines. The exported function's body is 305 of them (lines 114-418)
and mixes argument validation, `.Rdata` loading, parameter selection, CI computation over two
different designs, masking, scaling, PNG writing, two z-map methods and parallelism.

**Nothing about the public function changes.** Not the signature, not the argument meanings,
not a single number it returns. This is an internal reorganisation, and the golden master
(`test-regression-baseline.R`) plus the release gate are what hold that claim to account.

## What is actually in the 305 lines

Measured on `3ed6349`, the base of this branch:

| lines | what | already extractable? |
|---|---|---|
| 114-163 | write-path checks, tibble coercion, length check | partly — `missing()` must stay in the frame |
| 165-203 | `captureArgs`/`load`/`list2env`, four `exists()` checks, `s` -> `p` | yes |
| 205-214 | base image lookup + error | yes, pure |
| 216-222 | response aggregation (collapse repeats) | yes, pure |
| 224-252 | parameter row selection + 4096 -> 4092 truncation | yes, pure |
| 254-259 | `antiCI` sign flip | trivial, leave inline |
| 261-324 | CI computation: single-shot vs per-participant `foreach` | yes, riskiest |
| 326-340 | mask / scale / combine / save | already helpers, leave |
| 342-410 | two z-map methods + `plotZmap()` call | yes |
| 412-418 | return assembly | trivial, leave inline |

The four helpers the issue names (`applyMask`, `applyScaling`, `combine`, `saveToImage`)
already exist at lines 420-589, along with `hasMask`. **They are not moving and not being
renamed.** `R/plotZmap.R:136-137` calls `hasMask()` and `applyMask()`, and
`test-plotZmap.R`, `test-error-paths.R` reach both through `rcicr:::`. Renaming any of them
breaks those call sites for no gain.

## Target shape

Follows the per-concern file convention `R/rdata.R` and `R/parallel.R` already set, rather
than growing one file:

**`R/ci-inputs.R`** (new) — pure, no I/O, no globals:

- `coerceTrialVectors(stimuli, responses, participants)` -> list of three vectors; does the
  `unlist()` coercion and the equal-length check.
- `selectBaseImage(base_faces, baseimage)` -> matrix; raises the "did not contain any
  reference to base image label" error naming the available keys.
- `aggregateResponses(stimuli, responses)` -> list(stimuli, responses); the `aggregate()`
  collapse of repeated presentations.
- `selectStimulusParams(stimuli_params, baseimage, stimuli)` -> matrix or vector; row
  selection, the empty check, and the 4096 -> 4092 pre-0.3.0 truncation in both its
  multi-trial and single-trial forms.

**`R/rdata.R`** (extend) — `loadStimulusParams(rdata)` -> list(p, base_faces,
stimuli_params, img_size). Runs `load()` in *its own* frame, validates the four required
objects, converts `s` -> `p`, and calls `rdataWriterNote()` on that frame.

**`R/ci-compute.R`** (new) — `computeParticipantCIs(...)` -> list(ci, pid_cis). The
`foreach` block, including the individual-CI save branch.

**`R/zmap-compute.R`** (new) — `computeZmapQuick(ci, sigma, threshold, img_size)` and
`computeZmapTTest(ci, params, responses, p, pid_cis, img_size, n_cores)`.

`generateCI()`'s body lands at roughly 100 lines: validate -> load -> select -> compute ->
present -> return, each step one call.

## What I verified, and how

- **Baseline is green.** `devtools::install()` then `testthat::test_local()` on `3ed6349`,
  R 4.5.2, rcicr 1.2.3.9000. No failures, no errors.
- **`pid.cis` is a real cross-branch dependency.** `R/generateCI.R:392` (`noiseimages <-
  pid.cis`) reads a variable created inside the *participants* branch of the CI computation,
  60 lines earlier and in a different `if`. Any extraction has to make that an explicit
  argument rather than a shared frame.
- **That line is not covered by any test.** Checked every `test_that` block in
  `tests/testthat/` for one setting both `participants` and `zmapmethod = 't.test'`. The
  single file mentioning both, `test-parallel-progress.R:46`, uses them in *separate*
  `generateCI()` calls (`participant_ci()` and `zmap_ci()`), so the combined branch never
  runs.
- **The combined path does work today, and I have its numbers.** 32px, 12 trials, 3
  participants, seed 1, `nscales = 1`: `zmap[1:3, 1]` is
  `-0.0379431539017, -0.2251794128693, -0.3224847810780`; `sign(zmap) == sign(ci)`; no NAs.
  Those three cells are a smoke check that the path runs, **not** the pin: stage 0 stores the
  whole matrix, for the reason given there.

## Staging

One PR per stage, squashed, so a red gate bisects to a stage rather than to a 300-line
rewrite. Ordered by rising risk:

- **Stage 0 — cover the uncovered branch.** Add the `participants` + `t.test` test above.
  No `R/` change. This is the safety net stages 3 and 4 rest on, so it goes first.

  **It pins the whole z-map, not a sample of it.** Three cells plus signs would pass while
  every other magnitude moved — and this is the one branch where stage 3 could change
  researchers' numbers with nothing else watching. The test compares the full 32x32 matrix
  element-wise against a reference captured from the pre-refactor tree and stored as an
  `.rds` fixture, alongside the whole-matrix summaries that
  `test-regression-baseline.R:137-150` uses for the same job (`sum(abs())`, `sd`, `min`,
  `max`, spot cells). The summaries alone are the house idiom and cover every cell, but they
  are not an element-wise comparison; the fixture makes it one, and `fixtures/` already holds
  references of this kind (`zmap-raster-reference-input.rds`).
- **Stage 1 — `R/ci-inputs.R`.** The four pure helpers, with direct unit tests of each
  helper's contract. The 4096 -> 4092 truncation is **already pinned end to end** by
  `test-fixed-bugs.R:196-233`, which exercises both the single-trial and multi-trial branches
  and compares the single-trial CI against the exact expected noise image. That regression
  stays as it is and is the real guard; the helper-level test asserts
  `selectStimulusParams()`'s own contract on the inputs it can be handed directly, and must
  not be written as if the behaviour were newly covered.
- **Stage 2 — `loadStimulusParams()`.** Removes `captureArgs()`/`list2env()` from
  `generateCI()` entirely: with `load()` confined to a helper's frame there is no user
  argument in scope for a `.Rdata` field to clobber. `captureArgs()` itself stays — its
  other callers are `computeInfoVal2IFC()` (`R/computeInfoVal2IFC.R`) and
  `computeCumulativeCICorrelation()` (`R/computeCumulativeCICorrelation.R`), and
  `test-load-argument-guards.R:16-90` pins both against exactly this hazard.
  `generateReferenceDistribution2IFC()` is *not* one of them: it hand-rolls an explicit
  `.args` list at `R/generateReferenceDistribution.R:78-89` and restores each field by name.
- **Stage 3 — `computeParticipantCIs()`.**
- **Stage 4 — the two z-map helpers.** Carries its own parallel comparison, for a gap stage 0
  does not close: the `t.test` z-map's `foreach` loop runs only on the *participant-free*
  path, because with `participants` set the branch short-circuits to `noiseimages <- pid.cis`.
  Every existing exercise of that loop forces one core — `test-regression-baseline.R:105-119`
  and `tools/compare-harness.R:295-308` both pass `n_cores = 1`, and
  `test-parallel-progress.R` runs two but asserts only that a progress bar appeared. So stage
  4 must add an `n_cores = 1` vs `n_cores = 2` **z-map value** equality test on the
  participant-free `t.test` path, or it could change parallel z-map numbers with nothing to
  catch it.

This branch carries stage 0 and stage 1, then deletes this file.

**Stages 2-4 each change `R/` behaviour, so each begins with its own plan commit on its own
branch and its own plan-only draft review** — `AGENTS.md` "Git and merge strategy" and
`CONTRIBUTING.md` "Plan first, in the same pull request" apply per PR, not once per issue.
Those plans are short: this one's analysis does not have to be repeated, only the stage's own
target shape, what it verified, and how it will be shown correct. Copying the staging list to
issue #184 records *what is left* and does not substitute for that procedure.

**The staging list still goes onto issue #184 before this PR is squashed**, since this file
does not survive the branch and the tracker is where status lives.

## The step most likely to fail: stage 3

`foreach::getexports()` scans the `%dopar%` body for free variables and `get()`s each one
from the enclosing frame. This has already broken once here — #235, where `targetpath` was
unbound in a branch that could not execute, and every `participants` call with the default
`n_cores > 1` aborted anyway. Line 136-146's comment exists because of it.

Moving the loop body into a package function changes the free-variable set at the loop site.
Today `targetpath`, `mask`, `img_size`, `individual_scaling`, `individual_scaling_constant`,
`baseimage`, `antiCI`, `participants` and `save_individual_cis` are all free in that body and
are exported to workers by name. Afterwards they become `computeParticipantCIs()` arguments,
and the body's free set is drawn from *its* frame instead. That is the improvement, and it is
exactly where an unbound name would reappear.

So stage 3 is not merged on a serial run. It must show:

- `n_cores = 1` and `n_cores = 2` returning identical CIs. `test-parallel-equivalence.R:12`
  has the shape to copy but **does not currently cover this**: it varies `ncores` for
  *stimulus generation*, and its two `generateCI()` calls neither set `n_cores` nor pass
  `participants`, so the participant loop never runs and both calls use the same core count.
  This is a new comparison, not an extension of an existing one;
- the `save_individual_cis` branch writing its PNGs under `n_cores = 2`, which is the
  configuration #235 broke;
- `targetpath` absent entirely (`save_as_png = FALSE`, `n_cores = 2`), which is the
  configuration that made #235 visible.

Stage 4 inherits a smaller version of the same risk in the `t.test` z-map loop.

## Open for review

1. **Five PRs for one issue, or fewer?** Stages 0+1 are genuinely low-risk and could carry
   stage 2. I have kept stage 2 separate because it deletes the `captureArgs()` call that a
   documented bug (the `sigma` clobber) exists to prevent, and that deserves its own diff.
2. **File names.** `ci-inputs.R` / `ci-compute.R` / `zmap-compute.R` follow `rdata.R` and
   `parallel.R`. An alternative is one `generateCI-internals.R`, which keeps the grep locality
   but recreates a 300-line file.
3. **Should stage 2 keep `captureArgs()` in `generateCI()` as belt-and-braces?** My view is
   no: two mechanisms guarding the same hazard is how one of them rots unnoticed. But it is a
   deliberate removal of a guard, so it is worth a second opinion.
