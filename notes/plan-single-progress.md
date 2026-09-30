# Plan: a single participant or a single trial must not fail (#351)

## What is wrong

Three progress bars are created with `min = 1` and a maximum that can be 1. `txtProgressBar()` requires
`min < max`, so all three stop with `must have 'max' > 'min'` before any work is done:

- `R/generateStimuli2IFC.R:201`: `max = n_trials`, so `generateStimuli2IFC(n_trials = 1)` fails;
- `R/ci-compute.R:26`: `max = npids`, so `generateCI(participants = )` fails with one participant;
- `R/zmap-compute.R:34`: `max = n_observations` for the `t.test` z-map.

Every other progress bar in `R/` already starts at 0.

## Measured with the three bars set to `min = 0`, on this branch before any other change

**Now works:**
- one generated trial: the stimuli, their PNGs, `generateCI()` and the reference;
- one participant with six trials, with `n_cores = 1` and `n_cores = 2`. Its CI is `identical()` to
  the pooled CI of the same trials, the direct calculation the issue asks for.

**Still fails:**
- one participant with one trial: `incorrect number of dimensions`. `generateCI()` hands the
  participant path a parameter vector for a single selected trial, and `params[pid.rows, ]` then
  indexes it as a matrix;
- the `t.test` z-map with one participant: `not enough 'x' observations`, from `t.test()` at every
  pixel;
- the `t.test` z-map with one trial and no participants: `incorrect number of dimensions`.

## Changes

1. **All three bars start at `min = 0`.** The progress ticks are unchanged (1 to n), so every bar
   now starts from 0% rather than from 1/n. That is display only.
2. **`computeParticipantCIs()` accepts a single-trial parameter vector** by making it a one-row
   matrix first. A participant with one trial then takes the same one-row path it already takes
   when another participant has more trials, which `generateCINoise()` computes identically for a
   vector and a one-row matrix. No number changes for any input that works today.
3. **The `t.test` z-map stops with a clear error when its stack has fewer than two images**, before
   any work: fewer than two participants with `participants` given, or fewer than two distinct
   stimuli without. The error names the cause and suggests `zmapmethod = "quick"`. A one-observation
   t-test is not a supported case, so this replaces two cryptic errors with one clear one; nothing
   that worked changes.

## Tests (`test-fixed-bugs.R`, per its convention: they assert the intended behaviour)

- `generateStimuli2IFC(n_trials = 1)` succeeds and writes its `.Rdata` and PNGs, and a CI and a
  reference can be computed from that file.
- One participant with several trials, `n_cores = 1` and `n_cores = 2`, gives a CI `identical()` to
  the pooled CI of the same trials.
- One participant with one trial gives the CI of that single trial, `identical()` to
  `generateCINoise()` on its parameters.
- The `t.test` z-map with one participant, and with one trial and no participants, stops with the new
  message; with two participants it still runs.
- Mutations that must fail a test: any of the three bars back to `min = 1`, and dropping the
  one-row conversion.

## Risk

The participant path's parallel branch serialises `params` to workers. The one-row matrix must reach
them with its dimensions intact. The parallel single-participant test covers it.
