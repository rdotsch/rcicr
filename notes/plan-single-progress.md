# Plan: a single participant or a single trial must not fail (#351)

## What is wrong

Three progress bars are created with `min = 1` and a maximum that can be 1. `txtProgressBar()` requires
`min < max`, so all three stop with `must have 'max' > 'min'` before any work is done:

- `R/generateStimuli2IFC.R:201`: `max = n_trials`, so `generateStimuli2IFC(n_trials = 1)` fails;
- `R/ci-compute.R:26`: `max = npids`, so `generateCI(participants = )` fails with one participant;
- `R/zmap-compute.R:34`: `max = n_observations` for the `t.test` z-map.

Every other progress bar in `R/` already starts at 0.

## Measured with the three bars set to `min = 0`, on this branch before any other change

These are the calls, on a synthetic 32px PNG base labelled `face`, `nscales = 1`. `r6` is a six-trial
stimulus file; `r1` is a one-trial file from `generateStimuli2IFC(list(face = b), n_trials = 1,
img_size = 32, nscales = 1, stimulus_path = d, ncores = 1)`. The output is as captured:

```
generateStimuli2IFC(n_trials = 1)                         ok (its .Rdata and 2 PNGs written)
generateCI(1, 1, "face", r1)                              ok
generateReferenceDistribution2IFC(r1, iter = 5)           ok
generateCI(1:6, ..., participants = rep("a", 6)), n_cores = 1   ok
  the same, n_cores = 2                                   ok
generateCI(3, 1, ..., participants = "a")                 ERROR: task 1 failed - "incorrect number of dimensions"
generateCI(c(3, 4), c(1, -1), ..., participants = c("a", "b"))  ok
one participant's CI vs the pooled CI of the same 6 trials: identical() TRUE
t.test z-map, participants = rep("a", 6)                  ERROR: not enough 'x' observations
t.test z-map, generateCI(3, 1, ...) with no participants  ERROR: task 1 failed - "incorrect number of dimensions"
```

**Now works:**
- one generated trial: the stimuli, their PNGs, `generateCI()` and the reference;
- one participant with six trials, with `n_cores = 1` and `n_cores = 2`. Its CI is `identical()` to
  the pooled CI of the same trials, the direct calculation the issue asks for.

**Still fails:**
- one participant with one trial. `generateCI()` hands the participant path a parameter vector for a
  single selected trial, and `params[pid.rows, ]` indexes it as a matrix;
- the `t.test` z-map with one participant: `t.test()` fails at every pixel;
- the `t.test` z-map with one trial and no participants.

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
   that worked changes. After this guard, the z-map's own bar can no longer see a maximum of 1.
   Its `min = 0` is therefore for consistency with every other bar, and changes only what the bar
   displays.
4. **`NEWS.md`**, in the development version's "Bug fixes": the single-trial and single-participant
   calls that stopped, and the new `t.test` z-map message.

## Tests (`test-fixed-bugs.R`, per its convention: they assert the intended behaviour)

- `generateStimuli2IFC(n_trials = 1)` succeeds and writes its `.Rdata` and PNGs, and a CI and a
  reference can be computed from that file.
- One participant with several trials, `n_cores = 1` and `n_cores = 2`, gives a CI `identical()` to
  the pooled CI of the same trials.
- One participant with one trial gives the CI of that single trial, `identical()` to
  `generateCINoise()` on its parameters, with `n_cores = 1` and with `n_cores = 2`. The second is
  the case where the new one-row matrix is sent to workers.
- The `t.test` z-map with one participant, and with one trial and no participants, stops with the new
  message; with two participants it still runs.
- Mutations that must fail a test:
  - the stimulus bar or the participant bar back to `min = 1`;
  - dropping the one-row conversion.

  The z-map bar is not listed. Once the guard is in place it cannot fail, so a mutation there has
  nothing to detect.

## Risk

The participant path's parallel branch serialises `params` to workers. The one-row matrix must reach
them with its dimensions intact. The single-participant, single-trial test with `n_cores = 2` covers
it.
