# Plan — #301: reference generation re-reads the base images

## The defect

`generateReferenceDistribution.R:165` rebuilds the stimuli by calling `generateStimuli2IFC()`
with the saved `base_face_files`, so it reopens every original image. The `.Rdata` already holds
everything the calculation needs — `stimuli_params`, the noise basis `p`, `seed`, `n_trials` —
and the reference is a function of noise alone. Two failures follow:

- **A moved or archived experiment cannot produce a reference**, because the paths no longer
  resolve.
- **A uniform base generated with `maximize_baseimage_contrast = FALSE` cannot either**: the
  rebuild does not forward that argument, so contrast maximization is re-applied to an image it
  rejects.

## What was measured, and how

At `1345b08`, R 4.3.3, package installed from the working tree.

**All three committed legacy fixtures already fail**, with exactly this issue's error:

```
Base image "base" does not exist: /tmp/rcicr-legacy-120235/base.png
```

Their `base_face_files` point at temp directories on whichever machine wrote them. So the
archived-experiment case is not hypothetical — the repository's own fixtures cannot be scored
today, while carrying every field the reference needs.

**A file whose image still resolves gives the values this change must preserve.** Modern file, 4
trials, `nscales = 1`, `seed = 1`, `iter = 6`: norms `2.334, 1.847, 1.738, 1.548, 1.741, …`, and
the RNG state the call leaves behind, sampled as the next three `runif()` draws after a
`set.seed(99)` and one call: `0.913, 0.294, 0.459`.

**The RNG draw count can come from the saved parameters.** The replication below needs to consume
what `generateStimuli2IFC()` consumed. Three candidate sources agree everywhere:

| fixture | saved `nscales` | `ncol(stimuli_params)` | nparams from `nscales` | `max(patchIdx)` |
|---|---|---|---|---|
| legacy 1.0.1 | absent | 4092 | 4092 | 4092 |
| legacy 1.0.1 gabor | absent | 4092 | 4092 | 4092 |
| legacy 1.1.0 | 1 | 12 | 12 | 12 |

`ncol(stimuli_params)` is the direct answer and needs no `nscales`, which the 1.0.1 files do not
carry — today's code assumes 5 for them. Where the two would disagree (a file whose `nscales` was
non-default and unsaved, issue #81's case) `ncol` is the correct one, and no such file can
deviate numerically anyway because files of that vintage error out today.

## The change

**#305 already built this.** `generateBaseReference()` computes reference noise from
`stimuli_params` and the saved basis via `independentReferenceNoise()` — reading no image — and
then replicates the RNG the old path consumed: `set.seed(seed)`, one `runif(nparams)` per trial.
The shared-parameter path at `:165` is the only caller still rebuilding stimuli, so this makes it
do what the independent path does, and the two collapse into one.

1. Replace the `generateStimuli2IFC()` call with noise built from the saved parameters and basis,
   reusing `independentReferenceNoise()` rather than adding a second copy: it already resolves
   `p` versus the older `s`, validates `n_trials`, parallelises and drives the progress bar.
2. Replicate the RNG consumption exactly where the rebuild used to sit, using
   `ncol(stimuli_params)` as the draw count, so the simulated responses that follow land on the
   same stream. The default-seed guarantee documented in `?generateReferenceDistribution2IFC` is
   the thing this protects.
3. **The `nscales`, `noise_type` and `sigma` warnings go.** They exist because the basis was
   rebuilt and would otherwise be rebuilt wrong; once the saved basis is used, all three fields
   are unused on this path and warning that "the resulting infoVal will be wrong" would be false.
   This is the one part that is a judgement call rather than a consequence — it is included
   because leaving a warning that is no longer true is worse than removing it.

## Expected impact

- **No numeric change for any file that works today.** The saved basis is the one the stimuli
  were built from, so rebuilding it reproduces it; the norms and the RNG state above must come
  back identical.
- **Files that error today start working.** They have no previous value to preserve, which is
  also why the `ncol` choice above carries no reproducibility risk.
- **Release gate: no new `EXPECTED` entry, and none stops firing.** The battery runs four InfoVal
  configurations, three of them on this exact path, including non-default `nscales` and `gabor`
  with non-default `sigma`. The two listed deviations are about *v1.0.1's* behaviour and this
  touches only the current side, so they must still fire. A new deviation means the basis is not
  being reproduced; a stale entry means one stopped, and both fail the gate.

## Tests

New `tests/testthat/test-reference-from-saved-noise.R`:

- With the source image deleted, a reference is produced, and it **equals** the one from the same
  file with the image intact — the point is not that it runs but that it is the same number.
- A uniform base written with `maximize_baseimage_contrast = FALSE` produces a reference.
- The modern-file norms above are pinned, and the RNG state after the call is pinned, both
  against the values measured on the unfixed tree.
- Each legacy fixture produces finite norms, on a copy so `save_rdata` cannot touch the original.
- `computeInfoVal2IFC()` on a file with a moved image returns the same value as with it in place.
- Each will be checked to fail without the change with `git stash push -- R/`.

## Most likely to fail

The RNG replication. It is the whole of the "no numeric change" claim, and a wrong draw count is
invisible until the norms move — which is why the modern-file norms and the post-call RNG state
are both pinned rather than assumed, and why the release gate is the check that settles it.

## NEWS.md and DECISIONS.md

`NEWS.md`: a bug-fix entry — this only ever produced errors, so it sits below the behaviour
changes. It says an archived or moved stimulus set can now be scored, and that no existing value
changes.

`DECISIONS.md` is at exactly 5200 of 5200 words, and this needs no new entry: it makes an
existing one stale, and rewriting it in place is the whole job. "`set.seed()` in
`generateStimuli2IFC()` is load-bearing far beyond stimulus generation" opens
"`generateReferenceDistribution2IFC()` rebuilds the stimuli through it", which stops being true
here. The guarantee it protects does not change — the reference still consumes the stimulus
seed's stream before its response draws — but after this change it does so by replicating that
consumption rather than by re-running the generator, and *that* is the thing a future reader
must not "simplify" away. Rewriting the entry to say so should be close to word-neutral; any
shortfall comes out of the same entry.
