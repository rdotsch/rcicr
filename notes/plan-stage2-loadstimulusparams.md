# Plan: stage 2 of the `generateCI()` split — `loadStimulusParams()` (#184)

Second of four PRs for #184. Stages 0 and 1 landed in #256; the staging for all five stages
is on [#184](https://github.com/rdotsch/rcicr/issues/184#issuecomment-5309729168). This plan
covers stage 2 only and does not repeat #256's analysis.

**Nothing about the public function changes.** Same signature, same argument meanings, same
numbers. `test-regression-baseline.R` and the release gate hold that claim to account.

## The problem this stage removes

`generateCI()` calls `load(rdata)` in its own frame (`R/generateCI.R:167`). `load()` assigns
straight into that frame, so any object in the `.Rdata` silently overwrites an argument of
the same name. This is not hypothetical: since 1.1.0 the file carries the *noise* `sigma`
(25 by default), which replaced the z-map blur `sigma` of 3, so every z-map built from a
1.1.0 stimulus set was smoothed with the wrong constant and the caller's `sigma` was ignored
— found by `tools/compare-release-output.R` on 2026-07-28, pinned by
`test-fixed-bugs.R:288-325`.

The current defence is a snapshot-and-restore: `.args <- captureArgs(environment())` before
the `load()`, `list2env(.args, envir = environment())` after it. That works, but it is a
*guard around* the hazard rather than removal of it, and it has to be maintained — `AGENTS.md`
and `CONTRIBUTING.md` both carry a standing rule to check every new argument name against
every saved field, precisely because the exposure is structural.

## Target shape

`loadStimulusParams(rdata)` in `R/rdata.R`, beside `captureArgs()` and `rdataWriterNote()`
which are already there. It runs `load()` in *its own* frame — a frame containing no user
argument — validates, and returns what `generateCI()` actually uses:

```
loadStimulusParams(rdata) -> list(p, base_faces, stimuli_params, img_size)
```

It absorbs `R/generateCI.R:165-192`: the `captureArgs()`/`load()`/`list2env()` trio, the four
`exists()` checks with their `rdataWriterNote()` calls, and the pre-0.3.3 `s` -> `p`
conversion. `rdataWriterNote()` already takes an environment argument, so it is called on the
helper's frame with no change to it.

At the call site those 28 lines become:

```r
loaded <- loadStimulusParams(rdata)
p <- loaded$p
base_faces <- loaded$base_faces
stimuli_params <- loaded$stimuli_params
img_size <- loaded$img_size
```

After which `generateCI()` has no `load()` in its frame, so **`captureArgs()` and
`list2env()` are deleted from it** — there is no longer an argument in scope for an `.Rdata`
field to clobber. The hazard is gone structurally rather than guarded.

`captureArgs()` itself **stays in `R/rdata.R`, unchanged**. Its remaining callers are
`computeInfoVal2IFC()` and `computeCumulativeCICorrelation()`, both of which still `load()`
in their own frames; `test-load-argument-guards.R:16-90` pins both. (`generateReferenceDistribution2IFC()`
hand-rolls its own explicit `.args` list at `R/generateReferenceDistribution.R:78-89` and has
never called `captureArgs()` — noted because #256's plan got this wrong once.)

`R/zzz.R`'s `globalVariables()` list is **not** touched. The helper will reach the loaded
objects through an explicit environment (`get()` / `exists()` on its own frame) rather than as
bare symbols, so it adds no new undefined-global references — and the existing entries stay
needed regardless, because `computeInfoVal2IFC()` and `computeCumulativeCICorrelation()` still
`load()` the same names into their own frames.

## What must stay green, unmodified

Both are existing tests. Neither may be edited to accommodate this change; if either needs
touching, the extraction is wrong.

- **`test-fixed-bugs.R:288-325`** — the `sigma` clobber. Today it passes because
  `captureArgs()` restores `sigma`; after this stage it passes because `sigma` was never in
  the frame `load()` wrote into. **The test does not change, but what it proves does** — it
  goes from testing the guard to testing the structure. That is the whole point of the stage,
  and it is why this test staying byte-identical and green is the primary success criterion.
- **`test-fixed-bugs.R:422-439`** — `batchGenerateCI()` forwarding a *missing* `targetpath`.
  This is the forwarded-missing-promise case `captureArgs()` was written for. The
  `targetpath <- if (missing(targetpath)) NULL else targetpath` binding at
  `R/generateCI.R:147` is what actually handles it and **stays**; only its ordering comment
  changes, because the "must stay above `captureArgs()`" constraint disappears with
  `captureArgs()`. The `missing()` checks at `R/generateCI.R:120-134` still must run before
  anything assigns to those names.

## How I will show it correct

- Full suite, installed (not just `load_all()`), plus `lintr::lint_package()`.
- **Mutation check on the new helper**, in the idiom stage 0 used: make
  `loadStimulusParams()` return the file's `sigma` as a fifth element and assign it into
  `generateCI()`'s frame, and confirm `test-fixed-bugs.R:288` goes red. If it stays green the
  test is not watching what this stage claims it watches.
- A direct unit test of `loadStimulusParams()`'s own contract: the four fields returned, the
  pre-0.3.3 `s` -> `p` conversion, and each of the four missing-object errors carrying its
  `rdataWriterNote()` suffix.
- `test-regression-baseline.R` and `test-load-argument-guards.R` green, untouched.

## Open question

Whether to fold stage 3 in. Against: stage 3 is the `foreach` extraction, the riskiest of the
five (`foreach::getexports()` scans the `%dopar%` body for free variables, and moving that
body changes the free set — the failure mode of #235). Keeping it alone means a revert of it
reverts nothing else. This plan assumes **not** folded.
