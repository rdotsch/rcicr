# Plan: stage 4 of the `generateCI()` split — the two z-map helpers (#184)

Last of four PRs for #184. Stages 0–1 landed in #256, stage 2 in #257, stage 3 in #258; the
staging is on [#184](https://github.com/rdotsch/rcicr/issues/184#issuecomment-5309729168).

**Nothing about the public function changes.** Same signature, same argument meanings, same
numbers. `test-regression-baseline.R` and the release gate hold that claim to account.

## Target shape

`R/zmap-compute.R` (new), following the per-concern file convention `R/rdata.R`,
`R/parallel.R`, `R/ci-inputs.R` and `R/ci-compute.R` set:

```
computeZmapQuick(ci, sigma, threshold, img_size)                     -> zmap
computeZmapTTest(ci, params, responses, p, pid_cis, img_size, n_cores) -> zmap
```

`computeZmapQuick()` absorbs the blur / rescale / threshold block; `computeZmapTTest()` absorbs
the noise-image `foreach`, the `pid_cis` short-circuit, the per-pixel `t.test` and the
`sign(ci) * abs(qnorm(pmap / 2))` conversion.

## The cross-branch dependency becomes an argument

Today the `t.test` branch reads `pid.cis` — a variable created 60 lines earlier, inside the
*other* arm of the CI computation — via `if (all(is.na(participants))) { ...build... } else {
noiseimages <- pid.cis }`. That reach across branches is the one construct stage 0 exists to
pin, and stage 3 already made it an explicit return value.

Stage 4 finishes the job by making it an explicit **parameter**: `pid_cis` is `NULL` when the
single-shot path was taken, and the helper branches on `is.null(pid_cis)` rather than on
`all(is.na(participants))`. The two conditions are equivalent — `pid.cis` is created if and
only if `participants` was given — but the second states the actual dependency instead of
re-deriving it from an argument the helper would otherwise not need.

`generateCI()` therefore sets `pid.cis <- NULL` in the single-shot arm, so the name is always
bound before the z-map section reads it.

## The second `foreach`, and why the plan does not reuse the first

`computeZmapTTest()` contains the package's *other* `%dopar%` loop, over trials rather than
participants. The `getexports()` reasoning from stage 3 applies unchanged: every free variable
of that body (`weightedparameters`, `p`, `cl`, `pb`) must be built inside the helper or be one
of its formals.

It is **not** merged with `computeParticipantCIs()`'s loop despite the similar shape. They
iterate different things (trials vs participants), combine into different dimensions, and one
has an individual-CI save branch the other does not. Merging them would need a callback
parameterised over both, which puts a closure in the `%dopar%` body — precisely the kind of
free variable `getexports()` handles least predictably.

## The one new test

An `n_cores = 1` vs `n_cores = 2` **z-map value** comparison on the participant-free `t.test`
path — the only path that executes this loop. It does not exist today: every current exercise
of it forces one core (`test-regression-baseline.R:105-119` and
`tools/compare-harness.R:295-308` pass `n_cores = 1`), and `test-parallel-progress.R` runs two
but asserts only that a progress bar appeared.

Stage 3's lesson carries over directly and is the reason this is worth writing carefully:

- The comparison must be **path-specific-mutation testable**. A perturbation inside the shared
  helper reaches both core counts equally and cannot redden an equivalence test; only a
  mutation gated on `!is.null(cl)` can.
- Stage 0's note applies to any mutation of the *input*: the t statistic is invariant to a
  uniform rescaling, so `noiseimages * k` is a no-op by construction. Perturb the mean.

## What must stay green, unmodified

- **`test-generateCI-paths.R`** — stage 0's element-wise pin of the participants + `t.test`
  z-map against `fixtures/zmap-participants-ttest.rds`. This is the guard for the
  `pid_cis` short-circuit that stage 4 turns into a parameter, so it is the primary check.
- **`test-regression-baseline.R`** — the golden master, including its z-map values.
- **`test-fixed-bugs.R:288-325`** — the `sigma` clobber, which is a z-map value test.

## How I will show it correct

- Full suite, installed, plus `lintr::lint_package()`.
- Mutation checks, one per property, path-specific where the property is an equivalence:
  perturb the mean of `noiseimages` (stage 0's pin should redden); gate a perturbation on
  `!is.null(cl)` (the new comparison should redden).
- The release gate, which runs the actual v1.0.1 code. Note its battery is chosen by the
  *reference* version, so z-maps below 512px and undecorated z-maps are not in it — the
  in-repo z-map tests are what cover those.

## After this

`generateCI()`'s body should be roughly 100 lines: validate → load → select → compute →
present → return, each step one call. Closing #184 belongs to this PR.
