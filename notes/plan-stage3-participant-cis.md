# Plan: stage 3 of the `generateCI()` split — `computeParticipantCIs()` (#184)

Third of four PRs for #184. Stages 0 and 1 landed in #256, stage 2 in #257; the staging is on
[#184](https://github.com/rdotsch/rcicr/issues/184#issuecomment-5309729168). This plan covers
stage 3 only.

**Nothing about the public function changes.** Same signature, same argument meanings, same
numbers. `test-regression-baseline.R` and the release gate hold that claim to account.

## Why this is the riskiest of the five

`foreach::getexports()` scans the `%dopar%` body for free variables and `get()`s each one from
the environment the loop is evaluated in. **Moving that body changes which environment is
scanned**, so a variable that resolves today can stop resolving — which is exactly the failure
mode of #235, where an absent-and-required `targetpath` aborted every participant call at the
default core count even though the branch referencing it could not run.

The loop is `R/generateCI.R:198-231` and its body `:203-230`. That body has **sixteen** free
variables besides the loop index. Three are created immediately above the loop (`pids`, `pb`,
`cl`). The other thirteen come from `generateCI()`'s frame: nine of its arguments
(`responses`, `participants`, `mask`, `save_individual_cis`, `targetpath`,
`individual_scaling`, `individual_scaling_constant`, `antiCI`, `baseimage`) and four computed
earlier in the body (`params`, `p`, `base`, `img_size`). After the move every one of them must
be a formal of the helper or created inside it — there is no enclosing `generateCI()` frame to
fall back on.

## Target shape

`computeParticipantCIs()` in a new `R/ci-compute.R`, following the per-concern file convention
`R/rdata.R`, `R/parallel.R` and `R/ci-inputs.R` set:

```
computeParticipantCIs(params, responses, participants, p, base, baseimage,
                      img_size, mask, n_cores, save_individual_cis, targetpath,
                      individual_scaling, individual_scaling_constant, antiCI)
  -> list(ci, pid_cis)
```

It absorbs `R/generateCI.R:184-239`: the `pids`/`npids` derivation, the progress bar, the
backend start and teardown, the `foreach` loop with its individual-CI save branch, the
`dim()<-` reshape, and the average across participants.

Fourteen arguments is a lot, and that is the point rather than a smell: each one is a free
variable of the loop body, and passing them explicitly is what makes the `getexports()`
lookup a checked thing instead of an accident of scoping.

**`pid_cis` must be returned, not just `ci`.** The `t.test` z-map short-circuits to
`noiseimages <- pid.cis` (`R/generateCI.R:308`), reaching across branches into a variable this
loop created. That is the cross-branch dependency stage 0 exists to pin, and it is why stage 3
cannot return the group CI alone.

The `cl <- NULL` after `parallel::stopCluster(cl)` stays: `on.exit(stopClusterSafely(cl))`
resolves `cl` at exit time, and blanking it is what stops a second teardown of an
already-stopped cluster. Inside the helper the `on.exit` now fires at *helper* exit rather than
at `generateCI()` exit, which is strictly earlier and still covers the loop it guards. The
z-map branch builds its own `cl` and `pb` further down and is unaffected.

## What must stay green, unmodified

- **`test-fixed-bugs.R:383-420`** — #235. Covers `targetpath` legitimately absent at
  `n_cores = 1` and `n_cores = 2`, and `save_individual_cis = TRUE` writing PNGs at
  `n_cores = 2`. This is precisely the `getexports()` hazard the move disturbs, so it is the
  primary guard and must not be touched.
- **`test-regression-baseline.R`** — the golden master.
- **`test-parallel-progress.R`** — the progress bar still ticks from the parent.

## The one new test

An `n_cores = 1` versus `n_cores = 2` **value** comparison on the participant path. This does
not exist today, and the gap is easy to miss because the nearest test looks like it closes it.
`test-parallel-equivalence.R:12` is titled "ncores = 1 and ncores = 2 produce identical stimuli
and CIs" and its `:46` comment says "check the CI as well, since `generateCI()` has its own
parallel loop" — but the two calls at `:47-48` set neither `n_cores` nor `participants`, so
they take the single-shot `generateCINoise()` path, which contains no `foreach` loop at all.
The `ncores` variation in that test applies only to `generateStimuli2IFC()`. So nothing
currently checks that the participant loop returns the same numbers on two cores as on one.

Worth fixing that comment in passing, since it is the reason the gap reads as covered.

It must compare **both** returned elements — the group `ci` and the per-participant `pid_cis`
that feed the z-map — element-wise, not summaries. Comparing the group CI alone would be
materially weaker than it looks: it is a mean across participants, so it is blind to any change
in their *order*, which is exactly what a `.combine` bug in a parallel loop produces. Guarded
with `skip_on_cran()` in the idiom of `test-fixed-bugs.R:403`, since it spawns a real cluster.

Neither loop in this package draws random numbers, so there is no worker RNG stream to
diverge; equality here should be exact rather than tolerance-based, and the test will assert it
that way. If that turns out to be false the assumption is wrong and worth knowing.

## How I will show it correct

- Full suite, installed (not just `load_all()`), plus `lintr::lint_package()`.
- **Mutation checks — one per test, and they cannot be the same mutation.** An earlier draft of
  this plan proposed a single perturbation inside the helper and claimed it would redden both
  the new comparison and stage 0's z-map pin. It would not: the helper is shared by both core
  counts, so an alteration there reaches the `n_cores = 1` and `n_cores = 2` runs *equally* and
  they stay element-wise identical. An equivalence test can only be reddened by a
  **path-specific** mutation.

  - *For stage 0's z-map pin*: perturb one participant's CI inside the helper. Stage 0's note
    applies — the t statistic is invariant to a uniform rescaling of its input, so this must
    shift the mean rather than scale it.
  - *For the new equivalence test*: gate the perturbation on the backend, `if (!is.null(cl))`,
    so only the parallel run is altered. Both returned elements should go red.
  - *And a second, sharper one*: iterate the parallel path in reverse order only. Because the
    group CI is a **mean across participants** it is invariant to their order, so this leaves
    `ci` identical and reddens `pid_cis` alone. That is the concrete demonstration of why the
    comparison must cover both returned elements rather than the group CI — and it is the
    mutation most likely to model a real `foreach` `.combine` ordering bug.
- The release gate (`tools/compare-release-output.R`), which runs the actual v1.0.1 code rather
  than values this repo computed for itself.

## Open question

Whether `computeParticipantCIs()` should take the fourteen arguments individually or a single
list. Individually, in this plan: a list would move the `getexports()` problem rather than
solve it, since the loop body would then reference `args$targetpath` and the free variable
becomes `args` — one opaque object whose contents are no longer checkable at the call site.
