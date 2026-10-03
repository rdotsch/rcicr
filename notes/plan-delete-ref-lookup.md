# Plan: delete `ref_lookup` and the three imports only it uses (#392)

Decided on #392: delete. DECISIONS.md said to repopulate or delete, never half; this does all of the deleting.

## Why nothing returned changes

`ref_lookup` has had no rows since 2018, so `filter()` returns zero rows, `count() == 1` and `count() > 0` are both false, `ref_median` is never bound, and every call takes the `resolveReferenceNorms()` branch. That branch is all that remains. The interactive prompt sat behind `count() > 0`, so it could never run.

## The change

1. **`sharedReference()`**: delete the `tribble()`, both `filter()` blocks, the `summarise()`, the `yesno::yesno()` prompt and its non-interactive message, and the `exists("ref_median")` guard. Keep `requireStimulusSeed()`, which sits inside the same `if (!force_gen_ref_dist && !exists("reference_norms"))` block. It moves to that condition on its own, so a seedless file without a stored reference still stops with the same error *before* any simulation, as now. Drop the `ref_seed`/`ref_img_size`/`ref_n_trials` placeholders, which only silenced R CMD check notes about the table.
2. **Imports**:
   - remove `dplyr` and `yesno` from `DESCRIPTION` Imports, and the `@importFrom dplyr ...`, `@importFrom tibble tribble` and `@import yesno` tags (NAMESPACE regenerated);
   - **move `tibble` to Suggests**: three test files build tibbles to check tibble input, a supported use, and they get `skip_if_not_installed("tibble")` as CONTRIBUTING.md requires for a Suggests package.
3. **Tests**: the first test in `test-computeInfoVal2IFC.R` mocks `yesno()` "in case ref_lookup gets populated again". It keeps its assertion (a pre-seeded reference is used) and drops the mock. The `helper-fixtures.R` comment drops its mention of `yesno()`.
4. **Docs**:
   - AGENTS.md: the `ref_lookup` convention bullet goes.
   - DECISIONS.md: the "Repopulating `ref_lookup`" entry is replaced by one saying it was deleted and why. Its key (seed, `img_size`, `n_trials`, `iter`) cannot describe a reference that also depends on `nscales`, `noise_type`, `sigma`, `reference_stimuli`, the mask, `rng_kind` and `reference_method`.
   - `.lintr`: the pipe comment says `%>%` comes from an import "regardless"; after this the package has no pipe at all, so the comment is updated and the linter kept to stop a mix creeping in.
   - NEWS.md: "Performance and dependencies": dplyr and yesno are no longer imported, and tibble is only suggested.
5. **Container setup stays**: `tools/setup-container-r.sh` still installs the `yesno` stub and dplyr/tibble, because the gate's reference versions (v1.0.1 through 1.5.0) import them. Its comment and `CLAUDE.md`'s `yesno` bullet now say that is the only reason.

## Verification, and the step most likely to fail

- Full release gate against `main` (`--ref="$(git rev-parse origin/main)"`): every InfoVal key must show 0 deviations.
- Full test suite; `R CMD check` in CI (an undeclared or unused import fails there).
- **Most likely to fail:** step 1's guard. Today `requireStimulusSeed()` runs only when there is no stored reference and `force_gen_ref_dist` is FALSE. Moving it must keep exactly that condition. The existing seedless-file tests (#334) cover both the stored and the unstored case, and must pass unchanged.
