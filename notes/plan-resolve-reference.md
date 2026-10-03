# Plan: resolveReferenceNorms(), decision apart from reporting (#390)

Part of #384. 68 lines, cyclomatic complexity 35.

## The change

- `useStoredReference(entry, rdata, baseimage, seedless, subset, masked)`: the stored-reference branch (the "Using reference distribution found" line, the seedless message, the return).
- `reportSimulatedReference(readonly, save_rdata, rdata, baseimage, response_seed, stale, norms, entry)`: the four outcome lines after a simulation.
- `resolveReferenceNorms()` keeps the decision itself (`forced`, `stale`, `automatic`, `readonly`, `save_rdata`) and the simulation, in the same order.

## Proving nothing changed

Output order and text are the behaviour here.

- Before the change, a test captures the complete `stdout` and `message` output of each path, and must pass on `main`: stored, stored and seedless, simulated and saved, simulated read-only, seeded, stale and rebuilt with different values.
- **Full** release gate against `main` (`Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`): `0 expected deviations, 0 unexpected`.
- Full `testthat::test_local()` passing with **no existing test modified**.
- `lintr::lint_package()`: 0 lints.

## Step most likely to fail

The read-only branch's label depends on `baseimage`; the capture test covers both a shared and a per-base file.
