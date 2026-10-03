# Plan: generateCI() without the zmapbool rename or the duplicated return (#388)

Part of #384.

## The change

1. The z-map matrix gets its own name, `zmap_matrix`, so the `zmap` flag keeps its meaning and `zmapbool` and its comment go.
2. The return is built once: `out <- list(ci, scaled, base, combined)`, `out$zmap <- zmap_matrix` when a z-map was made, then one `structure(out, trial_design = design, scaling = scaling_record)`. Same element names and order.

## Proving nothing changed

- A new test asserts `identical()` between the result and a list built the old way, with and without a z-map, so the element order and attributes are pinned before the change.
- **Full** release gate against `main` (`Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"`): `0 expected deviations, 0 unexpected`.
- Full `testthat::test_local()` passing with **no existing test modified**.
- `lintr::lint_package()`: 0 lints.

## Step most likely to fail

Appending `$zmap` to a list that already has attributes: building the list first and adding the attributes last avoids it.
