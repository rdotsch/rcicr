# Plan: autoscale() by position, not by name (#375)

## The defect

`autoscale()` iterates `names(cis)` and indexes with `cis[[ciname]]`. With duplicate names only the first element of each name is read and written: the constant ignores the others and they come back without `$scaled`. With no names it fails with "dimensions must be a positive quantity".

## The change

1. Both loops iterate `seq_along(cis)`, so every element enters the constant and gets `$scaled` and its `scaling` attribute.
2. With `save_as_pngs = TRUE`, names are file names, so the call stops before anything is computed if any name is missing, empty or duplicated, naming them. One PNG silently overwriting another is the alternative.
3. An empty list stops with a clear message.

## Behaviour change

For a list with unique names (everything the batch functions return, and the gate's `autoscale` keys) the constant, `$scaled` and the PNGs are identical: same elements, same order. Only lists that were mishandled change.

## Tests

Duplicate names: constant covers both, both scaled. Unnamed with `save_as_pngs = FALSE`: works. Unnamed or duplicated with `save_as_pngs = TRUE`: errors naming the problem, writes nothing. Unique names: output `identical()` to the current implementation's (pinned against values computed on main).

## Step most likely to fail

None of the numbers should move; the gate run against `main` confirms it.
