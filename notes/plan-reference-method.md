# Plan: let the user choose how the InfoVal reference is computed (#354, revised)

This revises the plan already reviewed on #358. It replaces the automatic routing between the Gram
and rendered calculations with an argument the user controls. Nothing else in that plan changes.

## Why

Automatic routing on trials against pixels can pick the Gram route where it needs more memory than
before. Codex's example on #358 is 64px with 4,000 trials: the noise matrix (131 MB) and `G`
(128 MB) are held together. Routing on estimated peak memory instead would keep small sets on the
old calculation, and move numbers without the user seeing why. Making the method an argument gives
two things automatic routing cannot:

- **exact reproduction:** an old script adds one argument and gets its published InfoVal back, bit
  for bit;
- **no silent switch:** the package recommends, and the user decides.

## The argument

`reference_method = c("gram", "images")`, appended last to `computeInfoVal2IFC()` and
`generateReferenceDistribution2IFC()`, after `reference_stimuli`, and validated with `match.arg()`.

- **`"gram"` (default)** computes the norms from the stimulus Gram matrix, whatever the sizes. `G`
  is still built whichever way holds less, from the sparse-rendered noise or the parameters'
  cross-product. That choice is internal and changes nothing beyond rounding.
- **`"images"`** is the calculation of rcicr 1.5.0 and earlier, bit for bit: the noise images
  rendered in parallel over `ncores`, then one pixels-sized product per draw.

The names describe what each computes from. `"images"` is not called `"legacy"`, because with more
stimuli than pixels it is the better choice, not merely the old one. It is not called `"exact"`,
because both are exact arithmetic and differ only in rounding.

The routing threshold and its tests go. The build choice and its boundary test stay.

## The notice

A `message()`, only when a reference is actually computed (never on a cache hit, which produces no
new number), and silenced by `suppressMessages()`:

- **With `"gram"`**, always:
  > InfoVal reference computed with reference_method = "gram". References from rcicr 1.5.0 and
  > earlier used "images"; pass reference_method = "images" to reproduce them bit for bit.
- **When the estimate favours the other method**, one sentence is added to the same message, giving
  both estimates and the argument that needs less. With `"images"` there is otherwise no message.

It never warns or errors, and never changes the number.

## The estimate

Peak working memory in doubles, counting the dominant allocations. `P` is pixels (width times
height), `n` the selected stimuli, `K` the parameters, `L` the basis layers, and `k` the response
block, `min(1000, iter)`.

- `"images"`: `P * n`, the noise matrix. This is a lower bound: the parallel `cbind` currently
  combines it in copies.
- `"gram"`, as the sum of three terms:
  - the sparse basis: `1.5 * P * L`, 8 bytes per value plus 4 per index;
  - the larger of the build (`P * n + n^2` rendered, or `K^2 + n * K + n^2` cross-product, whichever
    the build picks) and the simulation (`n^2 + 2 * n * k`).

The implementation PR checks these against peak memory measured with `gc()` at a few sizes: 64px
with 4,000 stimuli, the 512px default, and 512px with 50 stimuli. The recommendation must point the
same way as the measured peaks. If an estimate misjudges one of them, the PR corrects the formula
and says so.

## Stored references

A reference already stored in the file is reused whatever `reference_method` says. Otherwise the
new default would silently regenerate every marked reference from before. To rebuild one exactly the
old way, pass `force_gen_ref_dist = TRUE, reference_method = "images"`.

Each newly stored reference records the method it was computed with, as append-only fields:
`reference_norms_method` for the shared reference, and `method` inside each `reference_norms_by_base`
and `reference_norms_by_stimuli` entry. `generateReferenceDistribution2IFC()` re-saves its frame, so
`reference_method` is removed before `load()`, exactly as `reference_stimuli` is, and a file's own
object of that name survives.

## Docs

- **`NEWS.md`.** The Reproducibility-impact entry is rewritten around the argument: `"gram"` is the
  default, differs by rounding, and `reference_method = "images"` reproduces earlier results bit for
  bit. The performance entry is unchanged.
- **`DECISIONS.md`.** The generalised entry records why this is an argument rather than automatic,
  within its word budget.
- **The help pages** document both values, the notice and the estimate.
- **The analysis.** Its text about routing is replaced by a note that the package lets the user
  choose. Its measurements stand, since it already measures both methods. One re-knit.

## Tests

- `"images"` is `identical()` to the rendered arithmetic, including with more stimuli than pixels
  and through the subset and independent-base paths.
- `"gram"` meets the existing 1e-12 parity on every basis layout, and also with more stimuli than
  pixels, which it now also handles.
- The notice:
  - appears with `"gram"` on computation, but not on a cache hit;
  - does not appear with `"images"` when the estimate agrees;
  - gains the recommendation sentence on each side of an estimate boundary, via a pure estimate
    function tested directly.
- A stored reference is reused under either method. A new one records its method, and a planted
  `reference_method` object survives the save.
- Mutations that must fail a test: ignoring the argument; dropping the cache-hit exemption from the
  notice; flipping the estimate comparison.

## Risk

The estimate is the step most likely to be wrong. It is a model of R's allocations, not a
measurement, so it is validated against `gc()` peaks before anything recommends from it. A wrong
recommendation costs memory or speed, never a number.
