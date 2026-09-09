# Plan — #303: a zero-signal classification image scales to NaN

## The defect

Exactly cancelling responses give a valid all-zero raw CI. Three scaling paths then divide by a
zero range and return NaN for every pixel:

- `applyScaling(scaling = 'independent')`, `generateCI.R:363-372` — the scaling constant is
  `max(abs(range(ci)))`, so `(ci + 0) / (2 * 0)`.
- `applyScaling(scaling = 'matched')`, `generateCI.R:358-361` — divides by
  `max(ci) - min(ci)`.
- `autoscale()`, `autoscale.R:59-69` — the same constant across a list, zero when *every* CI in
  it is zero.

`$combined` inherits the NaN. The raw `$ci` is correct and finite: no net signal is a result, not
a failure.

## What was measured, and how

Locally, R 4.3.3, package installed from the working tree at `be3facb`, on a
`stimuli = c(1, 1, 2, 2)`, `responses = c(1, -1, 1, -1)` CI whose raw `$ci` is all zero.

| path | scaled result |
|---|---|
| `none` | `0` — unchanged, correct |
| `constant` | **`0.5`** |
| `matched` | **NaN**, 1024/1024 pixels |
| `independent` | **NaN**, 1024/1024 pixels |
| `autoscale()`, every CI zero | **NaN**, 1024/1024 |
| `autoscale()`, one zero and one signal CI | **`0.5`** for the zero CI, signal CI finite |
| masked zero CI, `independent` | 512 intended NA kept, **512 unmasked pixels NaN** |

Two of the five paths already answer this question, and both answer `0.5`. The fix makes the
other two agree rather than inventing a policy.

Also measured: a base image's range after the default contrast maximization is exactly `[0, 1]`.

## The change

A zero-range CI renders **neutral** rather than NaN, and is not an error — a batch with
`save_individual_cis = TRUE` would otherwise abort over one participant whose responses happen to
cancel.

1. `independent`: a zero constant yields `0.5`. That is the limit of `(ci + c) / (2c)` as
   `c -> 0`, and what `constant` scaling already returns for the same input.
2. `autoscale()`: same, when the constant across the whole list is zero.
3. `matched`: the midpoint of the base image's range, `min(base) + (max(base) - min(base)) / 2`.
   **Not a literal `0.5`**: this method's contract is to map the CI onto the base image's
   intensity range, and `0.5` can fall outside it — a base spanning `[0, 0.3]` would render a
   no-signal CI brighter than any pixel in the base. With the default contrast maximization the
   base range is `[0, 1]`, measured above, so this *is* `0.5` in the default case and differs
   only for a base that was not normalised.
4. Masked pixels keep their NA in all three: the neutral value fills only unmasked pixels.
5. A warning names the degenerate CI, so a silently uniform image is never mistaken for a
   rendering bug.

Unchanged: `none` and `constant`; every non-degenerate result; `autoscale()`'s policy of leaving
`$combined` alone; the all-NA error in `autoscale()`, which is a masking mistake rather than a
zero-signal result.

## Expected impact

- **No change to any CI with signal**: every edited branch is reached only when a range is
  exactly zero, which is currently NaN. Nothing that produces a number today produces a
  different number after.
- **Release gate: expected clean, no new `EXPECTED` entry.** The battery has no cancelling
  response set, so no configuration reaches these branches. A flag there means the guard is
  catching a non-degenerate case and the change is wrong.

## Tests

New `tests/testthat/test-zero-range-scaling.R` — a focused file, as #308 and #305 did, not
`test-fixed-bugs.R`, which is the modernization P0 inventory:

- A cancelling response set gives a raw `$ci` that is all zero and finite, and `$scaled` and
  `$combined` that are finite under `matched` and `independent`.
- `independent` and `autoscale()` give `0.5` on every unmasked pixel; `matched` gives the base's
  midrange, asserted against a base whose range is **not** `[0, 1]` so the two policies are told
  apart rather than coinciding.
- A half-masked zero CI keeps 512 NA and has no NaN in the other 512.
- `autoscale()` with one zero and one signal CI is unchanged, and `$combined` is still untouched.
- Non-degenerate scaling is pinned for all four methods, so the guard cannot alter a CI that has
  signal.
- Each will be checked to fail without the change with `git stash push -- R/`.

## Most likely to fail

The warning, on the individual-CI path. `applyScaling()` is called from inside `ci-compute.R`'s
`foreach` body, and a warning raised in a worker does not reliably reach the user the way a
serial one does. The assertion is written serial-only; if the parallel case proves it cannot
surface, the warning stays serial rather than being replaced by something noisier, and the
limitation is stated where the argument is documented.

## NEWS.md and DECISIONS.md

`NEWS.md`: under the existing **Behaviour changes** heading in the development section, below
#302's entry — this one only ever produced NaN, where #302 produced missing files.

`DECISIONS.md` gets the reasoning for neutral-not-error, and for `matched` differing from the
other two. It is at 5197 words against a 5200 budget, so something comes out in the same commit;
`AGENTS.md` requires the trim rather than the overflow.
