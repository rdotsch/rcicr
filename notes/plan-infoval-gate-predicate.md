# Plan: give the InfoVal expectations a check predicate (#318)

## The gap

`tools/compare-release-output.R:58-62` states the rule against itself:

> `check`: optional predicate(ref_value, cur_value). Without one the key excuses *any*
> deviation in that output, including a later, unrelated regression in the same place.

Four InfoVal entries carry none — three scoped to `v1.0.1`
(`sinusoid-64-nscales3-infoval`, `gabor-64-sigma10-infoval`, `sinusoid-64-twobase-indep-infoval`)
and one to `v1.3.0` (`sinusoid-64-twobase-indep-infoval`), eight keys in total counting each
entry's `infoval` and `infoval_twice`. A fifth, the `v1.1.0` z-map blur entry, is a different
case; see "The dormant entries".

Nothing is wrong today. The corrected values are pinned at `1e-12` by
`tests/testthat/test-independent-reference.R:44-49`, against a value recomputed by
`reference_oracle()`, which rebuilds trial pixels from the saved parameters and basis rather
than calling the reference implementation. The test suite is a required check, so a regression
that moved these values fails there. This is about the *gate* being able to say it too, rather
than about missing verification.

## What the code allows

Verified rather than assumed:

- **A predicate runs in the gate's own process**, and recomputes expected values from source
  data on disk: `correct_alpha_base()` reads the base PNG through `alpha_bases()` and asserts
  *both* sides independently. It uses `png::readPNG` and base R only.
- **That process has no rcicr.** Each side is run by `system2()` on a separate `Rscript` with
  its own `libpath()` (`:385`, `:462`), the reference from a throwaway library and the working
  tree from its own. So a predicate cannot call `generateNoiseImage()`.
- **`compare-harness.R` runs inside each side, with that side's rcicr**, and already emits a
  derived scalar: `out$infoval_refresh_delta <- out$infoval_twice - out$infoval` (`:420`).
- **The harness must execute on the reference too** — `CLAUDE.md` warns a crash there aborts
  the whole gate. `generateNoiseImage(params, p)` exists in `v1.0.1` with the same signature,
  and `v1.0.1`'s `generateStimuli2IFC()` saves the basis as `p`, so an oracle built on them
  runs on both sides.
- **Entries are filtered to the resolved reference SHA** (`:349`), so an entry names exactly
  one reference and CI resolves only `v1.0.1` and the newest tag.

## The change

The oracle belongs in the harness, where rcicr and the saved `.Rdata` are both in hand, not in
the predicate, which would need `generateNoiseImage()` reimplemented in base R.

1. `compare-harness.R` gains an oracle beside the InfoVal call: rebuild each trial's pixels
   with `generateNoiseImage()` from the saved parameter matrix and basis, replay the same
   response draws, and compute `(norm(ci) - median(norms)) / mad(norms)` directly. Emit
   `out$infoval_oracle_delta` as the difference from `out$infoval`.
2. That output deviates by construction — nonzero on a reference that scores against the wrong
   base's null, zero on this tree — so it needs its own `EXPECTED` entry per reference, and
   that entry carries the predicate the existing ones cannot: the **current** side's delta is
   zero within tolerance. A regression moving the corrected InfoVal moves that delta off zero
   and is reported.

The existing `infoval` / `infoval_twice` entries then keep their reason and gain a cross
reference to the oracle output, rather than a predicate of their own: a predicate over those
two values alone cannot see the oracle, since `check` receives only its own key's pair.

## Most likely to fail

- **The reference side.** An oracle that errors on `v1.0.1` aborts the gate rather than failing
  one check. Two known wrinkles to handle rather than discover: `patchIdx` starting at 0 in
  older files, which `generateNoiseImage()` itself warns about, and the parameter-matrix width
  — the test oracle trims 4096 columns to 4092 before use. Measure by running the full battery
  against `v1.0.1` before touching any entry.
- **The configurations that do not deviate.** `sinusoid-64-infoval` has an InfoVal and no
  entry, so its new `infoval_oracle_delta` must agree across both sides — zero on each — or a
  clean configuration starts reporting a difference.
- **`infoval_twice`.** It exercises the cache-refresh path, so its oracle relationship is not
  obviously the same as `infoval`'s. Decide from a measured run whether one delta covers both
  or each needs its own.

## The dormant entries

The `v1.1.0` z-map blur entry (`defaults-512-sinusoid/{zmap_quick,zmap_plain}`,
`sinusoid-128-nscales3/zmap_plain`) and the `v1.2.3` individual-CI entry name references CI
never resolves, so nothing exercises them and no run can report them stale either. Tightening
them buys nothing while that holds. The useful question is whether they stay at all: keeping
them documents a real historical deviation and costs nothing; dropping them removes text that
cannot rot because it is never read. Settle it in this branch either way, and record the reason
in the entry or in its removal.

## Verification

- A deliberate mutation of the **corrected** InfoVal — not of the refresh path — makes the gate
  report a difference rather than absorb it. `--ref=$(git rev-parse origin/main)` is the sharp
  instrument here: against a branch no entry names, nothing is excused.
- The full battery against both references reports `0 unexpected` and none stale.
- `tests/testthat/` unchanged and still green: this touches `tools/` only.
