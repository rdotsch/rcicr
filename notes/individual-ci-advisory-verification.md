# Independent re-verification of the individual-CI advisory

Written before the advisory was announced publicly, to confirm from the running code — not from
the earlier write-ups — that the bug is where [`individual-ci-mislabelling.md`](individual-ci-mislabelling.md)
says it is, and that the alternative routes really are safe.

Every claim below was measured, not inferred. The versions were installed into throwaway
libraries and executed; where a conclusion rests on source alone, it says so.

## How the images were identified

A filename assertion cannot see this bug — that is why the original test stayed green for nine
years. So identification had to run off pixels.

Twelve participants, `p1 … p12`, in collection order, each given a **disjoint** block of ten
stimuli and a constant `+1` response. Each participant's CI is then the mean noise of their own
ten stimuli: a unique signature, near-orthogonal to every other participant's.

Ground truth came from `generateCINoise()` called once per participant **outside** the loop
under test. Each written PNG was then matched back by inverting the scaling exactly — the
individual CIs use `independent` scaling, so `combined = (scaled + base) / 2` is linear in the
CI and `2 * png - base` recovers it up to 8-bit quantisation — and correlating against all
twelve ground truths.

The true owner scores **r = 1.000**; the runner-up never exceeds **0.23**. There is no judgement
call in the matching.

## What was confirmed

**The bug is real in 1.2.3**, the last affected release. `generateCI(participants = ...,
save_individual_cis = TRUE)` put **11 of 12** images under another participant's name. Only `p1`
was right, and only because it is first in both orderings. The permutation is exactly the
documented one: the file named `unique(participants)[i]` holds `sort(unique(participants))[i]`.

**Only the names were wrong.** The set of MD5 hashes of the twelve PNGs is *identical* between
1.2.3 and 1.3.0. No pixel moved; the twelve names were permuted among the same twelve images.

**The recovery recipe in the advisory works.** Applied to the 1.2.3 output, the two-step rename
restores all twelve to the names 1.3.0 gives them, byte for byte.

**The group CI is untouched**, on both versions: `max |group_ci - mean(per-participant truth)|`
= 6.9e-18, and `all.equal()` is `TRUE`.

**`generateCI()` is the only affected function.** Two independent lines agree:

- *Source.* Across 1.2.3's `R/`, exactly one site pairs a `factor()`-ordered index with an
  `unique()`-ordered name — `R/generateCI.R:269` — and exactly one site writes into
  `individual_cis/`.
- *Measurement.* At 1.2.3, `batchGenerateCI()`, `batchGenerateCI2IFC()` and `generateCI2IFC()`
  each got **0 of 12** wrong, and none created `individual_cis/`. So did `batchGenerateCI()`
  under its **default** `scaling = 'autoscale'`, which writes through `autoscale()` rather than
  `saveToImage()` — a distinct path, and the one most callers actually take.

**The affected range is continuous.** `git log -S` on the defective expression finds exactly one
introduction (`8a74974`, 2017-08-15) and no removal before the 2026-08-17 fix. The expression is
byte-identical at every tag from `v1.0.1` to `v1.2.3` and present at the 0.4.1 and 1.0.0
boundary commits. CRAN-era code contains no mention of the option at all.

**0.4.0, on both sides of 2017-08-15.** Installing it needed a `NAMESPACE` patch — `spatstat` no
longer exports `blur`/`as.im` and `purrr::rbernoulli` is gone — but `R/generateCI.R` was left
byte-identical.

- From `8a74974` on: **11 of 12** wrong, the same permutation as 1.2.3.
- Before it: the call **errors** with `argument "filename" is missing, with no default` and
  writes nothing. The pre-2017-08-15 claim is right, and stronger than stated — the feature was
  not merely correct there, it was non-functional. `saveToImage()` was called with three
  arguments against a five-argument signature.

Versions 1.0.1, 1.1.0, 1.2.0, 1.2.1 and 1.2.2 were **not** executed. They carry the identical
expression with no intervening change, which is taken as sufficient.

## What it corrected

The one substantive finding the earlier write-ups missed.

**`generateCI2IFC()` and `batchGenerateCI2IFC()` could not run at all between 2016-10-28 and
2021-12-28.** `generateCI2IFC()` called `generateCI()` positionally:

```r
generateCI(stimuli, responses, baseimage, rdata, saveaspng, filename, targetpath, ...)
```

That was correct until `cd04b28` (2016-10-28) inserted `participants = NA` into `generateCI()`'s
fifth slot without updating the call. From then on `participants` received `saveaspng`. Since
`all(is.na(TRUE))` is `FALSE`, every call entered the per-participant branch with one "level",
reaching `txtProgressBar(min = 1, max = 1)` and stopping with `must have 'max' > 'min'`.
Confirmed at 0.4.0 for both functions, with `saveaspng` `TRUE` and `FALSE`; it fails either way.
Fixed by `90fee07` (2021-12-28).

**This window is not the same as "0.4.0".** 0.4.0 ran to 2021-09-23 and 0.4.1 took over from
there, but the misalignment lasted until 2021-12-28 — so 0.4.1 installs taken from
`development` between those two dates failed in exactly the same way. The boundary that matters
is the date, not the version string, and both the advisory and
[`individual-ci-mislabelling.md`](individual-ci-mislabelling.md) state it as a date.

**It was reported at the time and sat open for four and a half years.**
[Issue #79](https://github.com/rdotsch/rcicr/issues/79), opened 2017-07-19, is a user hitting
precisely this: `Error in txtProgressBar(min = 1, max = npids, style = 3) : must have 'max' >
'min'`, with the remark "I ran this function with the same data in the old package and it
worked" — which is consistent with the call having been correctly aligned before `cd04b28`. It
was closed by `90fee07` on 2021-12-28.

**There is no pull request for the fix.** `90fee07` has a single parent, `77aa698`, so it was
pushed straight to `development`; the work was tracked by issue #79 rather than by a PR. This
predates the current squash-merge convention.

This does not weaken the advisory — a function that stops cannot emit a mislabelled file, so
"not affected" holds. It changes what to tell anyone on a pre-2021-12-28 `development` install:
`batchGenerateCI()` is the only per-participant route that both ran and named its output
correctly there, and the advisory now says so.

**Which branch this was broken on.** The breakage lived on `development`, not on the default
branch, and the default branch never served it:

- `master` was frozen at an R-Forge layout with no root `DESCRIPTION`, so a plain
  `install_github()` could not install anything from it at all.
- `90fee07` fixed the misalignment at 12:27 on 2021-12-28. `69e3933`, the merge that made
  `development` the new default, is nine minutes later at 12:36 and has `90fee07` as its first
  parent.

So the fix was already in at the moment the switch happened. `install_github('rdotsch/rcicr')`
never returned broken 2IFC functions; only `ref = 'development'` did, and only before
2021-12-28.
