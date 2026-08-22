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

**`cd04b28` broke every caller of `generateCI()` on 2016-10-28.** "Every" is exact, not
rhetorical: at that commit `generateCI()` had two direct callers, `batchGenerateCI()`
(`R/rcicr.R:557`) and `generateCI2IFC()` (`R/rcicr_2IFC.R:177`), both calling positionally and
both broken, with `batchGenerateCI2IFC()` broken through `generateCI2IFC()`.
`R/rcicr_simulations.R` had none. The 2IFC pair stayed broken until 2021-12-28;
`batchGenerateCI()` was caught and repaired eleven days later. Each called `generateCI()`
positionally:

```r
generateCI(stimuli, responses, baseimage, rdata, saveaspng, filename, targetpath, ...)
```

That was correct until `cd04b28` (2016-10-28) inserted `participants = NA` into `generateCI()`'s
fifth slot without updating the call. From then on `participants` received `saveaspng`. Since
`all(is.na(TRUE))` is `FALSE`, every call entered the per-participant branch with one "level",
reaching `txtProgressBar(min = 1, max = 1)` and stopping with `must have 'max' > 'min'`.
Confirmed at 0.4.0 for both 2IFC functions, with `saveaspng` `TRUE` and `FALSE`; it fails either
way. Fixed by `90fee07` (2021-12-28).

**`batchGenerateCI()` had the identical call, and so the identical failure, for eleven days.**
At `8e44cfb^` it called `generateCI()` positionally in exactly the same shape. Reproducing R's
argument matching against that commit's real formals binds `participants` to `saveaspng` —
`TRUE` — giving `all(is.na(participants))` `FALSE`, `npids` 1, and `must have 'max' > 'min'`
from `txtProgressBar(min = 1, max = 1)`. `8e44cfb` ("Bugfixes to make batchGenerateCI work
again", 2016-11-08) repaired it by switching to named arguments with an explicit
`participants = NA` — the line that later write-ups cite as the reason `batchGenerateCI()` is
safe. That line is the *fix* for this defect, not an original design choice.

So no per-participant route ran at all between 2016-10-28 and 2016-11-08, and the advisory says
so rather than claiming `batchGenerateCI()` as an available alternative across the whole window.
Everything measured at 0.4.0 sits after `8e44cfb`, which is why `batchGenerateCI()` scored 0 of
12 there rather than failing.

**`cd04b28` is also the commit that implemented [issue #27](https://github.com/rdotsch/rcicr/issues/27)**,
the request for equally-weighted group CIs — the change that added `participants` in the first
place. Two days earlier, [issue #33](https://github.com/rdotsch/rcicr/issues/33) had been opened
asking to "check whether the functions that call `generateCI` still work". That check was never
run; all three callers broke on 2016-10-28, and #33 stayed open until 2026.

**This window is not the same as "0.4.0", at either end.** 0.4.0 ran to 2021-09-23 and 0.4.1
took over from there, but the misalignment lasted until 2021-12-28 — so 0.4.1 installs taken
from `development` between those two dates failed in exactly the same way. At the other end,
0.4.0 appeared on 2016-10-26 (`271751f`) and the misalignment only arrived on 2016-10-28, so
0.4.0's first two days are clean: `cd04b28^` is `a5921dc`, still 2016-10-26, and its positional
call still lines up. The boundary that matters is the date, not the version string, and both
the advisory and [`individual-ci-mislabelling.md`](individual-ci-mislabelling.md) state it as a
date.

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
"not affected" holds. It changes what to tell anyone on a `development` install predating
2021-12-28, and the answer has three windows rather than one: through 2016-10-27 all three
routes ran, from 2016-10-28 to 2016-11-07 none of them did, and from 2016-11-08 to 2021-12-27
`batchGenerateCI()` is the only one that both ran and named its output correctly. The advisory
says so.

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

## Do not date anything from an issue's close

**Issue close dates in this repository do not track when the code landed**, because the 2026
migration from `BACKLOG.md` to GitHub Issues closed a long tail of requests retroactively.

The trap is [issue #43](https://github.com/rdotsch/rcicr/issues/43), which asked for
`save_individual_cis` in the first place. It was opened 2016-11-08 and closed 2026-08-08 with a
comment beginning "Done." — phrasing that reads as a fresh implementation. Nothing was
implemented that day. The feature shipped in 2017: `a7b265b` (2017-07-31) added the argument,
`8a74974` (2017-08-15) made the saving work and introduced the mislabelling, and `d969ed2`
(2017-08-17) added `individual_scaling` and `individual_scaling_constant`. The only nearby 2026
change is `a2ccbd5` (2026-08-07), which made `targetpath` required in answer to CRAN's review.

This matters for the advisory specifically. A researcher who finds #43, reads a 2026 close on a
feature request, and concludes the option is two weeks old would then dismiss their own 2018 to
2023 output as impossible to affect — the exact wrong conclusion, and one nothing else on the
page would correct. A comment on #43 now gives the real dates and links the advisory.

The general rule: take dates from `git log`, and prefer `git log -S <the code itself>` over a
commit subject or an issue, both of which describe intent rather than what changed.
