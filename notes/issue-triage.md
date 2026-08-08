# Issue triage — the 25 open issues not touched by the modernization pass

Prerequisite for the backlog's migration to GitHub Issues (done; see `backlog-migration.md`). Item 19 closed the
eight issues the modernization fixed directly; these 25 are the remainder. All date from
2016–2017, none is a bug report.

**Approved and posted 2026-08-08.** Every comment below went out as written: eight closes and
the correction on #87, which stays open. **One of those closes was wrong** — #9, see below.
It rested on a claim about `generateCI()`'s return value that is false; it was reopened the
same day with a public correction. Verdicts were verified against the working tree at `1.2.3.9000`; this file is the record
of what was posted and why, not a live view of the tracker.

## Summary

| verdict | issues | n |
|---|---|---|
| **Close — implemented, verified in code** | #43, #27, #13 | 3 |
| **Close — moot, the premise no longer holds** | #92, #52 | 2 |
| **Close — no actionable specification** | #68 | 1 |
| **Close — nothing in the package to attach it to** | #53 | 1 |
| ⚠️ **Closed in error, since reopened** | #9 | 1 |
| **Keep open, rescope** — most of it shipped | #87 | 1 |
| **Keep open — real, specified, unimplemented** | #71, #76, #85, #15, #35, #37, #38, #46, #47, #54, #22, #69 | 12 |
| **Keep open — thin, but the scope is clear** | #6, #7, #10 | 3 |
| **Keep open — aspirational, and still intended** | #74 | 1 |

Eight closes and one corrective comment (#87) — nine notifications, best sent in one sitting.

---

## Close — implemented, verified in code

### #43 — save single-subject CIs alongside group CIs

**Done.** [`generateCI.R:83`](../R/generateCI.R#L83) takes `save_individual_cis=FALSE`, and
[`generateCI.R:259-269`](../R/generateCI.R#L259-L269) writes each participant's CI to
`<targetpath>/individual_cis`, scaled by its own `individual_scaling` /
`individual_scaling_constant` arguments — exactly the "one single command for single subjects
as well as the group level" the issue asked for.

> Draft comment:
> Done. `generateCI(..., participants = <col>, save_individual_cis = TRUE)` writes each
> participant's CI to an `individual_cis` subfolder of `targetpath` in the same call that
> produces the group CI, and scales them independently via `individual_scaling` and
> `individual_scaling_constant`. As this issue anticipated, they were already being computed —
> the group CI is their average — so this only added the saving.

### #27 — group CIs in which every participant has equal weight

**Done, and both of the follow-up checkboxes with it.** When `participants` is supplied,
[`generateCI.R:229-282`](../R/generateCI.R#L229-L282) computes one CI per participant and then
takes the *unweighted* mean across participants
([`generateCI.R:282`](../R/generateCI.R#L282)) — so trial count no longer confers weight. It
arrived as the `participants` argument rather than as a separate flag.

The two checkboxes in the thread are also satisfied:
- `stimuli_params$gender_neutral` is gone; the code uses the generic `params`
  ([`generateCI.R:189`](../R/generateCI.R#L189), [`:321`](../R/generateCI.R#L321)).
- The `t.test` z-map path uses `pid.cis` — the per-participant CIs — rather than per-trial
  noise images when `participants` is given ([`generateCI.R:349`](../R/generateCI.R#L349)).

**One checkbox is genuinely not done**, and I suggest saying so rather than closing over it:
`generateCI2IFC()` and `batchGenerateCI2IFC()` still do not accept `participants`
([`generateCI2IFC.R:61`](../R/generateCI2IFC.R#L61),
[`batchGenerateCI2IFC.R:51`](../R/batchGenerateCI2IFC.R#L51)). Both are the legacy 2IFC-named
wrappers, so the honest resolution is "use `generateCI()`", not a backport.

> Draft comment:
> Implemented, as the `participants` argument rather than a separate flag: with
> `participants` supplied, `generateCI()` computes one CI per participant and averages those
> unweighted, so a participant with fewer trials no longer carries less weight.
>
> Both follow-ups in this thread came with it — the hard-coded `gender_neutral` label is gone
> in favour of the generic `params`, and the `t.test` z-map path operates on the
> per-participant CIs (`pid.cis`) rather than on per-trial noise images when `participants` is
> given.
>
> Not done, deliberately: `generateCI2IFC()` and `batchGenerateCI2IFC()` still don't take
> `participants`. They are the legacy wrappers; `generateCI()` is the supported entry point
> and has taken `participants` for some time. Closing on that basis.

### #13 — "Update documentation" (2016, empty body)

**Superseded.** Every one of the 17 exported functions has a roxygen block and a generated
`.Rd`; two vignettes ship (`getting-started.Rmd`, `reverse-correlation-walkthrough.Rmd`); the
README carries the architecture walkthrough and the `.Rdata` anatomy; and the docs have since
been through a CRAN review. A 2016 issue with no body cannot be checked off against any of
that.

> Draft comment:
> Closing as superseded. The package now has roxygen documentation for every exported
> function, two vignettes (`getting-started`, `reverse-correlation-walkthrough`), a README
> covering the architecture and the `.Rdata` format, and it has been through CRAN review since.
> This issue has no body, so there is nothing left in it that is specific enough to act on —
> anything still missing is better filed fresh.

---

## Close — moot, the premise no longer holds

### #92 — add reference datasets

Both reasons this was *reopened* in 2017 are gone: there is **no `inst/extdata` directory at
all**, so the 60.5 Mb `installed size` NOTE cannot occur, and the undistributable `MNES.jpg`
is not in the package.

The original motivation — "standardize testing of new features" — is now met by other means:
`tests/testthat/helper-fixtures.R` builds synthetic base images on the fly (never a real
photo, so no licensing exposure), `test-regression-baseline.R` pins the numeric output of the
default pipeline, and `tools/compare-release-output.R` diffs the working tree against a
released version. That is a stronger guarantee than a shipped dataset would give, and it costs
nothing in tarball size.

What is *not* covered is a dataset for **user** convenience — a worked example someone can run
without generating stimuli first. That is a real want, but it is a different issue from this
one, and it collides with the same size and licensing constraints.

> Draft comment:
> Closing: both reasons this was reopened are resolved. There is no `inst/extdata` in the
> package any more, so the 60.5 Mb installed-size NOTE is gone, and the undistributable
> `MNES.jpg` is no longer shipped.
>
> The testing motivation is met a different way — the test suite generates synthetic base
> images on the fly, a golden-master test pins the numeric output of the default pipeline, and
> a release gate diffs every output against a released version before each release. That
> catches more than a bundled reference dataset would, without adding to the tarball.
>
> A small dataset for *users* to run examples against is still a reasonable idea, but it is a
> different request with the same size and licensing constraints attached. Worth a fresh issue
> if it is wanted.

### #52 — should we use random seeds? (question)

**Answered by the design that shipped.** `generateStimuli2IFC()` takes `seed = 1`
([`generateStimuli2IFC.R:46`](../R/generateStimuli2IFC.R#L46)) — a documented, fixed default
that the user can vary. That takes the "con" (no recovery for a user who lost the `.Rdata`,
reduced comparability) as the default, while leaving the "pro" (avoiding stimulus-set
monoculture) one argument away. The `.Rdata` file also records the seed, so a set is
reproducible even when a non-default one was used.

Note there is a **separate, live** defect nearby — `set.seed(seed)` is never undone, so the
user's RNG stream is left where stimulus generation put it. That is now #189; it is not what
this one asks.

> Draft comment:
> Closing as answered by what shipped. `generateStimuli2IFC()` has a documented `seed`
> argument defaulting to `1`: reproducible by default — which is the side this issue's "cons"
> argue for — and one argument away from a random set for anyone who wants to avoid stimulus-set
> monoculture. The seed is also stored in the `.Rdata` file, so a non-default set stays
> recoverable.
>
> Unrelated but nearby, and tracked separately: the seed is set and never restored, so the
> user's RNG stream is left where stimulus generation put it.

---

## Keep open, rescope

### #87 — refactor batch CI generation (deprecates `batchGenerateCI`)

Three of the four things in this thread shipped; the checkboxes are stale in both directions.

| item | state |
|---|---|
| `autoscaling` parameter on `generateCI` for use with `participants` | **Done** — and named `individual_scaling`, the rename asked for in the 2017-08-10 comment ([`generateCI.R:85`](../R/generateCI.R#L85)) |
| the target line of code `generateCI(..., participants=, individual_scaling=, save_individual_cis=)` | **Valid today**, argument for argument |
| single function vs `generateCI` + `scaleCI` | **Settled** — 2017-09-25, stick with `generateCI()`, matching the European Review specification |
| rewrite `batchGenerateCI` as a wrapper around `generateCI` | ticked in 2017, but it is **not** a wrapper today — [`batchGenerateCI.R:51`](../R/batchGenerateCI.R#L51) has its own loop |
| "function deprecated" warning on `batchGenerateCI` | ticked in 2017, but **there is no such warning in the code today** |
| `condition` argument, participants nested in condition | **Not done** |

So the residue is: no deprecation warning, no `condition` argument, and `batchGenerateCI` is
still a parallel implementation rather than a wrapper. Worth keeping, with the checkboxes
corrected — the two ticked-but-absent items are the sort of thing that makes a tracker
untrustworthy.

> Draft comment:
> Updating rather than closing — most of this shipped, and two checkboxes are ticked for work
> that is not in the code.
>
> Done: `generateCI()` takes `participants`, `save_individual_cis`, `individual_scaling` and
> `individual_scaling_constant`, so the target call in the first comment is valid today, with
> `individual_scaling` being the rename of `autoscale` asked for there. The
> one-function-vs-several question was settled in favour of `generateCI()` doing both.
>
> Still open, contrary to the ticked boxes: `batchGenerateCI()` is not a wrapper around
> `generateCI()` — it still has its own loop — and it emits no deprecation warning. The
> `condition` argument (participants nested in condition, group CIs per condition) was never
> implemented.

---

## Keep open — real, specified, unimplemented

Verified as genuinely absent from the code; each has enough of a specification to act on.

| issue | what | note |
|---|---|---|
| #71 | z-map comparing **two** CIs | Today's z-maps test one CI against zero — `quick` scales the blurred CI, `t.test` tests per-pixel across trials or participants ([`generateCI.R:303-356`](../R/generateCI.R#L303-L356)). A two-CI contrast is new work. |
| #15 | smoothing option on `generateCI` | `sigma` exists but smooths only the **z-map** ([`generateCI.R:307`](../R/generateCI.R#L307)); the CI itself is never blurred. The issue's exact proposal — `generateCI(smoothing = ...)` — is unimplemented. Overlaps #37 and #38. |
| #85 | infoVals for many CIs at once | `computeInfoVal2IFC()` takes a single `target_ci` ([`computeInfoVal2IFC.R:68`](../R/computeInfoVal2IFC.R#L68)). Cheap to add and it would amortise the reference distribution across CIs, which is the expensive part. |
| #76 | oblong base images | Still square-only. The maintainer's 2017 caveat stands: white noise is easy, sinusoid noise needs rework to avoid stretching. |
| #47 | ROI analysis + small-volume correction | The `mask` argument gets partway (restrict the CI to a region) but does no ROI statistics or correction. |
| #46 / #38 | cluster-based permutation test / Random Field Theory | The two competing routes to the same end: valid statistics on smooth CIs. Neither implemented. #38 carries the better reference material; #46 is likelier to be tractable. Worth deciding between them rather than carrying both. |
| #37 | blurring optimized to the base image's spatial frequency composition | Unimplemented. Related to #15. |
| #54 | 2IFC with randomly paired stimuli rather than ori+inv | The only one of these with a **methodological** claim behind it (the simulations in the thread: noise/anti-noise recovers the sign of a pixel deviation but not its magnitude). If that holds it also affects significance testing. Highest scientific value in this group. |
| #35 | CI timelapses | Unimplemented. The thread already contains the answer to its own difficulty: dump frames, make video assembly optional. |
| #22 | CI database + community building | Large, unimplemented, and a scope question rather than a coding one. |
| #69 | small tutorials | Two vignettes now ship, but neither covers the noise-only recipe in this issue's body. Narrow enough to keep as "document generating noise-only images". |

## The thin issues — title only, empty body, no discussion

Five issues name a plausible feature with nothing that can be verified as done or not-done.
Split on whether the title alone is a specification:

**Kept.** #6 (direct support for 4AFC and other variants) — the package is 2IFC-only end to
end, so the scope is unambiguous even with no body. #7 ((masked) correlations between lists of
CIs) — the title is nearly a function signature, and `mask` already exists to build on. #10
(meta-reverse correlation) — its body does specify the function: list of CIs + weights → a
weighted-average CI, plus an optional anti-image.

**Closed.** #68 below. #9 was also closed, in error — see below.

### #68 — adaptive reverse correlation

A research direction rather than a feature: adaptive RC is a decision about the whole task
loop — how stimuli are selected from trial to trial — not something checkable against the
code. Nine years with no body and no discussion.

> Draft comment:
> Closing as unspecified. This is a title with no body and no discussion, and adaptive reverse
> correlation is a research direction rather than a feature — it is a decision about the whole
> task loop (how stimuli are selected from trial to trial), not something that can be checked
> off against the code. Whatever eventually gets built will not come from this issue.
>
> Reopening or filing fresh with a concrete design in mind is welcome.

### #9 — save analysis metadata to a log file? ⚠️ **Closed on a false premise — since reopened**

**The verdict below was wrong, the comment stating it was posted, and the issue has since
been reopened with a correction.** It claimed
`generateCI()` returns the scaling method and constant. It does not:
[`generateCI.R:369-371`](../R/generateCI.R#L369-L371) returns `ci`, `scaled`, `base`,
`combined` and optionally `zmap`, and `applyScaling()`
([`:463-499`](../R/generateCI.R#L463-L499)) discards both the method and the computed constant
— `autoscale()` *prints* its constant to the console rather than returning it. The `.Rdata`
carries stimulus-*generation* parameters, not analysis settings, so it does not hold them
either.

So the metadata #9 asks for is genuinely not recoverable from the results, and the issue
describes real missing functionality. Caught by the Codex review on PR #172, after posting.

The lesson is the one this file was supposed to embody: every *other* verdict here cites a
file and line, and this one asserted a return value from memory. An unverified claim reads
exactly like a verified one once it is written down.

> Correction posted on reopening:
> Reopening — I closed this on a claim that is simply wrong. `generateCI()` returns `ci`,
> `scaled`, `base`, `combined` and optionally `zmap`; neither the scaling method nor the
> computed scaling constant is among them, and `applyScaling()` discards both.
> `autoscale()` prints its constant to the console rather than returning it, and the `.Rdata`
> file holds stimulus-generation parameters, not analysis settings.
>
> So the metadata this issue asks for is not recoverable from the results today, and the
> request stands as originally filed.

## Aspirational since 2017

Neither is stale in the sense of being done or wrong; stale in the sense that nothing moved in
nine years. Both were the maintainer's call, not a code check.

- **#74 — GUI. Kept.** Wanted in 2017, with a jamovi module floated. Still intended.
- **#53 — MDS plots with CIs on tooltip in Shiny. Closed.** A neat exploratory tool with a
  screenshot attached, but no MDS functionality exists in the package to attach it to.

> Draft comment for #53:
> Closing. The package has no MDS functionality for this to attach to, and adding one is a
> larger question than the viewer itself — rcicr computes classification images and leaves the
> downstream analysis to the user's own tooling, where `ggplot2`/`plotly` already do image
> tooltips well. The screenshot is still a good demonstration of the idea for anyone building
> it outside the package.

---

## Sequencing note for item 44

All nine postings went out together on 2026-08-08 — one sitting, as intended, and **before**
the 21 new issues opened from the backlog as #174–#194. That prerequisite was met.

The tracker afterwards: **17 open** — #6, #7, #10, #15, #22, #35, #37, #38, #46, #47, #54,
#69, #71, #74, #76, #85, #87 — and 18 after #9 was reopened.
