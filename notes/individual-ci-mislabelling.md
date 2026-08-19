# Which rcicr you had, and whether it carried the individual-CI bug

Reference for anyone checking whether a stored analysis is affected by the mislabelling fixed
in 1.3.0 (issue #261). The short answer is that **the version number alone does not tell you** —
for anything installed from GitHub before 2026, the date and the *branch* matter more than the
version.

## What the bug was

`generateCI()` builds one classification image per participant when given a `participants`
vector, and writes them to `<targetpath>/individual_cis/` when `save_individual_cis = TRUE`.

Inside that loop, two different orderings of the same participants were used:

```r
pids <- as.numeric(factor(participants))          # selects trials: SORTED order
...
saveToImage(..., unique(participants)[obs], ...)  # names the file: APPEARANCE order
```

`factor()` numbers its levels by `sort(unique(x))`. `unique()` preserves order of first
appearance and never sorts. So `obs = 1` selects the trials of the *alphabetically first*
participant but writes them under the name of the participant who appears *first in the data*.
When those two are not the same person, the file is misnamed.

The images themselves were always computed correctly. Only the filenames were wrong, and the
group-level classification image — a mean across participants — is unaffected, because a mean
does not depend on order.

**Sorting is lexical, so ordinary data triggers it.** `"p10"` sorts before `"p2"`. A study whose
participants appear in the natural order `p1, p2, ..., p12` is affected from the tenth
participant onward. Measured against the released 1.2.3 with exactly that design, **11 of 12
files carried the wrong participant's image**.

| `participants`, in collection order | affected |
|---|---|
| `"p1"` … `"p9"` | no |
| `"p1"` … `"p10"` or beyond | **yes** |
| `"p01"` … `"p12"` (zero-padded) | no |
| `"1"` … `"12"` (character) | **yes** |
| `1` … `12` (numeric) | no |

**Only one call is affected**: `generateCI(participants = ..., save_individual_cis = TRUE)`.
`batchGenerateCI()`, `batchGenerateCI2IFC()` and `generateCI2IFC()` cannot reach it — they call
`generateCI()` once per group with `participants = NA`, so the loop never runs, and they name
each image from the group value they are already iterating over. Verified against the released
1.2.3, where the direct call gets 11 of 12 wrong: both batch functions got **0 of 12** wrong and
never created `individual_cis/` at all.

Introduced 2017-08-15 (`8a74974`, merged in PR #88). Fixed 2026-08-17 (`706882e`), released in
1.3.0.

## Every version, and where it lived

| version | date | channel | carries the bug |
|---|---|---|---|
| 0.2.1 | 2014-07-31 | CRAN | no — predates the feature |
| 0.2.4 | 2014-09-30 | CRAN | no |
| 0.2.5 | 2014-11-17 | CRAN | no |
| 0.2.6 | 2014-12-19 | CRAN | no |
| 0.3.0 | 2015-01-23 | CRAN | no |
| 0.3.2 | 2015-07-20 | CRAN | no |
| 0.3.2.1 | 2016-02-05 | CRAN | no |
| **0.3.4.1** | **2016-07-13** | **CRAN — last release there, archived 2021-06-08** | **no** |
| 0.3.4.2 | 2016-09-22 | GitHub `master` only, under `pkg/` — never on CRAN | no |
| 0.4.0 | `development` branch only, from 2016-10-26 | GitHub, `ref = 'development'` | **only from 2017-08-15** (see below) |
| 0.4.1 | `development` from 2021-09-23; default branch from 2021-12-28 | GitHub | yes |
| 1.0.0 | default branch, from 2022-09-02 | GitHub | yes |
| 1.0.1 | commit 2023-01-13, tagged 2026-07-28 | GitHub tag | yes |
| 1.1.0 | 2026-07-27 | GitHub release | yes |
| 1.2.0 | 2026-07-28 | GitHub release | yes |
| 1.2.1 | 2026-07-28 | GitHub release | yes |
| 1.2.2 | 2026-08-07 | GitHub release | yes |
| 1.2.3 | 2026-08-08 | GitHub release | yes |
| **1.3.0** | **2026-08-18** | **GitHub release + CRAN submission** | **fixed** |

**`0.4.0` is the trap.** That version string sat on the `development` branch from 2016-10-26 to
2021-09-23 — five years — and the bug entered partway through, on 2017-08-15. Someone who
reports "I used rcicr 0.4.0" cannot be classified from the version alone; the install date
separates clean from affected, and the branch decides whether they could have had it at all.

## What "the latest version" gave you, by date

### From CRAN — `install.packages('rcicr')`

| you installed | you got | affected |
|---|---|---|
| 2014-07-31 → 2014-09-30 | 0.2.1 | no |
| 2014-09-30 → 2014-11-17 | 0.2.4 | no |
| 2014-11-17 → 2014-12-19 | 0.2.5 | no |
| 2014-12-19 → 2015-01-23 | 0.2.6 | no |
| 2015-01-23 → 2015-07-20 | 0.3.0 | no |
| 2015-07-20 → 2016-02-05 | 0.3.2 | no |
| 2016-02-05 → 2016-07-13 | 0.3.2.1 | no |
| 2016-07-13 → 2021-06-08 | 0.3.4.1 | no |
| 2021-06-08 → now | **nothing** — archived, so `install.packages('rcicr')` fails. The tarball stays downloadable from [CRAN's Archive](https://cran.r-project.org/src/contrib/Archive/rcicr/), and 0.3.4.1 from there is clean | no |

**No CRAN user has ever received a mislabelled file**, at any point in the package's seven years
there. The last CRAN release predates the `save_individual_cis` option by thirteen months.

### From the GitHub default branch — `install_github('rdotsch/rcicr')`

**Until 2021-12-28 this did not work at all**, and that is the single most important row in
this note. The default branch was `master`, which still carried the R-Forge layout: the
package sat under `pkg/` with no `DESCRIPTION` at the repository root, so `install_github()`
had nothing to install. `master` also stopped receiving code after 2016-10-21 and was frozen
at a README edit on 2017-01-31 — **it never contained the bug**.

All development happened on the `development` branch, which is where the bug lived from
2017-08-15. On 2021-12-28 that branch became the default (`69e3933`, "Removed old pkg master,
development is the new master"), and only from then does a plain `install_github()` return
affected code.

| you installed | result | affected |
|---|---|---|
| 2016-06-23 → 2016-09-18 | **install fails** — repository created, no package code yet | — |
| 2016-09-18 → 2021-12-28 | **install fails** — default branch is `master`; package under `pkg/`, no root `DESCRIPTION`. `install_github('rdotsch/rcicr', subdir = 'pkg')` did work, giving 0.3.4.2 | no |
| 2021-12-28 → 2022-09-02 | 0.4.1 | **yes** |
| 2022-09-02 → 2023-01-13 | 1.0.0 | **yes** |
| 2023-01-13 → 2026-07-26 | 1.0.1 | **yes** |
| 2026-07-26 → 2026-08-17 | 1.x.y or 1.x.y.9000 | **yes** |
| 2026-08-17 → now | 1.2.3.9000, then 1.3.0 | **fixed** |

### From the development branch — `install_github('rdotsch/rcicr', ref = 'development')`

Before December 2021 this was the only way to get a working install from GitHub, and it is
what the README of the time told people to use — with a warning attached:

> **Using the dev version (AT YOUR OWN RISK!)** … We do not recommend that you use the
> development version for analyses meant for publication. Please wait until we make a new
> version available on CRAN for that.

So anyone who took the bug from GitHub before December 2021 had to ask for this branch
deliberately, having been told not to publish from it. That warning stood from 2017-01-09 —
seven months before the bug landed — until it was removed on 2022-09-02.

Half of it expired in June 2021: CRAN archived the package, so "wait until we make a new version
available on CRAN" pointed at a release that was never coming. Unaffected code remained
obtainable, though — the archived 0.3.4.1 tarball is still downloadable from
[CRAN's Archive](https://cran.r-project.org/src/contrib/Archive/rcicr/), and
`install_github('rdotsch/rcicr', subdir = 'pkg')` returned the clean 0.3.4.2 on `master`. What
was unavailable after June 2021 was a *maintained* clean version, not a clean version.

| you installed | you got | affected |
|---|---|---|
| before 2016-10-26 | **install fails** — no root `DESCRIPTION` yet | — |
| 2016-10-26 → 2017-08-15 | 0.4.0 | no |
| **2017-08-15 → 2021-09-23** | **0.4.0** | **yes** |
| 2021-09-23 → 2021-12-28 | 0.4.1 | **yes** |
| after 2021-12-28 | the branch became the default; see the table above | **yes** |

### From the latest GitHub release — `install_github('rdotsch/rcicr@*release')`

| you installed | you got |
|---|---|
| before 2026-07-27 | **nothing** — no releases existed to resolve against |
| 2026-07-27 → 2026-07-28 | 1.1.0 |
| 2026-07-28 → 2026-08-07 | 1.2.0, then 1.2.1 the same day |
| 2026-08-07 → 2026-08-08 | 1.2.2 |
| 2026-08-08 → 2026-08-18 | 1.2.3 |
| 2026-08-18 → now | **1.3.0** |

Every GitHub release before 1.3.0 carries the bug. `v1.0.1` has a tag but no GitHub *release*,
and it was tagged retroactively in 2026 against a 2023 commit, so `@*release` never resolved to
it.

## How to check your own analysis

The affected call needs **both** things: a `participants` vector *and* `save_individual_cis =
TRUE`. So the question to answer is:

> **Did a single `generateCI()` call produce one classification image per participant?**

If the answer is no, nothing here touches you, and there are several ways to reach it that do
not depend on any file surviving.

**You did not make per-participant images at all.** A group-level classification image is
unaffected — it is a mean across participants and does not depend on their order.

**You made them with the batch functions.** `batchGenerateCI()` and `batchGenerateCI2IFC()`
cannot reach the defect: they call `generateCI()` once per group with `participants = NA`, so
the loop never runs, and they name each image from the group value they are iterating over.
This is also the route the documentation has recommended since the CRAN era. `generateCI2IFC()`
does not expose either argument.

Checked across the whole affected range, not only the latest release: `batchGenerateCI()`
passes `participants = NA` explicitly in every version from 2017-08-15 to 1.2.3, and neither
batch function mentions `save_individual_cis` in any of them. Measured too — run against the
released 1.2.3 with participants `p1 … p12` in collection order, where the direct call gets 11
of 12 filenames wrong, both batch functions got 0 of 12 wrong and created no `individual_cis`
directory at all.

**You still have the output.** The images live in a directory called `individual_cis`, and
nothing else in the package writes there — the only such write is in `computeParticipantCIs()`,
reachable only from `generateCI()`. That has held in every version: exactly one file contains
the write, `R/generateCI.R` from 2017 through 1.2.3 and `R/ci-compute.R` after the loop was
extracted. Individual files are named `ci_<participant>.png`, or
`antici_<participant>.png` for an antiCI, so the naming outlives a renamed folder or moved
contents — as a clue, not a verdict. `filename` is a free-form argument, so a group CI written
by `generateCI(save_as_png = TRUE, filename = "p3")` lands at `<targetpath>/ci_p3.png` and is
indistinguishable by name alone (`R/generateCI.R:207`, `saveToImage()` at 392-406). Once the
directory structure is gone, confirm a loose file from the script or the data before treating
it as an individual CI.

**Telling one folder of PNGs from another.** Every other writer in the package puts its images
*flat* into the path you gave it — `individual_cis/` is the only subdirectory any of them
creates, which is what makes its presence such a clean signal. The filenames differ too:

| what wrote it | where | filename |
|---|---|---|
| `generateCI(save_individual_cis = TRUE)` | `<targetpath>/individual_cis/` | `ci_<participant>.png` — bare participant id |
| `generateCI(save_as_png = TRUE)`, `generateCI2IFC()` | `<targetpath>/` | `ci_<filename>.png`, defaulting to the base image name |
| `batchGenerateCI()`, `batchGenerateCI2IFC()` | `<targetpath>/` | `ci_<baseimage>_<by>_<unit>.png` |
| `autoscale(save_as_pngs = TRUE)` | `<targetpath>/` | `<name>_autoscaled.png` |
| `plotZmap()`, `generateCI(zmap = TRUE)` | `<zmaptargetpath>/` | `<filename>.png` |
| `generateStimuli2IFC()` | `<stimulus_path>/` | `<label>_<base>_<seed>_<NNNNN>_ori.png` and `_inv.png` |

An `antiCI` swaps the `ci_` prefix for `antici_` throughout. So a bare `ci_p3.png` came from the
affected call, while `ci_face_participant_p3.png` came from a batch function and is fine — the
batch name always carries the base image and the grouping column, because it is built from them.

**You still have the script but not the output.** Read every `generateCI()` call that passes
participant identifiers. Searching for the string `save_individual_cis` is *not* sufficient on
its own: it is the sixth formal, so `generateCI(stim, resp, "face", rdata, pids, TRUE, ...)`
sets it positionally, and `do.call()` with an assembled list hides it too. Search for
`individual`, then read by hand any call passing six or more arguments positionally.

**You still have the data but neither.** Recompute with 1.3.0 and compare. This is the only
fully conclusive answer, and it hands you corrected output as a side effect: the responses and
the stimulus `.Rdata` are all that is needed, and a classification image is deterministic given
them.

**You have only the published figure.** Then it cannot be verified directly, and the ordering
question below is the best available evidence.

### If it was the affected call: was the ordering safe?

```r
identical(unique(participants), sort(unique(participants)))   # TRUE = the names were correct
```

"Sorted" means *lexically* sorted for text labels, per the table earlier: sorting a data frame
by a `"p1"`-style column does not make it safe once you reach ten participants.

Recovery is a rename, not a recomputation: the file named `unique(participants)[i]` holds the
image of `sort(unique(participants))[i]`.

## What it did to a second-stage analysis

The advisory tells readers that the mislabelling mostly costs effects rather than inventing
them. That is a claim about a *permutation*, and the two summaries it invites — "it only adds
noise" and "it could have produced spurious findings" — are each half right in ways that matter,
so it is measured rather than argued. What comes out is that the aggregate false positive rate
is untouched while individual verdicts flip both ways, and that in some designs the shuffle
raises apparent support rather than only destroying it. `tools/simulate-mislabelling-impact.R` runs 20,000
iterations per cell at α = 0.05, seeded; the numbers below are its output.

The design simulated is the one that matters: raters judge each individual CI, and those
judgments are related back to something about the participant the filename names. The bug
permutes that join.

**Almost nothing survives the permutation.** With IDs `p1 … pN` in collection order, the
number of files still carrying their own participant is 1 at N = 12, 20, 30 and 50 (and all
N at N ≤ 9, where lexical order still matches collection order). Zero-padded `p01 … p12`
is unaffected, as the version table above says.

**With no true association, the mislabelling does not raise the false positive rate.**
Rejection stays at α in
every design and at every N — 0.046 to 0.053 against a nominal 0.05. Those individual rates
should not be read too closely: at nsim = 20,000 the Monte Carlo SE at 0.05 is 0.0015, and each
test brings its own size (the Welch *t* in scenario B is mildly conservative at N = 12, 0.0471
in the correctly labelled column, where scenarios A and C use `cor.test()` and sit nearer
0.050). The statistic that carries the claim is therefore the **paired** difference, both
analyses seeing the same outcome vector within an iteration; the script prints it for all
twelve null cells across the three designs, and it runs from −0.005 to +0.003, each within 2.1
paired SEs of zero. The result is structural rather than a lucky set of draws: permuting one
side of a pair that is independent of the other leaves it independent, so the null distribution
is unchanged by construction. **The bug does not raise the false positive rate.**

It does change individual verdicts, which is a different question and the one a researcher is
actually asking. The `decisions_changed` column counts datasets where the two analyses disagree
about significance: 0.075 to 0.099 across the null cells, roughly half of them a rejection the
correctly labelled data would not have produced. A result that became significant only on
mislabelled files is a false positive handed over by the bug, even though the rate across many
studies is untouched.

**With a true association and an ordinary design — a covariate unrelated to the order
participants happened to be labelled in — the association is destroyed**, not merely weakened.
Detection collapses to α; the runs that then fail to reject — 95% of them at ρ = 0.5,
N = 50 — are the false negatives:

| N | ρ | detected, correct labels | detected, mislabelled | mean *r* recovered |
|---|---|---|---|---|
| 12 | 0.3 | 0.160 | 0.047 | 0.001 |
| 50 | 0.3 | 0.572 | 0.049 | −0.001 |
| 12 | 0.5 | 0.399 | 0.047 | −0.000 |
| 50 | 0.5 | **0.970** | **0.047** | −0.001 |

Significant results in the affected column split evenly by sign (≈ 0.024 each way), which is
what the permutation null looks like.

**How much damage depends on how far out of order the identifiers were, and A and C fix that at
close to the worst case.** `p1 … pN` in collection order leaves exactly one file correctly
paired past the tenth participant, so "destroyed" describes that scheme rather than every
affected study. Scenario E sweeps the severity using the real mechanism — build an appearance
order, let `mislabel_perm()` derive the permutation from it — at N = 50 with ρ = 0.5, which
correct labels detect 97% of the time:

| identifiers | share correctly paired | detected once mislabelled |
|---|---|---|
| one pair out of sequence | 0.96 | 0.949 |
| two pairs | 0.92 | 0.924 |
| five pairs | 0.80 | 0.818 |
| ten pairs | 0.60 | 0.560 |
| every pair transposed | 0.00 | 0.081 |
| `p1 … p50` in collection order | 0.02 | 0.047 |

Detection tracks the share still correctly paired almost one-for-one until the pairing is nearly
gone. Zero-padded identifiers entered with a couple of participants out of sequence are
therefore affected but barely dented; the `p1 … pN` scheme is the case where the association
disappears. That is the file-drawer case, and it is the common labelling scheme — but it is the
severe end of the range, not the whole of it.

**The exception is a covariate that tracks labelling order**, where the permutation is no
longer exchangeable with the design and the residual can come out either sign. With condition
assigned in blocks by collection order and d = 0.8:

| N | share of participants moved across conditions | detected, correct | detected, mislabelled | of those, in the predicted direction | opposite |
|---|---|---|---|---|---|
| 12 | 0.50 | 0.222 | 0.032 | 0.016 | 0.017 |
| 20 | 0.80 | 0.387 | 0.143 | 0.001 | **0.142** |
| 30 | 0.47 | 0.559 | 0.037 | 0.025 | 0.012 |
| 50 | 0.24 | 0.794 | 0.262 | 0.262 | 0.000 |

At N = 20, where four fifths of participants are swapped across the condition boundary, the
effect comes back **reversed**: 14.2% of runs are significant in the wrong direction against
0.1% in the right one. At N = 50, where only a quarter cross, it survives attenuated and
correctly signed. A covariate that *is* collection order (testing date, cohort, session
number) behaves the same way — at N = 50, ρ = 0.5 a mean *r* of 0.22 is still recovered, while
at N = 20 the mean *r* is −0.11.

**Can it ever hand back a result that supports the hypothesis?** Scenario B at d = 0 settles
the *rate* for an outcome which is exchangeable noise — permuting it leaves the null
distribution alone, so nothing is inflated on average, though as the `decisions_changed` figures
above show, individual datasets still flip in both directions. The harder case is an outcome
carrying real structure tied to collection order —
drift, practice, season — with no true condition effect, where the permutation has something to
move around. Scenario D runs it: condition blocked by collection order, `y` a real trend over
that order, so the trend is confounded with condition before the bug is involved at all.

**The assignment scheme decides the answer, so both are run.** Contiguous condition blocks put
every early observation in one condition and every late one in the other, which already
maximises the trend contrast under *correct* labels — at N = 12 the expected contrast in
position units is 6.00, and the permutation takes it to 0.00. Any permutation can only reduce
that, so testing blocks alone would make "the shuffle never helps" true by construction rather
than measured. Alternating assignment is the opposite case: the correct pairing spreads the
trend evenly across conditions, leaving 1.00 of contrast, and the permutation raises it to 2.33.

| N | trend | assignment | predicted direction, correct labels | mislabelled | |
|---|---|---|---|---|---|
| 20 | 0.3 | blocked | 0.177 | 0.006 | destroyed |
| 50 | 1.0 | blocked | 1.000 | 0.793 | reduced |
| 20 | 1.0 | blocked | 0.907 | 0.000 | destroyed, 0.128 lands opposite |
| 12 | 0.3 | alternating | 0.026 | 0.039 | **raised** |
| 12 | 0.6 | alternating | 0.020 | 0.047 | **raised** |
| 12 | 1.0 | alternating | 0.010 | 0.042 | **raised** |
| 50 | 0.3 | alternating | 0.025 | 0.028 | **raised** |

**So the mislabelling can hand back apparent support**, in 9 of the 12 alternating cells, and by
a factor of four at N = 12, trend = 1.0. An earlier draft of this note said it never could; that
was an artifact of testing only the blocked scheme, and it was wrong.

The sign the permutation pushes is fixed by the ID scheme and the assignment scheme together,
and nothing in lexical sorting refers to a hypothesis. How often that sign agrees with what a
study predicted is **not** answered here: every cell fixes the trend and the "predicted"
direction as positive rather than sampling them, so these runs cannot support a claim about
directional inflation across a literature. Naming conventions, condition coding and the usual
direction of practice or seasonal trends may well be related to one another in real research;
showing that they are not would take a different study from this one.

What this does establish is the consequence for a single analysis: it can be significant in the
predicted direction and still be no evidence for anything, having been computed on the wrong
pairing. That is why the advisory asks for a re-run whichever way an affected analysis came out,
rather than offering an average as an all-clear.

So the summary the advisory gives is the right one, with the scope it carries: on average an
effect is lost rather than invented, an unpublished null is the likeliest casualty, the false
positive *rate* is untouched — and no single affected analysis can be cleared on any of that,
because the verdict on an individual dataset can flip either way, and a design whose IDs encode
something can distort, reverse or amplify an effect instead of erasing it.
