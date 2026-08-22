# Which rcicr you had, and whether it carried the individual-CI bug

Reference for anyone checking whether a stored analysis is affected by the mislabelling fixed in 1.3.0 (issue #261). The short answer is that **the version number alone does not tell you**: for anything installed from GitHub before 2026, the date and the *branch* matter more than the version.

## What the bug was

`generateCI()` builds one classification image per participant when given a `participants` vector, and writes them to `<targetpath>/individual_cis/` when `save_individual_cis = TRUE`.

Inside that loop, two different orderings of the same participants were used:

```r
pids <- as.numeric(factor(participants))          # selects trials: SORTED order
...
saveToImage(..., unique(participants)[obs], ...)  # names the file: APPEARANCE order
```

`factor()` numbers its levels by `sort(unique(x))`. `unique()` preserves order of first appearance and never sorts. So `obs = 1` selects the trials of the *alphabetically first* participant but writes them under the name of the participant who appears *first in the data*. When those two are not the same person, the file is misnamed.

The images themselves were always computed correctly. Only the filenames were wrong, and the group-level classification image (a mean across participants) is unaffected, because a mean does not depend on order.

**Sorting is lexical, so ordinary data triggers it.** `"p10"` sorts before `"p2"`. A study whose participants appear in the natural order `p1, p2, ..., p12` is affected from the tenth participant onward. Measured against the released 1.2.3 with exactly that design, **11 of 12 files carried the wrong participant's image**.

| `participants`, in collection order | affected |
|---|---|
| `"p1"` … `"p9"` | no |
| `"p1"` … `"p10"` or beyond | **yes** |
| `"p01"` … `"p12"` (zero-padded) | no |
| `"1"` … `"12"` (character) | **yes** |
| `1` … `12` (numeric) | no |

**Only one call is affected**: `generateCI(participants = ..., save_individual_cis = TRUE)`. `batchGenerateCI()`, `batchGenerateCI2IFC()` and `generateCI2IFC()` cannot reach it: they call `generateCI()` once per group with `participants = NA`, so the loop never runs, and they name each image from the group value they are already iterating over. Verified against the released 1.2.3, where the direct call gets 11 of 12 wrong: both batch functions got **0 of 12** wrong and never created `individual_cis/` at all.

Introduced 2017-08-15 (`8a74974`, merged in PR #88). Fixed 2026-08-17 (`706882e`), released in 1.3.0.

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

**`0.4.0` is the trap.** That version string sat on the `development` branch from 2016-10-26 to 2021-09-23 (five years) and the bug entered partway through, on 2017-08-15. Someone who reports "I used rcicr 0.4.0" cannot be classified from the version alone; the install date separates clean from affected, and the branch decides whether they could have had it at all.

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

**No CRAN user has ever received a mislabelled file**, at any point in the package's seven years there. The last CRAN release predates the `save_individual_cis` option by thirteen months.

### From the GitHub default branch — `install_github('rdotsch/rcicr')`

**Until 2021-12-28 this did not work at all**, and that is the single most important row in this note. The default branch was `master`, which still carried the R-Forge layout: the package sat under `pkg/` with no `DESCRIPTION` at the repository root, so `install_github()` had nothing to install. `master` also stopped receiving code after 2016-10-21 and was frozen at a README edit on 2017-01-31. **It never contained the bug.**

All development happened on the `development` branch, which is where the bug lived from 2017-08-15. On 2021-12-28 that branch became the default (`69e3933`, "Removed old pkg master, development is the new master"), and only from then does a plain `install_github()` return affected code.

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

Before December 2021 this was the only way to get a working install from GitHub, and it is what the README of the time told people to use, with a warning attached:

> **Using the dev version (AT YOUR OWN RISK!)** … We do not recommend that you use the development version for analyses meant for publication. Please wait until we make a new version available on CRAN for that.

So anyone who took the bug from GitHub before December 2021 had to ask for this branch deliberately, having been told not to publish from it. That warning stood from 2017-01-09 (seven months before the bug landed) until it was removed on 2022-09-02.

Half of it expired in June 2021: CRAN archived the package, so "wait until we make a new version available on CRAN" pointed at a release that was never coming. Unaffected code remained obtainable, though: the archived 0.3.4.1 tarball is still downloadable from [CRAN's Archive](https://cran.r-project.org/src/contrib/Archive/rcicr/), and `install_github('rdotsch/rcicr', subdir = 'pkg')` returned the clean 0.3.4.2 on `master`. What was unavailable after June 2021 was a *maintained* clean version, not a clean version.

| you installed | you got | affected |
|---|---|---|
| before 2016-10-26 | **install fails** — no root `DESCRIPTION` yet | — |
| 2016-10-26 → 2017-08-15 | 0.4.0 | no |
| **2017-08-15 → 2021-09-23** | **0.4.0** | **yes** |
| 2021-09-23 → 2021-12-28 | 0.4.1 | **yes** |
| after 2021-12-28 | the branch became the default; see the table above | **yes** |

### From the latest GitHub release — `install_github('rdotsch/rcicr@*release')`

| you installed | you got | affected |
|---|---|---|
| before 2026-07-27 | **nothing** — no releases existed to resolve against | — |
| 2026-07-27 → 2026-07-28 | 1.1.0 | **yes** |
| 2026-07-28 → 2026-08-07 | 1.2.0, then 1.2.1 the same day | **yes** |
| 2026-08-07 → 2026-08-08 | 1.2.2 | **yes** |
| 2026-08-08 → 2026-08-18 | 1.2.3 | **yes** |
| 2026-08-18 → now | **1.3.0** | no |

Every GitHub release before 1.3.0 carries the bug. `v1.0.1` has a tag but no GitHub *release*, and it was tagged retroactively in 2026 against a 2023 commit, so `@*release` never resolved to it.

## How to check your own analysis

The affected call needs **both** things: a `participants` vector *and* `save_individual_cis = TRUE`. So the question to answer is:

> **Did a single `generateCI()` call produce one classification image per participant?**

If the answer is no, nothing here touches you. There are several ways to answer that question, and most of them work even if the output files themselves are long gone. Each of the headings below is one such way, starting from whatever you still have.

**You did not make per-participant images at all.** A group-level classification image is unaffected: it is a mean across participants and does not depend on their order.

**You made them with the batch functions.** `batchGenerateCI()` and `batchGenerateCI2IFC()` cannot reach the defect: they call `generateCI()` once per group with `participants = NA`, so the loop never runs, and they name each image from the group value they are iterating over. This is also the route the documentation has recommended since the CRAN era. `generateCI2IFC()` does not expose either argument.

Checked across the whole affected range, not only the latest release: `batchGenerateCI()` passes `participants = NA` explicitly in every version from 2017-08-15 to 1.2.3, and neither batch function mentions `save_individual_cis` in any of them. Measured too: run against the released 1.2.3 with participants `p1 … p12` in collection order, where the direct call gets 11 of 12 filenames wrong, both batch functions got 0 of 12 wrong and created no `individual_cis` directory at all.

**You still have the output.** The images live in a directory called `individual_cis`, and nothing else in the package writes there: the only such write is in `computeParticipantCIs()`, reachable only from `generateCI()`. That has held in every version: exactly one file contains the write, `R/generateCI.R` from 2017 through 1.2.3 and `R/ci-compute.R` after the loop was extracted. Individual files are named `ci_<participant>.png`, or `antici_<participant>.png` for an antiCI, so the naming outlives a renamed folder or moved contents. `filename` is a free-form argument, though, so a group CI written by `generateCI(save_as_png = TRUE, filename = "p3")` lands at `<targetpath>/ci_p3.png` and is indistinguishable by name alone (`R/generateCI.R:207`, `saveToImage()` at 392-406). Once the directory structure is gone, confirm a loose file from the script or the data before treating it as an individual CI.

**Telling one folder of PNGs from another.** Every other writer in the package puts its images *flat* into the path you gave it: `individual_cis/` is the only subdirectory any of them creates, which is what makes its presence such a clean signal. The filenames differ too:

| what wrote it | where | filename |
|---|---|---|
| `generateCI(save_individual_cis = TRUE)` | `<targetpath>/individual_cis/` | `ci_<participant>.png` — bare participant id |
| `generateCI(save_as_png = TRUE)`, `generateCI2IFC()` | `<targetpath>/` | `ci_<filename>.png`, defaulting to the base image name |
| `batchGenerateCI()`, `batchGenerateCI2IFC()` | `<targetpath>/` | `ci_<baseimage>_<by>_<unit>.png` |
| `autoscale(save_as_pngs = TRUE)` | `<targetpath>/` | `<name>_autoscaled.png` |
| `plotZmap()`, `generateCI(zmap = TRUE)` | `<zmaptargetpath>/` | `<filename>.png` |
| `generateStimuli2IFC()` | `<stimulus_path>/` | `<label>_<base>_<seed>_<NNNNN>_ori.png` and `_inv.png` |

An `antiCI` swaps the `ci_` prefix for `antici_` throughout. The composite names are the conclusive ones: `ci_face_participant_p3.png` came from a batch function and is fine, because a batch name always carries the base image and the grouping column, being built from them. A bare `ci_p3.png` is only a clue in the other direction: it is what the affected call writes, but `generateCI(save_as_png = TRUE, filename = "p3")` writes the same name for a group CI, one directory up, so confirm it against the script or the data before treating a loose file as an individual CI. Applying the recovery permutation to a set that is really something else would scramble output that was correct.

**You still have the script but not the output.** Read every `generateCI()` call that passes participant identifiers. Searching for the string `save_individual_cis` is *not* sufficient on its own: it is the sixth formal, so `generateCI(stim, resp, "face", rdata, pids, TRUE, ...)` sets it positionally, and `do.call()` with an assembled list hides it too. Search for `individual`, then read by hand any call passing six or more arguments positionally.

**You still have the raw data, but neither the output images nor the script.** Recompute with 1.3.0 and compare. This is the only fully conclusive answer, and it hands you corrected output as a side effect: the responses and the stimulus `.Rdata` are all that is needed, and a classification image is deterministic given them.

**You have only the published figure.** Then it cannot be verified directly, and the ordering question below is the best available evidence.

### If it was the affected call: was the ordering safe?

```r
identical(unique(participants), sort(unique(participants)))   # TRUE = the names were correct
```

"Sorted" means *lexically* sorted for text labels, per the table earlier: sorting a data frame by a `"p1"`-style column does not make it safe once you reach ten participants.

Recovery is a rename, not a recomputation: the file named `unique(participants)[i]` holds the image of `sort(unique(participants))[i]`.

## What it did to a second-stage analysis

The mislabelling scrambles the pairing between a classification image and the participant it came from, so any analysis joining a participant-level variable to those images was computed on a pairing it was not reported as. What that costs is measured rather than argued, in [`analyses/mislabelling-impact.md`](../analyses/mislabelling-impact.md) (five designs at 20,000 iterations per cell, generated by knitting [`analyses/mislabelling-impact.Rmd`](../analyses/mislabelling-impact.Rmd) so that the numbers in its prose cannot drift from the numbers in its tables).

The findings it establishes, kept separate rather than averaged, since nothing samples how common each design is:

- **A covariate unrelated to labelling order.** The association is lost, in proportion to how scrambled the ordering was; under the `p1 … pN` scheme that is close to total. The runs that then fail to reject are false negatives. This could lead to findings ending up in the file drawer: there may have been a real result, but it was very unlikely to be found in what the shuffle had turned into essentially random pairings.
- **The same, with no true association.** The false positive *rate* is unchanged, by exchangeability, which the iid null cells satisfy, and an outcome whose distribution drifts over collection order would not. Individual verdicts still flip in both directions, on roughly 9% of datasets.
- **Condition or covariate tracking labelling order.** The residual takes either sign: attenuated and correctly signed, or reversed.
- **A nuisance trend along collection order.** Contiguous blocks lose the effect; alternating assignment *can raise* hypothesis-consistent rejections above the correctly labelled analysis, in 9 of the 12 cells simulated. The three exceptions are all at N = 20, where that scheme's permutation happens to cut the contrast rather than add to it. That is a reminder that the direction belongs to the particular ordering, not to the design.

Two things hold across all of them. The permutation carries no reference to any hypothesis, so nothing in the mechanism aims it at a prediction, though how often it nonetheless agrees with one is not measured here. And none of these averages clears an individual analysis: the verdict on a single dataset can flip either way, so short of the check below it has to be re-run against corrected filenames whichever way it came out.

One thing does clear an analysis without a re-run: a property of the specific ordering, not something the averages above capture. If the permutation only exchanged participants sharing a value on the variable joined to the images, the join is unchanged and so is every number computed from it. The advisory gives readers `identical(predictor, predictor[misfiled])` for that. A `FALSE` leaves the question open; a `TRUE` means there is nothing to re-run.
