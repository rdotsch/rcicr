# Which rcicr you had, and whether it carried the individual-CI bug

Reference for anyone checking whether a stored analysis is affected by the mislabelling fixed
in 1.3.0 (issue #261). The short answer is that **the version number alone does not tell you** —
for anything installed from GitHub before 2026, the date matters more than the version.

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
| **0.3.4.1** | **2016-07-13** | **CRAN — last release, archived 2021-06-08** | **no** |
| 0.4.0 | branch only, from 2016-10-26 | GitHub branch | **only from 2017-08-15** (see below) |
| 0.4.1 | branch only, from 2021-09-23 | GitHub branch | yes |
| 1.0.0 | branch only, from 2022-09-02 | GitHub branch | yes |
| 1.0.1 | commit 2023-01-13, tagged 2026-07-28 | GitHub tag | yes |
| 1.1.0 | 2026-07-27 | GitHub release | yes |
| 1.2.0 | 2026-07-28 | GitHub release | yes |
| 1.2.1 | 2026-07-28 | GitHub release | yes |
| 1.2.2 | 2026-08-07 | GitHub release | yes |
| 1.2.3 | 2026-08-08 | GitHub release | yes |
| **1.3.0** | **2026-08-18** | **GitHub release + CRAN submission** | **fixed** |

**`0.4.0` is the trap.** That version string sat on the development branch from 2016-10-26 to
2021-09-23 — five years — and the bug entered partway through, on 2017-08-15. Someone who
reports "I used rcicr 0.4.0" cannot be classified from the version alone; only the install date
separates clean from affected.

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
| 2021-06-08 → now | **nothing** — archived, the install fails | — |

**No CRAN user has ever received a mislabelled file**, at any point in the package's seven years
there. The last CRAN release predates the `save_individual_cis` option by thirteen months.

### From the GitHub default branch — `install_github('rdotsch/rcicr')`

This is the one that bites: it installs whatever is at the tip of the branch, so what you got
depends entirely on when you asked.

| you installed | `DESCRIPTION` said | affected |
|---|---|---|
| 2016-06-23 → 2016-09-18 | **install fails** — repository created, no package code yet | — |
| 2016-09-18 → 2016-10-26 | **install fails** — code imported from R-Forge, but under `pkg/`, so there is no `DESCRIPTION` at the root for `install_github()` to find | — |
| 2016-10-26 → 2017-08-15 | 0.4.0 | no |
| **2017-08-15 → 2021-09-23** | **0.4.0** | **yes** |
| 2021-09-23 → 2022-09-02 | 0.4.1 | yes |
| 2022-09-02 → 2023-01-13 | 1.0.0 | yes |
| 2023-01-13 → 2026-07-26 | 1.0.1 | yes |
| 2026-07-26 → 2026-08-17 | 1.x.y or 1.x.y.9000 | yes |
| 2026-08-17 → now | 1.2.3.9000, then 1.3.0 | **fixed** |

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

Two steps, neither needing a re-run. `NEWS.md` under "Reproducibility impact" has the full
version.

1. **Look for an `individual_cis/` directory in your output.** This is the definitive check,
   because that directory is created by nothing else in the package. If it is not there, the
   affected call never ran.

   Prefer this over searching the script. `save_individual_cis` is the **sixth** formal of
   `generateCI()`, so a call can set it positionally —
   `generateCI(stim, resp, "face", rdata, pids, TRUE, ...)` — and write those files without the
   argument name appearing anywhere. A `do.call()` with an assembled argument list hides it the
   same way. Grepping for the name and finding nothing does *not* clear an analysis; finding no
   `individual_cis/` directory does.

2. **If the directory is there**, check whether the `participants` vector was in lexical order:
   `identical(unique(participants), sort(unique(participants)))`. `TRUE` means the names were
   correct.

If the output is long gone and only the script survives, search it for `individual` rather than
the full argument name, and read any `generateCI()` call with six or more positional arguments
by hand.

Recovery is a rename, not a recomputation: the file named `unique(participants)[i]` holds the
image of `sort(unique(participants))[i]`.
