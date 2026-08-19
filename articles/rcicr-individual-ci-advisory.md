# Individual-CI filename advisory (rcicr 0.4.0–1.2.3)

**This is an advisory, not a bug report. The bug is in versions 0.4.0
through 1.2.3 — 0.4.0 only if you installed it on or after 2017-08-15,
see “Which versions were affected” below — all of them GitHub-only, and
is fixed in
[1.3.0](https://github.com/rdotsch/rcicr/releases/tag/v1.3.0). It never
reached CRAN: the last CRAN release, 0.3.4.1, predates the option
involved by thirteen months, so nothing installed with
`install.packages('rcicr')` was ever affected.** It is published so that
anyone who finds it later can work out whether their own stored results
are affected.

## This is the bug

The bug is very specific to one function, which is
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md).
When called with `participants = ...` and `save_individual_cis = TRUE`,
it could save each participant’s classification image under a different
participant’s filename. The images themselves were computed correctly;
only the file names were wrong.

## Which versions were affected

It affects versions 0.4.1 through 1.2.3 outright, all of them
GitHub-only. **Version 0.4.0 is a special case**: it lived only on the
`development` branch, from 2016-10-26 to 2021-09-23, and the bug entered
partway through that window, on 2017-08-15 — so whether a 0.4.0 install
is affected depends on exactly when you got it, not the version string
alone (see
[`individual-ci-mislabelling.md`](https://github.com/rdotsch/rcicr/blob/main/notes/individual-ci-mislabelling.md)
for the install-date table). Any of these versions is only actually
affected whenever participant IDs weren’t already in sorted order
(e.g. “p10” sorting before “p2”, the ordinary case from the tenth
participant on).

However,
[`batchGenerateCI()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md),
[`batchGenerateCI2IFC()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md),
and
[`generateCI2IFC()`](https://rdotsch.github.io/rcicr/reference/generateCI2IFC.md),
the route most tutorials recommend for per-participant CIs, are not
affected by the bug. If that’s what your analysis used, this doesn’t
apply to you. The group classification image that
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
returns is also unaffected; only individual filenames from the direct
call could be wrong.

## How you can quickly check whether this affects you

If your output has an `individual_cis/` directory from a direct
[`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
call, run

``` r

identical(unique(participants), sort(unique(participants)))
```

using the participant vector in the order you originally passed it.
`TRUE` means the files were named correctly. If it’s `FALSE`, the images
are still correct; only the filenames are swapped, and the mapping to
fix them is exact (no need to recompute anything, see below).

## How to fix it

If you’re affected, the file named `unique(participants)[i]` actually
holds the classification image for `sort(unique(participants))[i]`, so
this recovers the correct filenames without recomputing anything (two
steps, so the rename doesn’t overwrite a file it hasn’t renamed yet):

``` r

targetpath <- "."  # set to whatever targetpath you originally passed to generateCI()

old <- unique(participants)         # names the files currently carry
new <- sort(unique(participants))   # participants they actually contain

dir <- file.path(targetpath, "individual_cis")
step1 <- file.rename(file.path(dir, paste0("ci_", old, ".png")),
                      file.path(dir, paste0("ci_", old, ".tmp")))
if (!all(step1)) stop("Rename failed for: ", paste(old[!step1], collapse = ", "))

step2 <- file.rename(file.path(dir, paste0("ci_", old, ".tmp")),
                      file.path(dir, paste0("ci_", new, ".png")))
if (!all(step2)) stop("Rename failed for: ", paste(old[!step2], collapse = ", "))
```

(This assumes [`sort()`](https://rdrr.io/r/base/sort.html) orders your
IDs the same way the original run did — always true on the same
machine/session, and true for numeric or lowercase-ASCII IDs like
`p1`/`p12` regardless. It can differ under a different locale if your
IDs mix case, accents, or punctuation; check
`Sys.getlocale("LC_COLLATE")` against the original run if you’re unsure,
or recover on the machine that produced the files.)

(use the `antici_` prefix instead of `ci_` if you generated anti-CIs).

If you plan to keep using the package, get the fixed version from
[GitHub](https://github.com/rdotsch/rcicr)
(`remotes::install_github("rdotsch/rcicr")`; a CRAN submission is in
progress but not live yet).

More detail:
[`individual-ci-mislabelling.md`](https://github.com/rdotsch/rcicr/blob/main/notes/individual-ci-mislabelling.md)
has the full technical background (every affected version, its dates,
and exactly which configurations trigger the bug), and [issue
\#267](https://github.com/rdotsch/rcicr/issues/267) is the pinned
advisory tracking this on the repository.

## Personal note

I wrote rcicr in 2014, then left academia in 2017 and stopped
maintaining it. I had other responsibilities and no time for something
that felt part of a previous life. Recently I decided to see if AI would
be able to maintain the package, and this week it found this bug. I’m
glad it did, and I consider it the right and scientific thing to do to
fix it and announce it publicly. I don’t have the time to do more than
that, so if you’re still active in this field and would like to maintain
this package, I’d gladly hand over the reins.
