# Backlog → GitHub Issues: the proposed issue set

`BACKLOG.md` item 44. Every remaining open item and sub-item, **verified against the working
tree at `1.2.3.9000`** before being proposed — the point of the migration is to stop carrying
claims nobody has rechecked, so minting a stale issue would defeat it.

**Nothing has been opened.** This is the proposal for the maintainer to approve.

Items 13, 14 and 15 are umbrellas holding 14 independent tasks; they are split, because an
umbrella issue cannot be closed while any box is open and so becomes the same drifting surface
the file is. The four settled non-fixes went to `DECISIONS.md` in PR #173 and are not here.

## Summary

| | n |
|---|---|
| **Propose opening** | 21 |
| **Drop — already done, box never ticked** | 2 |
| **Drop — a decision, not work** | 1 |
| **Fold into another issue** | 1 |

---

## Drop — verified done, the box was simply never ticked

| source | claim | what the tree says |
|---|---|---|
| 15 | "Add `CONTRIBUTING.md` — issue #24" | The file exists (314 lines) and #24 is closed. |
| 15 | "Expand the vignette set beyond getting-started: batch/multi-participant analysis, choosing a scaling method, interpreting InfoVal and z-maps" | `reverse-correlation-walkthrough.Rmd` shipped in item 22 and covers all four: scaling ×12, z-maps ×10, participants ×6, batch ×5, infoVal ×3. |

## Drop — a decision, already recorded

**`styler`.** [`DECISIONS.md`](DECISIONS.md#if-the-package-is-ever-run-through-styler-it-goes-in-as-a-commit-of-its-own)
already says that if it ever happens it goes in as a commit of its own. An issue would invite
someone to do the thing the decision defers.

## Fold together

**Item 1 (CRAN archived) into item 20 (reinstatement).** Item 1's three boxes are "decide CRAN
vs GitHub-only" (decided: CRAN), "if GitHub-only, fix the README" (moot), and "if returning,
see item 20". Nothing is left in it that item 20 does not carry.

---

## Propose opening — 21 issues

### P1 — waiting on CRAN

| # | title | source |
|---|---|---|
| 1 | Log CRAN's reply to the 1.2.3 submission and update the README once the outcome is known | 20 + 1 |
| 2 | Publish the announcement post once the CRAN outcome is known | 21 |

### P2 — correctness and behaviour

| # | title | source | verified |
|---|---|---|---|
| 3 | A uniform base image silently becomes all-`NaN` | 31 | `(img - min(img)) / (max(img) - min(img))` at `generateStimuli2IFC.R:130` is 0/0 on a constant image; no guard |
| 4 | A decorated z-map below 256px dies in base R with `figure margins too large` | 33 | `zmapdecoration = TRUE` is still the default (`generateCI.R:87`) |
| 5 | The progress bar reports nothing under parallel execution | 13 | `pb` is created in the parent (`:168`) and ticked at `:230` and `:236`, both inside the `%dopar%` body — root cause of the closed #82 |
| 6 | `generateStimuli2IFC()` writes the loop counter `trial` into the saved `.Rdata` | 13 | present in the `save()` call at `:261` |
| 7 | Validate base images up front, naming the file and its actual dimensions | 14 | no such check; the failure is #124's `non-conformable arrays` |
| 8 | Name the writing package version in the `.Rdata` validation errors | 14 | the guards at `generateCI.R:146-158` already name the *missing field*; none names the version — **rescoped from the original box, which is otherwise done** |

### P2 — maintainability

| # | title | source | verified |
|---|---|---|---|
| 9 | `@import matlab` masks `base::sum()` with MATLAB semantics in six files | 13 | 6 files carry `@import matlab` |
| 10 | Add `lintr` as a config plus a CI job, gated at "no new lints" | 13 | absent from `.pre-commit-config.yaml` and from CI |
| 11 | Split `generateCI()` — 533 lines mixing CI computation, masking, scaling, PNG writing, z-maps and parallelism | 13 | 533 lines, up from the ~440 the backlog recorded |
| 12 | De-duplicate the mask-import logic between `generateCI()` and `plotZmap()` | 13 | `generateCI.R:409` and `plotZmap.R:52` |
| 13 | `raster` costs four packages and a C++ toolchain for three plotting calls | 34 | three `raster::plot()` calls (`plotZmap.R:125`, `:131`, `:172`); `raster` in `DESCRIPTION` Imports |

### P2 — documentation and onboarding

| # | title | source | verified |
|---|---|---|---|
| 14 | Publish a pkgdown site | 15 | no `_pkgdown.yml`, no `docs/` |
| 15 | Add a `CITATION.cff`, kept in sync with `inst/CITATION` | 15 | `inst/CITATION` exists; no `CITATION.cff` |

### P3 — smaller cleanups and questions

| # | title | source | verified |
|---|---|---|---|
| 16 | `generateStimuli2IFC()` leaves the user's RNG stream where it landed | 41 | `set.seed(seed)` at `:82`, never restored. **Named publicly** in the #52 close comment as tracked separately, so it should exist |
| 17 | Vendor a small, licence-clean example dataset for users | 15 | **Named publicly** in the #92 close comment as worth a fresh issue |
| 18 | Retire `ChangeLog` as a live file | 40 | still present, last substantive entry 2014 |
| 19 | Decide the fate of the superseded Medium link in `README.md` | 42 | `README.md:68`, the only URL that has ever 403'd a checker; already described there as superseded |
| 20 | The vignette's `par()` bookend restores what nothing changed | 43 | `old_par` at walkthrough `:38`, restored `:422` |
| 21 | Adopt tidyverse style throughout, as a v2 breaking change | 36 | explicitly post-v1 |

---

## Issue bodies

Each issue takes its `BACKLOG.md` section as its body, with three edits: the item number
dropped, any stale line corrected against the verification column above, and cross-references
to other items rewritten as issue links once the numbers are known. The prose is kept — it
carries the measurements and the rejected alternatives, which is the part that is expensive to
reconstruct.

## Afterwards

`BACKLOG.md` is deleted and `AGENTS.md`'s "Backlog" section becomes a pointer to the tracker,
keeping the guiding constraint it records — never change existing call syntax, argument
meanings, or numeric output silently — which belongs to the project, not to the file.
