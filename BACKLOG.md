# rcicr modernization backlog

A prioritized backlog for bringing `rcicr` to a modern, maintainable state **without
breaking the API that researchers depend on**.

**Compiled:** 2026-07-26, against `main` @ `b6ab269` (v1.0.1).

> **State on 2026-08-07. CRAN reviewed 1.2.1 and asked for changes; 1.2.2 answers them and
> is released.** The reply was a request, not a rejection on the merits: seven points, all
> of them mechanical apart from one. `cran-comments.md` opens with a point-by-point
> response.
>
> **1.2.2 is tagged `v1.2.2` and published as a GitHub release**, with every external check
> in and recorded: win-builder 1 NOTE on both R-devel and R-release, R-hub `Status: OK` on
> all three platforms. **The one step left is the maintainer resubmitting to CRAN
> personally**, built from the tag rather than from `main` — `main` now carries `.9000`
> again, which a submitted tarball must never do.
>
> **The one that mattered: no function may write to a default path.** `stimulus_path`,
> `targetpath` and `zmaptargetpath` have lost their defaults (`./stimuli`, `./cis`,
> `./zmaps`) and are now required whenever the call actually writes. This is a **breaking
> change**, the first this package has made, and `NEWS.md` leads with it and the migration.
> It closes **item 24** for free. Numerically inert: the release gate reports 135 checks
> identical against v1.2.1, `max|d| = 0` everywhere.
>
> The rest: the `DESCRIPTION` description reworded and given two DOI references; all 11 bare
> `T`/`F` replaced (part of item 13); every `\dontrun{}` and `\donttest{}` removed so that
> **every example now runs**, in about nine seconds total; `plotZmap()` restores `par()`;
> the walkthrough vignette resets `par()` without `on.exit()`.
>
> **Item 34 stays held**, for the same reason as before: a CRAN request is answered by a
> version addressing what was asked, not one that also swaps the plotting backend under
> three rendered outputs. **Item 21** still announces availability, which is not yet true.
>
> Open items are **1, 20, 21, 25, 27, 30, 31, 33, 34, 40, 41, 42** in the table below, plus
> **13–15** and **36**, which are kept out of it.

**Last updated:** 2026-07-28 — P0 items 2–8, plus 9, 10, 11, 12, 16, 18, 19 and 22, fixed
and released as **v1.1.0**; see `NEWS.md`. All mechanical CRAN blockers are closed; what
remains in item 1 is the submission decision itself. **Items 23–25 were opened by a
test-intent audit on 2026-07-27**; 23 is now fixed (the mask is applied, and a second bug
in its boolean conversion was found and fixed with it), 24 is cosmetic, 25 is a logged
non-fix. **Items 28–31 were opened by a full source review on 2026-07-28**; 28 and 29 are
fixed, 30 is a documentation correction with optional follow-up, 31 is an unhandled edge
case. **Item 35 was opened by the first R-hub run on 2026-07-28** — the suite's only
non-portable assertion, failing on macOS — and is fixed and released as **v1.2.1**, which is
the tree submitted to CRAN on 2026-07-29. Items 13–15 are the only substantive untouched
work; they were held until after the submission and are now open.

**Reproducibility, verified 2026-07-28 — and re-checkable on demand.**
`tools/compare-release-output.R` installs v1.0.1 (`b6ab269`, the last release before 1.1.0) from its
own commit and runs it and the current tree through the same battery, then compares every
output. Coverage: 10 configurations across 64/128/512px, sinusoid and Gabor noise, `nscales`
3/5/6, `sigma` 10/25, one and two base images with shared and independent parameters,
`antiCI`, masks, non-contiguous stimulus subsets, and every analysis entry point —
`generateCI()`, `generateCI2IFC()`, `batchGenerateCI()` + `autoscale()`,
`computeCumulativeCICorrelation()`, both z-map methods and `computeInfoVal2IFC()`.

Result: noise basis, patch indices, RNG-drawn stimulus parameters, base images and the
stimulus PNGs written to disk are **bit-identical**; classification images and all four
scaling methods differ by at most **2.22e-16**, one ULP, which is the expected signature of
the documented `rowMeans(dims = 2)` change, and **0 pixels differ** once quantised to 8-bit.
Two deviations are on record as deliberate, both InfoVal at non-default `nscales`/`sigma`,
where v1.0.1 measured against the wrong null (item 2). The `.Rdata` contract held
append-only (gained `nscales`, `sigma`).

The gate is a **release blocker** — `CONTRIBUTING.md` → "Releasing" — and it earned that on
its first full run by catching item 32, a regression against v1.0.1 that no test in the
suite could have found.

Sources: the GitHub issue tracker (45 open issues), the published literature, and a
direct review of the codebase. Items marked **[verified]** were reproduced by running the
code, not inferred from reading it — each includes the observed behavior. Items marked
**[own review]** are not in the issue tracker; they came out of reading the source.

**Maintainer:** confirmed as **Ron Dotsch &lt;rdotsch@gmail.com&gt;** (`Authors@R`,
role `cre`). Git history shows the address was migrated from the old `r.dotsch@uu.nl`
university address to the current gmail one — relevant to item 1, since that migration
removes the stated cause of the CRAN archival.

## Guiding constraint: don't break existing analysis scripts

`rcicr` is used in published research, and analysis scripts are often re-run years later
to reproduce a figure or respond to review. So:

- **Never change the meaning of an existing argument, or the numeric output of an
  existing function, silently.** If a computation is wrong, fix it *and* make the change
  loud (warning + `NEWS.md` entry + version bump).
- **Deprecate, don't delete.** Keep old function names/arguments working with
  `lifecycle::deprecate_warn()` for at least one release cycle. `generateCI2IFC()` is
  already the model here: it's a thin wrapper kept purely for backward compatibility.
- **The `.Rdata` contract is sacred.** Every CI-computation function reads the `.Rdata`
  file written at stimulus-generation time. Old files must keep loading — that's what
  the `pre_0.3.0` / `generator_version` fields exist for. Add fields, never rename or
  remove them.
- **Anything that changes a published number** (most importantly `computeInfoVal2IFC()`)
  needs a documented migration note so researchers can tell which version produced a
  result.

Legend: **[P0]** correctness/blocking · **[P1]** high value · **[P2]** worthwhile ·
**[P3]** nice to have

## Suggested order of attack

**All seven P0 code bugs (items 2–8) are fixed, along with every P1 *code* item — the
dependency, toolchain, parallelism, test-coverage and vignette work — all released as
v1.1.0.** Items 1, 20 and 21 are not code: item 20's checklist has now been run twice and
**1.2.2 is tagged and released on GitHub**, so item 1 waits on the maintainer's resubmission
and then on CRAN's reply, and item 21 waits on the outcome. **Item 23 is fixed** — the `plotZmap()` mask is applied, under a "Behaviour change"
heading in `NEWS.md` — as is **item 32**, the `.Rdata` field that was capturing
`generateCI()`'s z-map `sigma`, caught by the release gate on its first full run. **Items 37,
38 and 39 are all fixed** — the submission unblocked them, and they went in that order: the
error message, the `load()` guards it exposed, then the tests that would have caught both.
The rest of the table is triage: items 27, 30, 31, 33 and 34. Items **13, 14 and 15** (modernize the
R code, better errors, docs and onboarding) are the only substantive work still untouched —
they are the backlog proper for after CRAN, and are deliberately kept out of the table.
**Item 36** (tidyverse style as a v2 breaking change) is out of the table too, and is not
scheduled at all — see "Beyond v1" at the end.

| # | Item | Why | Size |
|---|---|---|---|
| ~~9~~ | ~~Drop 13 unused `Imports`~~ | **Done** — the 13 unused declarations are gone; `DESCRIPTION` now imports 15 packages. Removed the install-failure risk and the `R CMD check` NOTE that would have blocked resubmission | S |
| ~~16~~ | ~~Pointless 1.5 GB array exported to every worker~~ | **Done** — the per-worker array copy is gone. Addresses, but does not fully close, issue #12: the `ncores == 1` path and `img_size` scaling remain, which is why #12 stays open | S |
| ~~10~~ | ~~Replace deprecated `progress_estimated()` / `rbernoulli()` / `citEntry()`~~ | **Done.** The `rbernoulli()` replacement had to preserve the random *stream*, not just the distribution — it is `runif(n) > (1 - p)`, and the obvious `rbinom(n, 1, p)` would have silently changed every infoVal computed from a given seed | S |
| 1 | CRAN archived | **1.2.1 was submitted 2026-07-29 and the review asked for seven changes; 1.2.2 answers them and is tagged `v1.2.2`.** All external checks are in and recorded in `cran-comments.md`. What remains is the maintainer resubmitting personally, from the tag's tarball — nothing else here is actionable until CRAN replies | M |
| ~~11~~ | ~~Cluster cleanup (`on.exit`), serial fallback~~ | **Done.** `on.exit` cleanup in #130; the `ncores == 1` serial fast path landed via `startBackend()`. Test suite 140s → 4s, and serial/parallel output verified bit-identical | M |
| ~~12~~ | ~~Widen test coverage (scaling methods, z-maps, `participants`)~~ | **Done** — suite at 180 tests, 0 skips. It found three real `plotZmap()` bugs in code that read fine, one of which made `zmapdecoration = FALSE` entirely dead since R 4.2 | M |
| ~~18~~ | ~~Codecov step fails for want of a token~~ | **Done** — `fail_ci_if_error: false`; a red `main` now means the package is broken | S |
| ~~19~~ | ~~Close the 8 issues already fixed in `main`~~ | **Done** — 7 closed, #12 commented and left open as only partly fixed. ~22 remaining issues still unswept | S |
| 20 | CRAN resubmission checklist | **Run twice now — for 1.2.1 and again for 1.2.2**, which is what the checklist is for. At 1.2.2: `--as-cran` 0 errors / 0 warnings / 2 expected NOTEs, win-builder 1 NOTE on both R-devel and R-release, R-hub `Status: OK` on all three platforms. Keep it: every further round of CRAN review means running it again | M |
| 21 | Announcement post | Drafted in `notes/`; hold until the CRAN outcome is known | S |
| ~~22~~ | ~~Move the Medium walkthrough into a vignette~~ | **Done** — `vignettes/reverse-correlation-walkthrough.Rmd`. It now executes at build time, which proved its own premise: two lines of the published tutorial had already stopped working | M |
| ~~23~~ | ~~`plotZmap(mask=)` validated then never applied~~ | **Done** — the mask is applied, under a "Behaviour change" heading in `NEWS.md`. It had never worked in any released version, so no published result depended on the old output. Two more bugs sat behind it, including a boolean conversion that set every cell `FALSE` | S |
| ~~24~~ | ~~`generateReferenceDistribution2IFC()` litters a `./stimuli` dir~~ | **Done in 1.2.2** — fell out of removing the default write paths CRAN objected to. The directory is now created only when something is written to it, and this caller writes nothing. Not merely cosmetic after all: writing to the working directory uninvited is a CRAN policy violation | S |
| 25 | InfoVal test oracle mirrors the implementation | **Deliberately left** — risk already covered by the hand-check against the erratum and the golden master. Logged so it is not mistaken for an independent check | S |
| ~~26~~ | ~~InfoVal's null is seeded by accident and cannot be varied~~ | **Done** — the determinism is now a documented guarantee with a comment guarding the `set.seed()` it rests on, and `response_seed` makes the null varyable on purpose. Default output verified byte-identical to before the change | S |
| 27 | `return_as_dataframe = TRUE` drops all but the first base image | **Documented, not fixed.** Correct under the default `use_same_parameters = TRUE`; silent only with `FALSE`. Widening the frame changes the return shape, so it needs a new argument rather than a redefinition — do it only if a user asks | S |
| ~~28~~ | ~~`generateCI()`'s single-trial 4096-parameter truncation is a no-op~~ | **Done** — the vector branch tested `length(params) == 4092` and then truncated to 4092, so it could never fire on the 4096-length input it existed for. Single-trial analysis of pre-0.3.0 files worked in no released version | S |
| ~~29~~ | ~~`autoscale()` aborts on masked classification images~~ | **Done** — a bare `range()` over the `NA`s that `generateCI(mask=)` writes by design. The single-CI `applyScaling()` path had always guarded this; only the batch path did not | S |
| 30 | InfoVal `ref_lookup` table empty since 2018 | **Open, triage.** Empty *correctly* — the erratum formula redefined the norms its rows summarised. Either repopulate (four measurements) or remove ~55 lines of matching machinery; the machinery is kept for now so repopulating stays cheap | S |
| 31 | A uniform base image silently becomes all-`NaN` | **Open.** `maximize_baseimage_contrast` computes 0/0 on a constant image and writes the `NaN` base into the `.Rdata` with no warning. Only bites synthetic or blank bases, but fails silently | S |
| ~~32~~ | ~~The `.Rdata`'s noise `sigma` overwrote `generateCI()`'s z-map blur `sigma`~~ | **Done** — `load()` assigns into the function's frame, so the `sigma` item 2 added to the file replaced the z-map argument of the same name from 1.1.0 on. Every argument is now kept across the `load()`. Caught by the release gate, not by the test suite | S |
| ~~35~~ | ~~`test-plotZmap.R:68` fails on macOS — blocked the CRAN submission~~ | **Done** — released as 1.2.1. The second distinct value was macOS quartz writing an **alpha channel** where cairo writes RGB, not an antialiasing artifact, so every option drafted before a macOS run was wrong. The count now ignores the alpha plane, and the value it compares is an ordering rather than a constant — the first fix pinned the background grey and macOS failed that too, at 0.573 vs 0.502. CI gained macOS and Windows runners, and the z-map is now pinned by the golden master on all three | S |
| 34 | `raster` costs 4 packages and a C++ toolchain for 3 plotting calls | **Open, post-CRAN.** `raster` → `terra` → GDAL/GEOS/PROJ is why R-hub's Linux and macOS jobs spent 30+ minutes installing dependencies while Windows, on binaries, took minutes. Used only by three `raster::plot()` calls in `plotZmap.R`. The release gate compares the z-map *matrix*, so it cannot catch a rendering regression here | M |
| 33 | A decorated z-map below 256px dies with `figure margins too large` | **Open.** `zmapdecoration = TRUE` is the default, so `generateCI(zmap = TRUE)` on a 128px stimulus set fails from inside base R, naming neither `rcicr` nor the cause. Not a regression — 256px and up are fine. Needs a clear early error, or a documented fallback | S |
| ~~37~~ | ~~Error paths are largely untested~~ | **Done** — `test-error-paths.R`, suite 291 → 323. Covers the length mismatch, every `.Rdata` "did not contain X" guard in both functions, all four mask-import failures, unreadable and non-square base images, and the `targetci` path that was never exercised. Each was confirmed to fire *its own* guard rather than an incidental error | M |
| ~~38~~ | ~~Two error messages paste the base image matrix, not its label~~ | **Done** — both sites now name `baseimage`. The defect was worse than logged: `paste0()` is vectorized, so it built one complete message *per pixel*, and `stop()` concatenated them — 1,024 copies and 8,190 characters at 32px, ~7 MB at 512px. Only error text changed; the gate reports this tree still reproduces v1.2.1 | S |
| ~~39~~ | ~~Two of the four `load(rdata)` sites did not guard their arguments~~ | **Done** — preventive, no live collision. `computeInfoVal2IFC()` restored 3 of its 5 arguments and `computeCumulativeCICorrelation()` none. The one that mattered: `target_ci` is read after a *second* `load()`, so a file carrying that name would have scored a different CI and returned a plausible number, not an error | S |
| 41 | `generateStimuli2IFC()` leaves the user's RNG stream where it landed | **Open, triage.** `set.seed(seed)` is never undone, so a script's next `runif()` differs depending on whether it generated stimuli first. Same family as CRAN's "do not change the user's state", though here the user asked, via a documented `seed` argument. Restoring `.Random.seed` on exit changes nothing this package computes — but it does change what a user's *next* draw returns, so it needs the gate run and a `NEWS.md` note | S |
| 40 | Retire `ChangeLog` as a live file | **Open, triage.** Not mandatory — R indexes `NEWS.md` for `news()` and ignores `ChangeLog` entirely — but it cannot simply be deleted: its 27 entries are the *only* record of 0.2.2 through 1.0.1. Freeze it as the pre-1.1.0 archive and drop the duplicated 1.1.0+ pointer entries, so it does one job and stops being a per-release chore | S |
| 42 | The superseded Medium link is the only URL that ever 403s a checker | **Open, triage.** `README.md:68`. Local `--as-cran` flags it; for 1.2.3 no external check did — clean on both win-builder runs and all three R-hub platforms — but win-builder *did* flag it on 1.2.1, so it tracks the checker's network, not the version. The README already calls the post superseded by the vignette. Not done now because `README.md` ships in the tarball, so removing it would force every external check to re-run | S |

Items 2, 3, 6 and 7 shared a shape worth remembering, because it will recur: **the
package failed silently or misleadingly rather than telling the user what went wrong.**
Item 13 (better errors) is the general form of that problem, and the fixes for 4 and 7
have started on it.

---

## P0 — Correctness and availability

> **Status: items 2–8 are FIXED** (see `NEWS.md` and
> `tests/testthat/test-fixed-bugs.R`, which holds a regression test for each).
> Item 1 (CRAN archival) remains open — it is a process/decision task, not a
> code fix. The golden-master test `test-regression-baseline.R` stayed green
> through every fix, which is the evidence that default-configuration results
> are unchanged.


### 1. Package is archived on CRAN — `install.packages('rcicr')` does not work
`rcicr` was **archived from CRAN on 2021-06-08**, reason: *"email to the maintainer was
undeliverable"* ([CRAN page](https://cran.r-project.org/package=rcicr)). The last CRAN
version is **0.3.4.1**, while GitHub is at **1.0.1** — so anyone following CRAN gets code
several years behind, and anyone following the README's `install.packages('rcicr')` gets
nothing at all.

This is the single highest-impact item: it affects every new user's first five minutes.

- [ ] Decide: return to CRAN, or commit to GitHub-only distribution?
- [ ] If GitHub-only: fix the README so it leads with `devtools::install_github()` and
      states plainly that CRAN is archived. *(Note: the README currently still says
      "Install the latest stable release from CRAN" — that instruction is broken today.)*
- [ ] If returning to CRAN: the stated archival reason (undeliverable maintainer email)
      **no longer applies** — `Authors@R` now carries `rdotsch@gmail.com`, replacing the
      old `r.dotsch@uu.nl` address that was presumably the one that bounced. Archived
      packages can be reinstated. **The concrete checklist is item 20 below**, built from
      an actual `R CMD check --as-cran` run rather than from guesswork.
- [x] Either way, add a `NEWS.md` so users can see what changed between 0.3.4.1 and 1.0.1.
      **Done** — `NEWS.md` now exists, including a "Reproducibility impact" section for
      researchers with published or in-progress results.

### 2. `nscales` and `sigma` are not saved in the `.Rdata` file → silently wrong InfoVal  ✅ **FIXED**
**Confirmed by reproduction.** Issue [#81](https://github.com/rdotsch/rcicr/issues/81).

`generateStimuli2IFC()` saves `base_face_files, base_faces, img_size, label, n_trials,
noise_type, p, seed, stimuli_params, stimulus_path, trial, use_same_parameters,
generator_version` — but **not `nscales` or `sigma`**.

This is not merely cosmetic. `generateReferenceDistribution2IFC()` re-generates the
stimulus set from the `.Rdata` file to build the InfoVal null distribution, and because
`nscales` isn't stored it silently falls back to the **default `nscales = 5`**. Verified:
stimuli generated with `nscales = 1` (12 noise parameters) produce a reference
distribution built from a `nscales = 5` basis (4092 parameters). **Anyone who used a
non-default `nscales` and reported an InfoVal has a number computed against the wrong
null distribution.**

- [x] Save `nscales` and `sigma` in the `.Rdata` file.
- [x] Have `generateReferenceDistribution2IFC()` read and forward them.
- [x] Fall back to the old default *with a loud warning* when loading a pre-fix `.Rdata`
      that lacks the fields (backward compatibility, but not silent).
- [x] Regression test asserting the regenerated basis matches the original.

### 3. `computeInfoVal2IFC(force_gen_ref_dist = TRUE)` silently does nothing  ✅ **FIXED**
**Confirmed by reproduction.** Issue [#113](https://github.com/rdotsch/rcicr/issues/113).

With a pre-existing `reference_norms` in the `.Rdata`, setting `force_gen_ref_dist = TRUE`
leaves the distribution completely untouched (verified: median identical before and
after). The flag short-circuits the *lookup-table* branch but never reaches the
regeneration branch, because `reference_norms` still `exists()` after `load()`.

The user gets no error — they simply believe they forced a recomputation when they did not.

- [x] Make `force_gen_ref_dist = TRUE` actually bypass the cached `reference_norms`.
- [x] Regression test.

### 4. `generateCI()` / `generateCI2IFC()` break on tibbles  ✅ **FIXED**
**Confirmed by reproduction.** Issues [#70](https://github.com/rdotsch/rcicr/issues/70)
and [#123](https://github.com/rdotsch/rcicr/issues/123) — the same root cause, reported
seven years apart, which suggests it keeps catching people.

Passing tibble columns yields `arguments must have same length` from
`aggregate.data.frame()`, because `tbl[, "col"]` stays a 1-column tibble where
`df[, "col"]` drops to a vector. Since `readr`/`dplyr` return tibbles by default, this is
now the *normal* path for a modern user, and the error message gives no hint of the cause.

- [x] Coerce inputs defensively at the top of `generateCI()`
      (e.g. `stimuli <- unlist(stimuli, use.names = FALSE)`).
- [x] Do the same in `batchGenerateCI()` / `computeCumulativeCICorrelation()`.
- [x] Tests covering `data.frame`, `tibble`, and bare-vector input.

### 5. `simulateNoiseIntensities()` is dead on arrival  ✅ **FIXED**
**Confirmed by reproduction.** Errors 100% of the time with
`object of type 'closure' is not subsettable`: it sizes its progress bar with
`data[, by]`, but neither `data` nor `by` is a parameter of the function (copy-pasted
from `batchGenerateCI()`); `data` resolves to `utils::data`, the function object. It also
ignores its own `img_size` argument, hardcoding `generateNoisePattern(img_size = 512)`.

- [x] Fix both bugs, or formally deprecate the function if it's no longer intended for use.
- [x] It is currently pinned by a `\dontrun{}` example and a test documenting the failure —
      update both when fixed.

### 6. The `mask` argument is completely unusable on R ≥ 4.2 **[verified] [own review]**  ✅ **FIXED**
Not in the issue tracker. `generateCI()` branches on `if (!is.na(mask))` at
`R/generateCI.R:178` and `:201`. When `mask` is a matrix — which is exactly what the
documentation tells users to pass — `is.na(mask)` returns a *matrix* of logicals, so the
`if` gets a condition of length > 1.

R 4.2.0 made that a hard error (it was previously a warning that silently used the first
element). Verified on R 4.3.3: passing a correctly-sized mask fails immediately with

```
Error: the condition has length > 1
```

So a documented feature of the package's main function is broken for every user on a
current R, and the error message points nowhere near the cause.

Behind it sits a **second** bug: `applyMask()` checks `!all(dim(mask_matrix) == 512)`
(`R/generateCI.R:315`) — a hardcoded `512` rather than `img_size`. So even after fixing
the `if`, any mask for a non-512px stimulus set is wrongly rejected, and the error message
misreports the expected dimensions as `img_size`.

- [x] Replace both `if (!is.na(mask))` guards with a scalar test
      (e.g. `if (!identical(mask, NA))` or `!is.null(mask)`).
- [x] Replace the hardcoded `512` with `img_size`.
- [x] Tests covering a matrix mask and a PNG-file mask at a non-512 `img_size`.
- [ ] Consider defaulting `mask` to `NULL` rather than `NA` (a cleaner sentinel), keeping
      `NA` accepted for backward compatibility. **Not taken:** the fix added a scalar
      `hasMask()` guard that accepts `NA`, `NULL`, a path, or a matrix, so changing the
      default would be churn for no user-visible gain. Still worth doing if the signature
      is ever revised.

### 7. `generateStimuli2IFC()` requires the base image to already be exactly `img_size` **[verified] [own review]**  ✅ **FIXED**
**Root cause of open issue [#124](https://github.com/rdotsch/rcicr/issues/124)**
("non-conformable arrays"), which has sat unexplained since 2024-09.

`generateStimuli2IFC()` validates only that the base image is *square*
(`R/generateStimuli2IFC.R:61`). It never checks that its dimensions match `img_size`,
because the resizing step is **commented out**:

```r
# Adjust size of base face
#base_faces[[base_face]] <- biOps::imgMedianShrink(img, x_scale=img_size/ncol(img), ...)
```

That line was disabled when the `biOps` dependency was dropped, and nothing replaced it.
So the base image is used at its original size while the noise is generated at `img_size`,
and they are added together at `R/generateStimuli2IFC.R:138`
(`combined <- (stimulus + base_faces[[base_face]]) / 2`) — which fails whenever the two
differ.

Reproduced: a 64x64 square base image with the default-style call `img_size = 32` gives
exactly the reported error:

```
Error: task 1 failed - "non-conformable arrays"
```

The reporter in #124 used a `.jpg` with the default `img_size = 512`; any base image that
isn't exactly 512x512 fails this way. The error surfaces from inside a `foreach` worker,
so it names neither the image nor the size mismatch — which is why it reads as inscrutable.

- [x] Validate up front that `dim(img) == c(img_size, img_size)`, with an error naming the
      file, its actual dimensions, and the requested `img_size`.
- [ ] Better: restore automatic resizing with a maintained package
      (`magick::image_resize()`, or `raster`/`terra`), so the original intent works again.
      Gate it behind an argument if silent resizing is considered too implicit.
      **Not taken:** the fix validates and errors clearly instead, which needs no new
      dependency and never silently alters a researcher's stimuli. Resizing remains the
      nicer end state.
- [x] Test covering a base image that doesn't match `img_size`.

### 8. Legacy `sinusoids`/`sinIdx` `.Rdata` support is unreachable  ✅ **FIXED**
`generateNoiseImage()` validates `length(params) != max(p$patchIdx)` **before** the block
that renames pre-0.3.3 `sinusoids`/`sinIdx` to `patches`/`patchIdx`. So an old-style `p`
hits `max(NULL)` → `-Inf` and always errors with "Stimulus generation aborted", despite
the code clearly intending to support it. This defeats the backward compatibility the
package explicitly promises for old `.Rdata` files.

- [x] Move the rename block above the validation.
- [x] Test with a genuine pre-0.3.3-shaped `p`.

### 23. `plotZmap(mask = ...)` is validated and then never applied **[verified] [own review]**  ✅ **FIXED**
Not in the issue tracker. Found 2026-07-27 while writing content tests for `plotZmap()`.

`R/plotZmap.R:34–69` does real work on the `mask` argument: reads it from a PNG if given a
path, checks it against the z-map's dimensions, verifies the values are 0/1 or TRUE/FALSE,
collapses identical RGB channels to one layer, and converts the result to boolean. Then the
function moves on to `zmap[abs(zmap) < threshold] <- NA` and **`mask` is never referenced
again** — it appears nowhere after line 69 (the mention at line 116 is inside a comment).

So the documented behaviour — *"If a cell evaluates to TRUE, the corresponding zmap pixel
will be masked"* — does nothing at all. A user who passes a correct mask gets a z-map with
no masking applied and no warning, which is the bad failure mode: the output looks
plausible and is silently wrong about which regions carry signal.

Blast radius is narrow. `generateCI()` calls `plotZmap()` **without** `mask`
(`R/generateCI.R:327`) and masks the CI itself via `applyMask()` beforehand, so z-maps
produced through the normal pipeline are unaffected. Only direct `plotZmap(mask = ...)`
calls are hit.

Why it survived the item-12 sweep: the only existing `mask` test asserts the *error* path
for mismatched dimensions, which exercises lines 34–69 and stops there. The tests added
alongside this entry check rendered content, but deliberately do not cover `mask` — see
below.

**Decided 2026-07-28: apply it.** The decision needed two facts that were not in hand when
this was written. First, `git log -S` traces the argument to commit `18e07cb` (2016-11-01),
whose message is *"add mask import and custom filename to plotZmap (todo: applying the
mask)"* — the import half was landed with an explicit todo that was never picked up, and two
weeks later the masking feature was implemented for `generateCI()` instead (PR #56, branch
`16-adding-masking-feature`). The pre-refactor monolith at `d93cb2b^` confirms the
application was never there to be lost in the 2017 file split. Second, GitHub code search
finds **no caller of `plotZmap()` outside this repository** — only a vendored copy in the
`rcpyci` Python port, which reimplemented `mask` as a translucent overlay, its author's
guess rather than evidence of intent. (Weak evidence: reverse-correlation analysis scripts
usually live in OSF supplements, not public repos.)

So the behaviour-change risk this entry was guarding against cannot materialise: no
published z-map has ever been masked, so no result depends on the current behaviour. The
only affected user is one passing a mask today and silently getting it ignored, whose intent
is unambiguously the masked version. Deprecating instead would have removed a documented
feature on the grounds that it was never built.

- [x] **Decided:** apply the mask. Landed as `zmap[mask] <- NA` after the threshold step.
- [x] Documented under a **"Behaviour change"** heading in `NEWS.md`.
- [x] Content tests added: a half-masked z-map must return that half to background and leave
      the other half painted, an all-zero mask must blank the map, matrix and PNG inputs must
      produce identical images, and one mask must remove the same half in `plotZmap()` as in
      `generateCI()`.
- [x] **A second bug found while fixing this.** The boolean conversion was
      `mask[mask == 0] <- TRUE` followed by `mask[mask == 1] <- FALSE` — a swap without a
      temporary. The first assignment coerces `TRUE` to `1`, so the second unsets exactly
      what the first set, and **every cell came out `FALSE` for any input**. Applying the
      mask without fixing this would have masked nothing. Verified in R before fixing, not
      inferred from reading.
- [x] **A third bug, and the one that nearly shipped inverted.** Both functions' roxygen said
      a *matrix* masks where the cell is `1`/`TRUE` while a *PNG* masks where it is black
      (`0`) — opposite conventions in a single sentence. The first version of this fix
      believed the docs and implemented the split, which would have made the same mask remove
      complementary halves in `plotZmap()` and `generateCI()`. `applyMask()` does
      `mask_matrix == 0` unconditionally (`R/generateCI.R`), verified by running it, so
      **`0`/black/`FALSE` is the masked region everywhere**. `plotZmap()` now matches, and
      both `@param mask` descriptions were corrected to the code rather than the reverse —
      users' existing masks were built against the behaviour. Caught by Ron asking whether
      the two input forms should really be opposite; the answer was in the sibling
      implementation, not in the documentation.
- [x] Checked the same argument in the sibling that accepts a mask. `generateCI()` applies it
      correctly, and now has a test covering the PNG import path and the masked region, which
      it previously lacked. The import logic remains duplicated between the two — issue #89,
      in item 13.

### Also found and fixed while working through the P0 items

None of these were in the original backlog; they surfaced only once the surrounding bug
was fixed and the code path actually ran. All three are covered by tests.

- **`write()` was creating a file called `data`.** `computeInfoVal2IFC()` called
  `write("Note that now that this simulated reference distribution has been saved...")`
  with no destination. Base R's `write()` defaults to `file = "data"`, so the message was
  never printed — it silently wrote a file named `data` into the user's working directory
  on every run that regenerated a reference distribution. Every other `write()` in the
  package passes `stdout()`; this one was missed. Invisible until the item 3 fix made the
  branch reachable.
- **`applyMask()` rejected masks it had just successfully converted.** In the
  RGB→greyscale branch, `stop()` sat *outside* the `else`, so it ran unconditionally —
  even when the channels were identical and the conversion had worked.
- **`applyMask()` compared against a hardcoded `512`** rather than `img_size`, so masks
  could never work at any other stimulus size, and the error message reported `img_size`,
  which is not in that function's scope.


---

## P1 — Modernize dependencies and toolchain

### 9. Remove 13 unused `Imports` (and shrink a 27-package dependency surface)  ✅ **FIXED**
`R CMD check` flags: `DBI`, `abind`, `assertthat`, `deldir`, `ggplot2`, `goftest`,
`gridExtra`, `iterators`, `munsell`, `plyr`, `polyclip`, `sp`, `tensor` — declared but
never used. These look like transitive dependencies of old `spatstat` that were pinned
directly and never cleaned up.

Every one is an install-failure risk for users and a CRAN-resubmission blocker. `sp` in
particular is in the retiring `sp`/`rgdal` ecosystem.

- [x] Delete all 13 from `Imports`. Low risk: nothing references them. **Done** — 27 to
      15, dropping `purrr` too once `rbernoulli` went. `scales` and `viridis` were kept:
      they are used via `::` and so never appear in `NAMESPACE`. This cleared both
      dependency `R CMD check` NOTEs.
- [ ] Consider `raster` → `terra` (raster is superseded and depends on `sp`).
- [ ] Move genuinely optional deps (`ggplot2`, `gridExtra` if reintroduced) to `Suggests`.

### 10. Replace deprecated function calls  ✅ **FIXED**
Verified deprecation warnings on current package versions:

| Call | Status | Used in |
|---|---|---|
| `dplyr::progress_estimated()` | deprecated in dplyr 1.0.0 | `batchGenerateCI`, `batchGenerateCI2IFC`, `computeCumulativeCICorrelation`, `generateReferenceDistribution`, `simulateNoiseIntensities` |
| `purrr::rbernoulli()` | deprecated in purrr 1.0.0 | `generateReferenceDistribution` |
| `citEntry()` / `citHeader()` | superseded by `bibentry()` | `inst/CITATION` |

These still work today but emit warnings on every run and will eventually be removed —
at which point five functions break at once.

- [x] Swap progress bars for `utils::txtProgressBar()` (already used elsewhere in the
      package, so no new dependency) or `cli`. **Done** in all five call sites.
- [x] Replace `purrr::rbernoulli(n, p)`. **Done — but NOT with `rbinom()`, despite what
      this item originally recommended.** `rbernoulli(n, p)` is literally
      `runif(n) > (1 - p)`; `rbinom(n, 1, p)` consumes the random stream differently.
      Verified across 150 seed/probability combinations that the two diverge, so swapping
      in `rbinom()` would have silently changed every reference distribution — and
      therefore every infoVal — computed from a given seed. Replaced with the `runif`
      form, verified bit-identical. **Anyone modernizing RNG calls in this package should
      check the stream, not just the distribution.** This did drop the `purrr` dependency.
- [x] Modernize `inst/CITATION` to `bibentry()`. **Done.** Also fixed a bug this exposed:
      the file read `meta$Date`, which broke when the stale `Date: 2023-01-13` field was
      dropped from DESCRIPTION — citations were rendering as "Ron Dotsch ()." with an
      empty year.

### 11. Make parallelism robust and user-controllable  ✅ **DONE**
Several related problems in one area:

- `generateReferenceDistribution2IFC()` **hardcodes** `ncores = parallel::detectCores() - 1`
  with no way to override. This already broke CI: under CRAN checks
  (`_R_CHECK_LIMIT_CORES_`) `makeCluster()` refuses >2 workers, so the test suite failed
  on the GitHub runner while passing locally. Currently worked around by mocking
  `detectCores()` in tests — the real fix is an `ncores` argument.
- Issue [#50](https://github.com/rdotsch/rcicr/issues/50): "closing unused connections"
  warnings — clusters not reliably stopped on error. Use `on.exit(parallel::stopCluster(cl))`.
- Issue [#66](https://github.com/rdotsch/rcicr/issues/66): no graceful fallback to serial
  execution on a single-core machine.
- `detectCores()` reports *physical* cores, not the cgroup/container limit — it
  over-subscribes on shared/HPC systems, which is exactly where researchers run big jobs.

- [x] Add `ncores` to `generateReferenceDistribution2IFC()` (default preserving current
      behavior). **Done** — this also removed the need for the `detectCores()` mock that
      the test suite had been using to stay under CRAN's 2-core check limit.
- [x] `on.exit()` cluster cleanup everywhere a cluster is created. **Done** — a
      `stopClusterSafely()` helper is registered via `on.exit()` at all three cluster
      sites, so workers are released even when an error interrupts the `foreach` loop.
      Verified: a mid-loop failure now leaks zero connections (issue #50).
- [x] **Skip the cluster entirely when `ncores == 1`. Done** via `startBackend()` in
      `R/zzz.R`, which registers `foreach::registerDoSEQ()` instead of building a
      one-worker PSOCK cluster. No loop body changed. `devtools::test()` went 140.1s →
      4.3s and `R CMD check`'s test phase `[8s/126s]` → `[8s/37s]`.
      Output is bit-identical, verified rather than assumed: neither parallel loop draws
      random numbers (parameters are drawn under `set.seed()` *before* the loop), so
      there is no per-worker RNG stream to diverge. Pinned by
      `tests/testthat/test-parallel-equivalence.R`.
**Known limitation, not planned.** `detectCores()` still reports physical cores rather
than a cgroup/container quota, so the default can over-subscribe on shared or HPC systems
— exactly where big jobs run. `default_ncores()` caps at 2 only under `R CMD check`.
Reading cgroup limits properly is platform-specific and not worth the complexity here;
users on such systems should pass `ncores` explicitly. Recorded so it is a known
trade-off rather than an oversight.

### 12. Fill out the test suite  ✅ **DONE**
A testthat suite now exists (**180 tests**) covering the pure functions, light I/O paths,
and an end-to-end smoke test. Every gap listed below is now closed — and closing the last
two turned up three previously unknown bugs in `plotZmap()`, which is the argument for
writing tests for code that "obviously works":

- [x] **Lock in the InfoVal formula.** Done — pinned in `test-regression-baseline.R`.
      See the note under item #17 — the implementation is
      currently correct per the published erratum, but *nothing tests it*, so a future
      refactor could silently regress a published metric.
- [x] **Cover the `scaling` methods. Done** — `test-scaling.R`, 15 assertions. Each method
      is tested against the property that defines it (`none` is the identity; `constant` is
      `(ci + k)/2k` and centres zero at mid grey; `matched` reproduces the base image's
      range; `independent` is `constant` with `k = max(abs(ci))` and touches exactly one
      boundary), plus the warning paths for a too-small constant and an unrecognised
      method. The load-bearing one is last: **scaling must change only `$scaled`, never
      `$ci`** — published numbers come from `$ci`, so a leak there would make results
      depend on a display choice.
- [x] Cover `mask` handling in `generateCI()` — **done**, now that item 6 is fixed; see
      `test-fixed-bugs.R`.
- [x] **Cover the `zmap = TRUE` paths. Done** — `test-generateCI-paths.R`. **Writing these
      found three real bugs, all now fixed** (see `NEWS.md`):
      1. `plotZmap()`'s undecorated branch used `if (bgimage != '')` on an image matrix —
         a length-`img_size^2` condition, an error on R >= 4.2. `generateCI()` always
         passes a background image, so `zmapdecoration = FALSE` was **entirely dead**.
         Same root cause as item 6; the decorated branch already used `identical()`.
      2. `plot.new()` ran *before* `par(mar = c(0,0,0,0))`, and it rejects a device too
         small for the current margins — so z-maps below roughly 100px failed with
         `figure margins too large`. Verified the ordering is the whole cause, and
         verified rendered output at 512px is **unchanged**, so this is not a visual
         change.
      3. `zmaptargetpath` was documented and accepted but never forwarded to `plotZmap()`,
         so every z-map went to `./zmaps` regardless.
- [x] **Cover `participants`. Done** — asserts the group CI equals the mean of
      independently computed per-participant CIs, that an **unbalanced** design makes
      per-participant averaging diverge from trial pooling (the reason the argument
      exists, and a check that the test is not vacuous), and that `save_individual_cis`
      writes one PNG per participant.

### 18. The Codecov CI step fails without a repository token  ✅ **FIXED**

`.github/workflows/test-coverage.yaml` runs `covr::package_coverage()` and uploads the
result to Codecov, which needs a `CODECOV_TOKEN` repository secret. That secret does not
exist, so the upload step fails. The coverage *computation* itself is fine — it is only
the upload that fails.

Note the exact trigger, because it is easy to misread: the step sets

```yaml
fail_ci_if_error: ${{ github.event_name != 'pull_request' || secrets.CODECOV_TOKEN }}
```

so on a **pull request** it is `false` and a missing token is tolerated, but on a **push
to `main`** it is `true` and the run goes red. That is why every PR so far has been green
while `main` shows a failing badge — the merge of #130 on 2026-07-27 is the current
example. Anyone judging the health of the repo from the front page sees a red X.

Three ways out were put to Ron, who chose the third:

- [ ] **Get a token.** Sign in to codecov.io with the GitHub account, add `rdotsch/rcicr`,
      copy the upload token into Settings → Secrets → Actions as `CODECOV_TOKEN`. Keeps
      the coverage badge and per-PR coverage comments. Codecov is free for public repos.
      *Still available later — see below.*
- [ ] **Or drop the upload.** Delete `.github/workflows/test-coverage.yaml` (and
      `codecov.yml`, plus its `.Rbuildignore` line). Coverage is deliberately partial
      here — I/O-heavy functions get lighter tests by design — so the reporting is of
      limited value, and `R CMD check` on release and devel is the check that actually
      matters.
- [x] **Keep the workflow but stop it failing.** `fail_ci_if_error: false`,
      unconditionally. Coverage is still computed and printed in the step log; only the
      upload is best-effort. Smallest change, and it keeps the token option open.

**Chosen and applied.** The conditional expression is gone; the setting is now a plain
`false` with a comment recording why, so nobody restores the old expression thinking it
was deliberate.

What this does *not* do: it does not make coverage reporting work. There is still no
Codecov token, so no badge, no per-PR coverage comments, and the `codecov.yml`
thresholds are inert. Coverage numbers are visible only in the workflow log. If that
reporting is wanted later, add the token and set this back to `true` — the first bullet
above is the recipe. The point of this change is narrower: a red `main` should mean the
package is broken, and it now does.

### 19. Eight issues fixed in `main` but still open on the tracker  ✅ **DONE**

Not a code task, but it is the largest gap between what the package *is* and what a
prospective user *sees*. As of 2026-07-27 the tracker has 30 open issues, and at least
these are already fixed:

| issue | what | fixed by |
|---|---|---|
| #70 | `generateCI()` cannot handle tibbles | #129 (P0 item 4) |
| #81 | `nscales` not saved in `.Rdata` | #129 (P0 item 2) |
| #113 | `force_gen_ref_dist` does nothing | #129 (P0 item 3) |
| #123 | `aggregate.data.frame` length error | #129 (P0 item 4, same cause as #70) |
| #124 | `non-conformable arrays` from `generateStimuli2IFC()` | #129 (P0 item 7) |
| #50 | "closing unused connections" warnings | #130 (`on.exit()` cluster cleanup) |
| #12 | large stimulus sets exhaust memory | #130 — **partially**, see below |
| #122 | `generateNoiseImage()` speed (a PR) | closed by #131 |

Why this matters more than it sounds: someone hitting the `non-conformable arrays` error
today searches the tracker, finds #124 open since September 2024 with no resolution, and
concludes the package is unmaintained — while the fix sits in `main`. The tracker is the
package's public health signal and it currently understates the state of the code badly.

- [x] **Done 2026-07-27.** #70, #123, #81, #113, #124, #50 each got a comment explaining
      the actual root cause — not just "fixed" — plus install instructions and a pointer to
      the `NEWS.md` reproducibility section. #122 closed automatically with #131.
- [x] **#12 commented, deliberately left OPEN.** `NEWS.md` says "addresses", not "fixes", and that
      wording is deliberate: #130 removed the ~1.5 GB array copied to every worker, which
      was the dominant cost, but the `ncores == 1` serial path in item 11 is still open and
      per-worker memory still scales with `img_size`. Comment with what changed and what
      did not, and leave it open, or close it explicitly scoped to the array copy.
- [ ] Sweep the remaining ~22 open issues the same way — some are likely stale or already
      resolved by earlier releases. This has not been done; the eight above are only the
      ones this modernization pass touched directly.

Best done as one batch after #131 merges, so #122 closes with it rather than by hand.

### 20. CRAN resubmission checklist  **[verified against `R CMD check --as-cran`]**

> **Submit the `v1.2.1` tag, and not the `v1.2.0` one.** The v1.2.0 tree fails `R CMD check`
> on macOS: `test-plotZmap.R` asserted properties of a rendered PNG that belong to the
> graphics device rather than to the drawing, and macOS quartz writes an alpha channel and
> applies colour management where cairo does not. Found 2026-07-28 by the first R-hub
> dispatch, which ran against the v1.2.0 tree. `tests/` is not in `.Rbuildignore`, so the
> tests ship and get run. The fix (item 35) landed on `main` **after** that tag, which is why
> 1.2.1 exists at all. CRAN's incoming checks are mainly Linux and Windows so it might not
> have blocked acceptance, but macOS binaries are checked on the farm after publication, and a
> package returning from archival should not arrive already red. See `DECISIONS.md` → Testing
> for the underlying rule.
>
> **All three external checks ran against the 1.2.1 tarball and are recorded in
> `cran-comments.md`** (2026-07-29): win-builder R-devel and R-release both `1 NOTE` — the
> incoming feasibility one — and R-hub `Status: OK` on Linux, Windows and macOS. **The
> submission was made the same day; this item is complete and the checklist below is kept
> because a rejection means running it again for the next version.**

**Current status — every box below is ticked.** The latest `--as-cran` run, on the
`rcicr_1.2.1.tar.gz` built at the repo root from the `v1.2.1` tree (2026-07-28), gives
**0 errors, 0 warnings, 2 NOTEs**, with `PDF version of manual ... OK`,
`HTML version of manual ... OK` and `checking for hidden files and directories ... OK` — the
last of those being the check the `^\.git$` `.Rbuildignore` entry exists for. The two
survivors are the archived/new-submission pair (carrying the Medium 403) and the sandbox
clock, both expected and both explained in `cran-comments.md`. `Version contains large
components` is gone and cannot return, because the tarball is built from a tag.
Nothing mechanical is left in this item beyond re-running win-builder and R-hub against the
1.2.1 tarball; the go/no-go and the submission itself are Ron's.

The run below is the earlier baseline, kept because it is what established the check
environment. **Re-run 2026-07-27 against `main` @ `5428b79`, in a check environment that
genuinely matches CRAN's.** `texlive` and `tidy` are installed, and the run set
`_R_CHECK_CRAN_INCOMING_=TRUE` / `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`, so the incoming
checks CRAN actually performs on submission ran for real.

Result: **2 NOTEs, no ERRORs, no WARNINGs.**

- `checking PDF version of manual ... OK`
- `checking HTML version of manual ... OK`
- `checking tests ... [8s/37s] OK`  (was `[8s/126s]`, see item 11)
- `checking examples with --run-donttest ... [12s/13s] OK`  (was `[15s/75s]`)

This supersedes an earlier run that reported 1 ERROR + 1 WARNING + 4 NOTEs. **Those were
all this sandbox, not the package** — no `pdflatex`, no `tidy`, and a leftover
`rcicr-manual.tex` from the failed PDF build. Installing the real toolchain was the only
change needed. Worth recording because an earlier note here claimed the ERROR/WARNING had
been *resolved* between runs; it had not — that run simply passed `--no-manual`, which
skips the check rather than passing it.

**The two remaining NOTEs:**

1. **`checking for future file timestamps`** — `unable to verify current time`. The
   sandbox cannot reach `worldclockapi.com`. Environmental; will not appear on CRAN.
2. **`checking CRAN incoming feasibility`** — the one that matters:

```
New submission
Package was archived on CRAN
Version contains large components (1.0.1.9000)
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2021-06-08 as email to the maintainer was undeliverable.
Found the following (possibly) invalid URLs:
  URL: https://medium.com/@rondotsch/...  Status: 403 Forbidden
    From: inst/doc/getting-started.html, README.md
```

Of these, only **"Version contains large components"** was ever an actual blocker, and the
1.1.0 bump closed it — see the first ticked box below; it no longer appears in the output.
"New submission" and "Package was archived" are expected for a reinstatement and are what
`cran-comments.md` exists to explain. The `codecov.io` URL that this NOTE previously also
flagged is gone.

Two details of the block above are pre-1.1.0 and should not be read as current. The
version string is now `1.1.0`. And the URL NOTE now reads `From: README.md` only — item 22
moved the walkthrough into a vignette, so `inst/doc/getting-started.html` no longer carries
the Medium link. That narrowing is what lets `cran-comments.md` tell a reviewer the content
ships *in the package* and the link is courtesy to nine years of citations.

#### Must do before submitting

- [x] **Bump the version to `1.1.0`. Done.** `1.0.1.9000` tripped "Version contains large
      components" — the `.9000` development suffix is not acceptable in a submission.
      `DESCRIPTION` now reads `1.1.0`, the `NEWS.md` heading is dated `2026-07-27`, and
      `ChangeLog` gained a 1.1.0 entry pointing at `NEWS.md` as the canonical changelog
      from this release onward. Minor rather than patch because some changes alter
      behaviour; not major because the public API is untouched.
      **This closes the last mechanical CRAN blocker** — what remains in this item is
      the go/no-go and the submission itself, which are Ron's.
- [x] **Fix the two flagged URLs — partly done.** The `codecov.io` one is **gone**: the
      badge rendered `unknown` (nothing has ever been uploaded, there being no token), so
      it was removed rather than repointed — a badge that reports nothing while looking
      like it reports something is worse than no badge. Re-verified by a second
      `--as-cran` run: it no longer appears.

      The **Medium link still flags, and that is fine.** It is out of `DESCRIPTION` now but
      remains in `README.md:51` and `vignettes/getting-started.Rmd:36`, where it is a
      genuinely useful pointer to the method walkthrough. **Do not delete it to silence the
      NOTE** — explain it in `cran-comments.md` instead.

      Be precise about *why* in that explanation. An earlier draft here said Medium "returns
      403 to non-browser user agents", which is wrong: retested with
      `curl -A "Mozilla/5.0"` and it still returns 403. Medium blocks by network origin
      (datacenter IPs), not by user agent, so spoofing a browser UA does not help and
      neither would any change on our side. The link resolves normally from a residential
      browser. Claiming the wrong cause to a CRAN reviewer who can check it is worse than
      claiming none.

      Item 22 (port the walkthrough into a vignette) is the durable fix: once the content
      lives in the package, the external link becomes a courtesy pointer rather than the
      only copy, and can be dropped from the vignette if a reviewer objects.
- [x] **Add `BugReports:` and a repo `URL:` to `DESCRIPTION`. Done.** There is currently no
      `BugReports` field at all, and `URL:` points only at the Medium article rather than
      the repository. Suggested:
      `URL: https://github.com/rdotsch/rcicr`,
      `BugReports: https://github.com/rdotsch/rcicr/issues`.
- [x] **Cap the core count under check. Done** via `default_ncores()` in `R/zzz.R`;
      verified to return `detectCores() - 1` normally and `2` when `_R_CHECK_LIMIT_CORES_`
      is set. Original text: `generateStimuli2IFC()`,
      `generateReferenceDistribution2IFC()` and `generateCI(n_cores=)` all default to
      `parallel::detectCores() - 1`. **CRAN policy allows at most 2 cores** in examples,
      tests and vignettes, and reviewers frequently object to `detectCores()` defaults on
      their own. The standard idiom keeps user-facing behaviour identical:
      ```r
      ncores = if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) 2L else parallel::detectCores() - 1
      ```
      Note this is *not* a behaviour change for researchers — only under `R CMD check`.
- [x] **Get the test time down. Partly done, and the diagnosis changed** —
      `skip_on_cran()` on three files. **Verified it fires**: under a real `R CMD check`
      the suite reports `SKIP 3 | PASS 119` (vs 145 in development).

      But the saving is only **148s → 126s**, not the large drop expected, because those
      three were not as dominant as assumed. The real finding is in the ratio:
      `checking tests ... [8s/126s]` — **8 seconds of CPU against 126 seconds elapsed.**
      Nearly all of it is waiting, not computing. `testthat.Rout` shows **22 PSOCK cluster
      spawns**, each starting a fresh R process that runs `library(rcicr)`. The same shape
      appears in `--run-donttest`: `[15s/75s]`.

      **This makes item 11 (the `ncores == 1` serial fast path) the single biggest lever on
      check time**, not a nicety — every test already passes `ncores = 1`, and every one of
      them still builds a one-worker cluster to run a sequential loop. Do item 11 before
      worrying further about check duration.
      Note the trap found while verifying this: `devtools::test()` and
      `testthat::test_local()` **set `NOT_CRAN=true` themselves**, so `skip_on_cran()` can
      never fire under them and they cannot be used to check that a skip works. Only a
      real `R CMD check` shows it. Original text: `checking tests` took **`[10s/148s]`**. CRAN wants the
      whole check comfortably under ~10 minutes on hardware slower than this. Put
      `skip_on_cran()` on the three slowest files — `test-recovery.R`,
      `test-smoke-pipeline.R`, `test-regression-baseline.R`. All three are development
      guards; none of them protects a CRAN *user*, and they keep running in GitHub CI.
- [x] **Guard the interactive prompt. Done.** Non-interactive callers now decline and
      regenerate with the requested iteration count rather than silently substituting a
      distribution built with a different one. Original text: `computeInfoVal2IFC()` calls `yesno::yesno()` at
      `R/computeInfoVal2IFC.R:118` with no `interactive()` check around it. In a
      non-interactive session — CRAN's checks, or anybody's batch script — that either
      hangs or errors instead of taking a sensible default.

#### Then, before hitting submit

- [x] **Check on platforms this machine cannot provide. Done 2026-07-29**, against the
      1.2.1 tarball: win-builder R-devel (r90304 ucrt) and R-release (4.6.1), `1 NOTE` each;
      R-hub Linux, Windows and macOS on R-devel, `Status: OK` on all three. CRAN tests those
      platforms and we do not.
- [x] **Write `cran-comments.md`. Done** — at the repo root, `.Rbuildignore`d. It states
      plainly that the 2021-06-08 archival was for an undeliverable maintainer address
      rather than any code or policy problem, that the address now works, and explains both
      remaining NOTEs. The win-builder and R-hub results are filled in as of 2026-07-29.
- [x] **Submit at <https://cran.r-project.org/submit.html>. Done 2026-07-29**, by Ron
      himself, from the 1.2.1 tarball — CRAN emails the maintainer address for
      confirmation, and that address working again is the entire point.
- [ ] Update `README.md` once the outcome is known.

### 21. Announcement post — drafted, waiting on the CRAN outcome

Ron asked for a write-up covering the important changes, the reproducibility guarantees,
and an honest note on how he used AI. **Drafts are written and live in `notes/`:**

- `notes/announcement-draft-medium.md` — long-form (~1800 words). Medium is the natural
  home: the existing Medium article is literally the `URL:` field in `DESCRIPTION`, so
  this reads as its sequel.
- `notes/announcement-draft-linkedin.md` — short teaser linking to the above.

Deliberately **not published yet** — hold until the CRAN question (item 1 / item 20) is
settled, so the post can say where the package actually lives instead of hedging.

- [ ] Ron to review both drafts. They are written in his voice and make claims about his
      work and his use of AI; they are a starting point, not a finished text.
- [ ] Verify every factual claim still holds at publication time (test counts, the 6x
      figure, dependency count). The drafts state specific numbers.
- [ ] Fill in the Medium link in the LinkedIn draft once published.

**One thing not to soften when editing:** the reproducibility section says default results
are unchanged *and* that two fixes genuinely change infoVal. Flattening that to "nothing
changed" would be both untrue and less useful to the researcher it is written for.

### 22. Move the Medium walkthrough into the package as a vignette  ✅ **DONE**

Ported as `vignettes/reverse-correlation-walkthrough.Rmd`. Covers the same ground as the
Medium post — install, generate stimuli, analyse, batch, scaling, online tasks, citation —
plus material the post did not have: what the four scaling methods actually do to an image,
how to read a z-map threshold, and a simulated observer with a known template so the
walkthrough shows a **recovered signal** (r = 0.55) rather than a grey smudge.

Every chunk executes at build time. Vignette rebuild is 15s for both vignettes together;
`--as-cran` still reports the same 2 NOTEs.

**Porting it immediately proved the premise.** Two lines of the published tutorial no longer
work against the current package: `autoscale(cis, saveasjpegs = TRUE)` (the argument is
`save_as_pngs`) and `install_github("rdotsch/rcicr", ref = "development")` (no such branch).
That is exactly the silent drift a vignette prevents, since a broken chunk now fails the
build.

Writing it also surfaced a real usability trap: after `batchGenerateCI()` the batch CIs
rendered as a nearly invisible overlay, because `autoscale()` rewrites `$scaled` and leaves
`$combined` alone. **That is intended** — Ron confirmed `$combined` is meant to survive the
call untouched so existing scripts keep plotting the same image — so it is now *documented*
in `?autoscale` and the vignette, and pinned by a test, rather than "fixed". An initial
change that rewrote `$combined` was reverted. Documentation that runs is a test.

- [x] Port the walkthrough into a vignette with executing chunks.
- [x] README and `getting-started.Rmd` now point at the vignette first; the Medium link is
      demoted to "further reading" rather than removed, so its inbound links keep working.
- [ ] **Ron:** add a pointer at the top of the Medium post to
      `vignette("reverse-correlation-walkthrough", package = "rcicr")`. Nothing in the repo
      can do this.
- [ ] The Medium URL still trips the CRAN URL NOTE (Medium blocks datacenter IPs). It is
      explained in `cran-comments.md`. If a reviewer objects, dropping it is now a one-line
      change in `README.md`, because the content no longer lives only there.

#### Original entry


The method walkthrough — the thing a new user is actually sent to read — lives at
`https://medium.com/@rondotsch/reverse-correlation-image-classification-using-r-a0701648fb0/`,
outside the repository. Three problems with that:

- It returns **403 to `R CMD check`** (Medium blocks non-browser agents), so it is a
  standing NOTE on every CRAN submission that has to be explained away each time.
- It cannot be versioned with the code. If an argument changes, the tutorial silently
  goes stale and nothing catches it — whereas a vignette is **executed** at build time,
  so a broken example fails the build.
- It is not available offline or via `vignette()`, which is where R users look first.

The plan was: port the walkthrough into a vignette with its chunks actually running;
keep the Medium post published (nine years of inbound links — a move of the canonical
copy, not a takedown); and demote the remaining `medium.com` references in `README.md`
and `vignettes/getting-started.Rmd` to further reading. All of that is done — see the
checklist at the top of this item, which is the live one.

---

## P2 — Usability and maintainability

### 13. Modernize the R code itself
- [x] **[own review]** Version metadata is inconsistent and stale. `generateStimuli2IFC()`
      hardcodes `generator_version <- '0.4.0'` (`R/generateStimuli2IFC.R:168`) into every
      `.Rdata` file, while `generateNoisePattern()` writes the real
      `utils::packageVersion('rcicr')`. So a file written today claims to come from 0.4.0
      — useless for the provenance purpose the field exists for, and it undercuts the
      `pre_0.3.0`-style compatibility checks that key off it. Covers issue
      [#29](https://github.com/rdotsch/rcicr/issues/29). Replace the literal with
      `utils::packageVersion('rcicr')`, but keep *reading* tolerant of old values.
      **Done 2026-07-27.** Nothing in the package reads the field, so "keep reads
      tolerant" had no code to change — the tolerance requirement is instead recorded at
      the write site and in `NEWS.md`, because it falls on *user* code: files in the wild
      say `'0.4.0'` regardless of what wrote them, and the field is now a
      `package_version` rather than a character string. The `pre_0.3.0` compatibility path
      does not in fact key off this field — it detects the old `sinusoids`/`sinIdx` layout
      structurally — which is just as well, given the field was lying.
- [ ] **[own review]** The `matlab` package exports its own `sum()` with MATLAB semantics
      (column sums for a matrix, not a single total), which masks `base::sum()` wherever
      `@import matlab` is in effect — six files. Package code is currently safe (its only
      two `sum()` calls are on vectors, where the two agree), but this is a silent trap
      for future edits, and it already bit the test suite. Prefer `@importFrom matlab ...`
      over `@import matlab`, importing only what is actually used.
- [x] ✅ **DONE in 1.2.2** — CRAN asked for it in the review of the 1.2.1 submission, so it
      was done there rather than waiting for this item. All 11 replaced; `tests/` and
      `vignettes/` were already clean.
      Replace bare `T`/`F` with `TRUE`/`FALSE` (11 occurrences across `autoscale.R`,
      `generateCI.R`, `generateStimuli2IFC.R`, `plotZmap.R`; two of them are *default
      argument values* in the public API — `generateCI(zmap = F, zmapdecoration = T)` and
      `plotZmap(decoration = T)`). **Re-checked 2026-07-29: in package code this is style,
      not the correctness hazard described here previously.** `T`/`F` are rebindable, but
      both defaults and function bodies resolve through the function frame → namespace →
      imports → `base`, which never reaches the user's global environment. A user
      redefining `T` breaks their own scripts, not this package. Still worth fixing —
      it is free, and the reasoning above has to be re-derived every time someone notices.
- [ ] **`lintr` — additive, and separable from `styler`.** It rewrites nothing, so it can
      land as a config plus a CI job with a zero-line diff to `R/`. Worth doing for the
      static analysis rather than the style: this codebase has form for bugs in code that
      reads fine (item 12 found three in `plotZmap()`), and `lintr` flags exactly that
      class — the bare `T`/`F` above, `seq_len()` vs `1:n`, unused variables, and the
      `if (matrix)` length-`n` conditions that caused items 6 and 23. Expect a large first
      report; take it as a to-fix list, not a gate, and set the gate at "no new lints".
- [ ] **`styler` — all-or-nothing, and still deliberately deferred.** It would reformat
      nearly every file in one sweep and destroy `git blame`. If it ever happens it goes in
      as a commit of its own, changing nothing else — see `DECISIONS.md`. Neither tool is
      in the pre-commit config, and that is why.
- [ ] `generateCI()` is ~440 lines mixing CI computation, masking, scaling, PNG writing,
      z-maps, and parallelism. Extract the internal helpers (`applyMask`, `applyScaling`,
      `combine`, `saveToImage` already exist) and cover them directly.
- [ ] De-duplicate mask-import logic — issue
      [#89](https://github.com/rdotsch/rcicr/issues/89) (it appears in both `generateCI()`
      and `plotZmap()`).
- [ ] **[own review]** The progress bar in `generateStimuli2IFC()` is created in the parent
      (`R/generateStimuli2IFC.R:113`) but ticked *inside* the `foreach` worker (`:163`), so
      it updates a copy in another process and the user sees no progress. This is the root
      cause of issue [#82](https://github.com/rdotsch/rcicr/issues/82). Use
      `doSNOW`-style progress handling, or report progress from the parent via `.combine`.
- [ ] **[own review]** `generateStimuli2IFC()` writes `trial` (a loop counter) into the
      saved `.Rdata` — harmless but meaningless, and it pollutes the documented file
      contract. Drop it when next touching the save call (item 2), noting it in `NEWS.md`.

### 14. Better errors for the most common user mistakes
The GitHub issues are dominated by confusing failure modes, not missing features:

- Issue [#124](https://github.com/rdotsch/rcicr/issues/124): `"non-conformable arrays"`
  from `generateStimuli2IFC()` — **root cause now identified, see item 7**; the underlying
  problem is a base-image size mismatch, but the message says nothing about that, and it
  surfaces from inside a parallel worker.
- Issue [#94](https://github.com/rdotsch/rcicr/issues/94): `object 'noise_type' not found`
  when an older `.Rdata` lacks the field.
- Issue [#123](https://github.com/rdotsch/rcicr/issues/123) / [#70](https://github.com/rdotsch/rcicr/issues/70):
  the tibble error above.

- [ ] Validate base images up front: exists, readable, square, expected size — with a
      message naming the offending file and its actual dimensions.
- [ ] Validate `.Rdata` contents on load, naming any missing fields and the package
      version that wrote the file.
- [x] Check `length(stimuli) == length(responses)` before `aggregate()`. **Done** as part
      of the item 4 fix.

### 15. Documentation and onboarding
- [ ] Publish a **pkgdown** site (docs are otherwise only readable via `?help` or an
      external Medium post).
- [ ] Expand the vignette set beyond getting-started: batch/multi-participant analysis,
      choosing a scaling method, interpreting InfoVal and z-maps.
- [ ] Add `CONTRIBUTING.md` — issue [#24](https://github.com/rdotsch/rcicr/issues/24).
- [ ] Vendor a small, licence-clean example dataset — issue
      [#92](https://github.com/rdotsch/rcicr/issues/92). Examples currently need a base
      face image users must supply, and the example data lives in a separate repo
      (`rcicr_examples`).
- [ ] Add a `CITATION.cff` so GitHub renders a "Cite this repository" button; keep it in
      sync with `inst/CITATION`.

### 24. `generateReferenceDistribution2IFC()` writes a stray `./stimuli` directory **[verified] [own review]**  ✅ **FIXED**
Fixed in 1.2.2, as a side effect of CRAN's request to remove default write paths.
`generateStimuli2IFC()` now creates its directory only when `save_as_png` or `save_rdata`
is set, and this caller passes neither — so it needs no path at all and writes nothing.
The `withr::with_dir()` workaround is gone from both affected examples, and
`test-required-paths.R` asserts the working directory is untouched. The empty
`tests/testthat/stimuli/` that a suite run used to leave behind was still present in the
checkout when this was fixed.


Found 2026-07-27. The function re-generates the stimulus set by calling
`generateStimuli2IFC()` (`R/generateReferenceDistribution.R:81`) **without forwarding a
`stimulus_path`**, so that call falls back to its own default and creates a `stimuli`
directory relative to whatever the working directory happens to be. Running the test suite
leaves a `tests/testthat/stimuli/` behind.

Harmless to results — `save_as_png = FALSE` and `save_rdata = FALSE` are passed, so the
directory ends up empty — but it litters the user's working directory, and it was invisible
here because `stimuli` is already in `.gitignore`, so it never appears in `git status`. The
function's own roxygen already documents the behaviour as a caveat (with a
`withr::with_dir()` workaround in the example), which is a note where a fix belongs.

- [ ] Forward a `stimulus_path` (a `tempdir()` default is fine — nothing is written to it)
      rather than documenting the side effect.
- [ ] Wrap the tests in `withr::local_dir()` so a suite run leaves no directory behind
      either way.

### 26. InfoVal's null is seeded by accident, and cannot be varied **[verified] [own review]**  ✅ **DONE**
Found 2026-07-27. `generateReferenceDistribution2IFC()` has no `seed` argument, yet its
output is fully deterministic: it rebuilds the stimuli via `generateStimuli2IFC()`, which
calls `set.seed(seed)` internally (`R/generateStimuli2IFC.R:53`) using the seed stored in
the `.Rdata` file, and that reset lands *before* the `runif()` draws in the simulation
loop. So the reference distribution — and every InfoVal computed from it — is fixed by the
stimulus file alone.

**Verified:** seeding 42 vs 99 around two calls gives byte-identical norms, and
`ncores = 1` vs `ncores = 2` also agree, so the null is portable across machines. Both are
now pinned by tests.

**The behaviour is worth keeping** — reproducibility is this package's guiding constraint,
and the Monte Carlo error it freezes in is small. Measured by splitting a 100,000-iteration
run into ten independent 10,000-iteration batches (`img_size = 32`, `n_trials = 300`):

| statistic | relative SD at `iter = 10000` |
|---|---|
| `median(reference_norms)` | 0.19% |
| `mad(reference_norms)` | 1.26% |

For a CI at infoVal ≈ 3 that is a spread of 2.95–3.06 (SD 0.036), i.e. below anything
interpretable. The existing `iter >= 10000` warning is well calibrated.

**The problem is that it is emergent, not designed.** Nobody chose to have the null inherit
the stimulus seed; it is a side effect of the stimulus rebuild. Two consequences:

- There is **no way to vary the null**, so a user cannot check the Monte Carlo error
  themselves — measuring the table above required batching one long run.
- Studies sharing `(seed, n_trials, img_size, nscales, noise_type)` share the *same draw*
  of the null, not merely the same null distribution, so their InfoVal errors are perfectly
  correlated rather than independent. For comparing InfoVals that is arguably a feature
  (one ruler), and it is evidently the intent behind the `ref_lookup` table in
  `computeInfoVal2IFC()` — but it means those errors cannot be treated as independent in a
  meta-analysis.

The real risk is silent breakage: anyone who moves or removes that `set.seed()`, or adds a
seed argument to `generateStimuli2IFC()`, changes every InfoVal ever computed without
touching `computeInfoVal2IFC()` at all.

**Done 2026-07-27.** All three boxes closed; default behaviour verified byte-identical
against norms captured before the change.

- [x] Document the determinism as an intended **guarantee** in
      `?generateReferenceDistribution2IFC`, and note that InfoVal is reproducible from the
      stimulus file alone.
- [x] Consider an explicit `seed` argument defaulting to the current behaviour, so the null
      *can* be varied deliberately. Purely additive — no change to the `.Rdata` contract or
      to any existing call. Shipped as **`response_seed`**, not `seed`: `seed` is an object
      in the `.Rdata` file, and since the function re-saves its own frame an argument of
      that name would overwrite the stimulus seed *and be written back*, corrupting the
      record of how the stimuli were generated. The seed is applied after the stimulus
      rebuild, not forwarded into it — forwarding would rebuild a different stimulus set, so
      the null would describe stimuli the participants never saw.
- [x] Add a comment at `R/generateStimuli2IFC.R:53` recording that the reference
      distribution depends on that `set.seed()` call, so it is not moved casually.

Also shipped: `save_rdata` and an invisible return on the generator, and a
`reference_norms_seed` provenance field in the `.Rdata` (append-only) so a file carrying a
varied null is distinguishable from one carrying the default.
`computeInfoVal2IFC(response_seed = )` forces regeneration — without that it would have been
silently ignored on every file that already had `reference_norms`, which is every file after
the first call, the same shape as item 23's dead `mask` argument — and never caches its
result.

### 25. `computeInfoVal2IFC`'s test oracle mirrors the implementation **[own review]**
Found 2026-07-27 during the test-intent audit. `test-computeInfoVal2IFC.R:24` recomputes
`(norm(ci, "f") - median(reference_norms)) / mad(reference_norms)` — the same expression as
the implementation. It therefore pins the *implementation*, not the published definition: a
formula that is wrong but consistently wrong in both places passes.

Deliberately **not** changed. The formula was hand-checked against the Schmitz et al.
erratum (item 17), and the golden master pins the resulting number, so the risk is already
covered from two other directions. Recorded here so that a future reader does not mistake
this test for the independent check it resembles — the genuinely independent oracle in the
suite is the one in `test-generateNoiseImage.R`.

- [ ] If ever revisited: assert against a worked example taken from the paper, or a
      hand-computed case with a fixed 2×2 CI and a known reference vector.

### 27. `return_as_dataframe = TRUE` silently drops all but the first base image **[verified] [own review]**
Found 2026-07-27 while fixing issue [#82](https://github.com/rdotsch/rcicr/issues/82). The
`return()` that hands back a trial's noise sits *inside* the per-base-image loop, so with
several base images it fires on the first one and the rest never run. The returned data
frame has one column per trial, not per trial × base image, so the shape cannot represent
them anyway.

With the default `use_same_parameters = TRUE` this is correct — every base image shares one
parameter set, so one noise image per trial is all there is. With
`use_same_parameters = FALSE` the caller silently receives only the first base image's
noise.

**InfoVal is not affected, checked rather than assumed.**
`generateReferenceDistribution2IFC()` is the only in-package caller of this path. It does
not pass `use_same_parameters`, so the rebuild always uses the default, and the first base
image's parameters are drawn from the same leading block of the RNG stream either way —
measured identical, max absolute difference 0 across both settings with two base images.

- [x] Document the restriction on `@param return_as_dataframe`: one noise image per trial,
      meaningful only under `use_same_parameters = TRUE`. **Done 2026-07-28**, with a
      `NEWS.md` entry under Documentation. Behaviour unchanged.
- [ ] Only if a user actually needs it: widen the returned frame to trial × base image.
      That changes the return shape, so it needs a new argument rather than a redefinition
      — the `.Rdata`/return contract is append-only.

### 16. Memory ceiling on large stimulus sets — concrete root cause found **[own review]**  ✅ **FIXED**
Issue [#12](https://github.com/rdotsch/rcicr/issues/12) reports that large stimulus sets
exhaust memory, with the suggested fix being "distribute stimulus matrices over multiple
.RData files". Reading the code suggests that's treating a symptom — there's a specific
and much cheaper fix.

`generateStimuli2IFC()` pre-allocates the **entire** stimulus array before the parallel
loop:

```r
stimuli <- matlab::zeros(img_size, img_size, n_trials)   # allocated up front
cl <- parallel::makeCluster(ncores, outfile = "")
...
stims <- foreach::foreach(trial = 1:n_trials, .packages = 'rcicr', ...) %dopar% {
  stimuli[,,trial] <- generateNoiseImage(...)            # writes to a worker-local copy
```

Two compounding problems:

1. At default settings that array is **1.5 GB** (512 × 512 × 770 × 8 bytes), and because
   it exists in the parent environment before `makeCluster()`, `foreach` exports a full
   copy to **every worker** — roughly 10 GB of pure overhead on an 8-core machine.
2. Each worker writes exactly one slice into its own private copy, which is then
   **discarded** when the worker returns. The assignment does nothing useful; only the
   local `stimulus`/`combined` values derived from it are actually used.

So the array is both enormous and pointless. Replacing `stimuli[,,trial] <- ...` with a
plain local variable should remove nearly all of the memory pressure without changing
behavior or the `.Rdata` contract.

- [x] Replace the pre-allocated array with a per-iteration local variable. **Done.**
- [ ] Re-measure the ceiling afterwards and document it, so users can size jobs up front.
- [ ] Only if still needed: chunking, or regenerating noise from `seed` on demand rather
      than storing every parameter vector.

### 28. `generateCI()`'s single-trial 4096-parameter truncation is a no-op **[verified] [own review]**  ✅ **FIXED**
Found 2026-07-28 in a full source review. `R/generateCI.R` truncates pre-0.3.0 parameter
matrices from 4096 columns to 4092. The matrix branch tests `ncol(params) == 4096`; the
vector branch, added by `4a6c58a` ("make generateCI work with a single input trial"), tested
`length(params) == 4092` and then truncated to `1:4092` — **a no-op that can never fire on
the 4096-length input it exists for**. Reproduced: with a 4096-parameter `.Rdata`,
multi-trial returns a CI and single-trial dies with `Stimulus generation aborted: number of
parameters doesn't equal number of patches!`.

Narrow (pre-0.3.0 files, single-trial CIs) but the failure is total, and it was present in
1.0.1 and every CRAN release. Fixed with a regression test that also pins the multi-trial
branch, so the two cannot be inverted.

### 29. `autoscale()` aborts on masked classification images **[verified] [own review]**  ✅ **FIXED**
Found 2026-07-28 in the same review. `generateCI(mask = ...)` sets masked pixels to `NA` by
design. `autoscale()` called a bare `range()` on `$ci`, so the scaling constant became `NA`
and the next line failed with `missing value where TRUE/FALSE needed`.

The telling detail: `applyScaling()` guards **every** reduction with `[!is.na(...)]`, so the
single-CI scaling path has always handled masked CIs and only the batch path did not.
Fixed with `na.rm = TRUE` plus an explicit error for a CI that is entirely `NA`, which would
otherwise produce an infinite constant.

### 30. The InfoVal reference lookup table has been empty since 2018 **[verified] [own review]**
`R/computeInfoVal2IFC.R` carries a `ref_lookup` tibble whose four rows are all commented
out. That was correct and deliberate: `01e547e` (2018-07-31) adopted the Euclidean norm and
scaling factor *k* from the erratum to Schmitz et al. (2019), which redefined the norms
those medians and MADs summarise. Reusing them would score CIs against a null built from the
wrong formula.

The consequence is that roughly 55 lines of lookup, matching and interactive-prompt code run
against an empty table and can never match, and the reference distribution is always
regenerated — correct, just slow. **`CLAUDE.md` described this as a working cache**, which
has been corrected. The machinery is deliberately kept rather than deleted so that
repopulating is a matter of measuring four numbers.

- [ ] Optional: re-measure `median(reference_norms)` and `mad(reference_norms)` under the
      current formula for seed 1, 512px, 10000 iterations at 100/300/500/1000 trials, and
      uncomment the rows. Each is one `generateReferenceDistribution2IFC()` run.
- [ ] If they are not going to be re-measured, delete the matching and prompt machinery
      instead — but not both halves independently, or the feature becomes unrecoverable.

### 31. A uniform base image silently becomes all-`NaN` **[verified] [own review]**
`maximize_baseimage_contrast = TRUE` (the default) computes `(img - min(img)) / (max(img) -
min(img))` in `R/generateStimuli2IFC.R`. On a constant image that is 0/0. Reproduced: the
`NaN` base image is written into the `.Rdata` with no error and no warning, and every CI
computed from that stimulus set inherits it.

A photograph is never uniform, so this only bites synthetic or accidentally-blank base
images — but it fails silently, which is the worst mode, and the symptom (all-`NaN` CIs)
appears far from the cause.

- [ ] Error when `max(img) == min(img)`, naming the file: a base image with no contrast
      cannot be used, and saying so at generation time costs nothing.

### 32. The `.Rdata`'s noise `sigma` overwrote `generateCI()`'s z-map blur `sigma` **[verified]**  ✅ **FIXED**
Found 2026-07-28 by `tools/compare-release-output.R` on its first full run — the first bug
the release gate caught, and one no test in the suite could have found, because it is a
regression *relative to 1.0.1* introduced by a change that was itself correct.

`generateCI()` reads the stimulus set with `load(rdata)`, which assigns into the function's
own frame. Item 2's fix started saving the noise `sigma` into the `.Rdata` — and
`generateCI()` has an argument of the same name, the z-map blur width. So from 1.1.0 on, the
saved value silently replaced the argument: z-maps were blurred with 25 instead of 3, and an
explicitly passed `sigma` did nothing. Measured at 512px: the number of pixels surviving the
`threshold = 3` cut moved from 1,157 to 2,439, and the map lost its entire negative tail
(range `-3.32 .. 4.07` became `3.00 .. 3.88`).

Only z-maps are affected — `sigma` is used for nothing else — and only for stimulus sets
generated with 1.1.0, which was a GitHub tag for about a day and never on CRAN. **Real-world
impact is therefore close to zero**, and `NEWS.md` says so rather than alarming anyone; the
value of the find is that the gate caught a live regression on its first run, in a code path
no test covered. Fixed by keeping copies of every argument across the `load()`, the
same guard `generateReferenceDistribution2IFC()` already carries, so a field added to the
`.Rdata` later cannot capture another argument. Regression test in `test-fixed-bugs.R`;
`NEWS.md` carries the reproducibility note.

### 33. A decorated z-map below 256px dies in base R with no useful message **[verified] [own review]**
Found 2026-07-28 while extending the release gate's battery. Measured on 1.1.0.9000:

| `img_size` | `zmapdecoration = TRUE` | `zmapdecoration = FALSE` |
|---|---|---|
| 128 | **`figure margins too large`** | ok |
| 256 | ok | ok |
| 512 | ok | ok |

`zmapdecoration = TRUE` is the **default**, so `generateCI(zmap = TRUE)` on a 128px stimulus
set fails outright, from `plot.new()` inside `.rasterImagePlot()`, with an error that names
neither `rcicr` nor the actual cause. Not a regression: v1.1.0 behaves identically, and
v1.0.1 was worse (it could not render a small z-map either way — that half was fixed in
1.1.0). Small stimulus sizes are common in simulations and tests, which is where this bites.

The margins the decoration needs genuinely do not fit on a 128px device, so the fix is not to
force it. Options, in order of how much they presume:

- [ ] Error early with a message naming the function, the size and the way out
      (`zmapdecoration = FALSE`, or a larger `size`).
- [ ] Or fall back to undecorated with a warning — silently changing what gets rendered,
      which is a behaviour change and needs Ron's call.

### 35. `test-plotZmap.R` fails on macOS — blocked the CRAN submission **[verified]**  ✅ **FIXED**
Fixed in #151 (`e0b74a7`) and released as **1.2.1**; the submitted tarball is built from
`v1.2.1`, not `v1.2.0`. What the fix turned out to be, and what it cost to get right:

- **The cause was the alpha channel, not antialiasing.** cairo (Linux, Windows) writes RGB;
  macOS quartz writes RGBA. An opaque alpha plane is a solid block of `1`s, so it adds a
  second distinct value to the array **without a pixel having been painted**. The count now
  runs over colour channels only, handling 1-, 2-, 3- and 4-channel layouts, verified against
  synthetic images in each: unpainted still collapses to one value, painted still yields two,
  so no detection power was given up.
- **The first fix over-tightened and macOS caught that too.** It also pinned the flat value to
  the background grey — which fails there, because quartz renders a `0.5` background at ~0.573
  where cairo gives 0.502 (colour management). That repeated the mistake being corrected, one
  line below it. The assertion is now an *ordering* — the same render over a darker background
  must come out darker — which survives any monotone transfer function.
- **The general rule, now in `DECISIONS.md`:** when a test reads pixels back from a graphics
  device, *every absolute property of those pixels belongs to the device* — channel count and
  value alike. Only relationships between renders are portable.
- **The class of failure cannot recur unseen.** `R-CMD-check.yaml` gained macOS and Windows
  runners (the matrix had varied only the R version against a hardcoded `ubuntu-latest`, so
  what read as two-platform coverage was one platform twice), and
  `test-regression-baseline.R` now pins z-map values for both `zmapmethod`s. The literals were
  captured on Linux and passed unchanged on macOS ARM64 and Windows — first direct evidence
  that the z-map arithmetic agrees across platform and architecture.

The original report follows.

Found 2026-07-28 by the **first R-hub run** (`gh run view 30370708485`, branch
`check/rhub-v1.2.0`, i.e. the `v1.2.0` tree). `linux (R-devel)` and `windows (R-devel)`
pass; `macos (R-devel)` returns `Status: 1 ERROR`:

```
── Failure ('test-plotZmap.R:68:3'): sub-threshold z-scores are not painted and supra-threshold ones are ──
Expected `unique(as.vector(below))` to have length 1.
Actual length: 2.
[ FAIL 1 | WARN 0 | SKIP 5 | PASS 220 ]
```

**This was a blocker.** CRAN checks on macOS, and a test ERROR there fails the submission —
so it had to be resolved before the tarball went in, even though nothing was wrong with the
package itself.

The failing assertion is the uniform-background trick from `CLAUDE.md`'s test conventions:
render a z-map that is entirely below threshold onto a flat grey background, and "nothing was
painted" becomes `expect_length(unique(as.vector(below)), 1)`. That is exact-equality over
every pixel of a rendered PNG, which is **not** a portable thing to assert. The default `png()`
device differs by platform — quartz on macOS (`capabilities("aqua")`), cairo on this box and
on the Linux runners — and they do not agree pixel-for-pixel on `rasterImage()` interpolation
and antialiasing. Two distinct values instead of one is consistent with an edge artifact of a
few pixels, not with the z-map being painted.

Unverified from here: the failure is macOS-only and cannot be reproduced on this Linux
container, so the exact pixel counts are unknown. **Confirm before fixing** — the artifact
(`gh run download 30370708485`) has the `.Rcheck` directory but not the rendered PNGs, since
the test writes to a `withr::local_tempdir()`. Print `table(as.vector(below))` from a macOS
run first, or add it to the failure message.

Options, cheapest first — **none of these was what the fix turned out to need**, because all
three assumed an edge artifact and the cause was a whole extra channel:
- [ ] **Assert the interior, not the whole frame** — exclude a one-pixel border and keep the
      exact comparison. Preserves the test's intent; costs nothing if the artifact is at the
      edge, which is what the two-value result suggests.
- [ ] **Assert a dominant value rather than a unique one** — e.g. ≥ 99% of pixels equal the
      modal value. Honest, but weaker, and it would no longer fail if a faint z-map were
      painted over the whole frame.
- [ ] **Pin the device** — `png(type = "cairo")` inside the test. Rejected unless the above
      fail: it tests a rendering path users on macOS do not have, which is the opposite of
      what a portability failure is telling us.
- [x] **Drop the alpha plane before counting** — what was actually done, once a macOS run
      showed the second value was a solid block of `1`s. Note the instruction above to get
      `table(as.vector(below))` from a macOS run before fixing was the right call: every
      option drafted from here was wrong.

Note the same file's other content tests (`the threshold is applied per pixel`, and the mask
tests from item 23) compare *renders against renders* and so passed — they are device-agnostic
by construction. Only line 68 asserts absolute uniformity. That is the pattern to copy.

Introduced by #138 (`887aea4`) and extended by #145 (`b2b0d9a`); it had never run on macOS
before, because CI was `ubuntu-latest` only. That separate question — whether
`R-CMD-check.yaml` should gain a macOS runner — was answered yes in the same PR: the whole
value of this R-hub run was catching what a single-platform CI cannot.

### 34. `raster` costs four packages and a C++ toolchain for three plotting calls **[verified] [own review]**
Found 2026-07-28, from watching R-hub: `windows (R-devel)` finished in minutes while
`linux (R-devel)` and `macos (R-devel)` sat in dependency installation for over half an hour.
Windows had CRAN binaries; the others built from source, and the expensive branch of the tree
is `raster` → `terra`, which compiles against GDAL/GEOS/PROJ.

`raster` is used in exactly three lines, all in `R/plotZmap.R` — `raster::plot()` at
[106](R/plotZmap.R#L106), [112](R/plotZmap.R#L112) and [147](R/plotZmap.R#L147), each on a
`raster(zmap)` built from a plain matrix. Nothing else in the package touches it. Dropping it
removes four packages nothing else here needs:

| removed | why it is expensive |
|---|---|
| `raster` | the API actually used |
| `terra` | C++ against GDAL/GEOS/PROJ — the bulk of the build time |
| `sp`, `Rcpp` | pulled in behind them; needed by no other `Imports` entry |

**A third argument, measured 2026-07-29: it is why every Linux CI job is at the mercy of an
apt mirror.** `terra` needs GDAL/GEOS/PROJ as *system* libraries, so
`r-lib/actions/setup-r-dependencies` fetches **333 apt packages** — `libgdal-dev`,
`gdal-bin`, `gdal-data`, `gdal-plugins`, `libgdal-grass`, `libgdal34t64` and their
dependencies — before `R CMD check` can begin. On a fast mirror that is ~6 minutes; on
2026-07-29 a slow `azure.archive.ubuntu.com` turned the same job into **35m23s**, with
individual gaps of 321s, 129s and 92s between `Get:` lines. An identical job on another PR
finished in 6m18s *inside that window*, so the variance is external — but the exposure is
ours, and it exists for three plotting calls. It also cost real diagnosis time: the stall
looked like a hung job and was called one before the log was readable.

**A second argument, found 2026-07-29 and not previously recorded: this is also the largest
security surface in the tree.** Of the 53 non-base packages rcicr pulls in recursively, 29
need compilation, and their vulnerabilities are not in R code — they belong to the C
libraries underneath. `terra` binds **GDAL, GEOS and PROJ**, which carry by far the most CVE
history of anything here, and they arrive solely so that three `raster::plot()` calls can
draw a matrix. (`png` → libpng and `jpeg` → libjpeg are the other native bindings, and unlike
this one they are load-bearing.) Since CRAN publishes no advisory database, shrinking the
native surface is the only lever available — see `DECISIONS.md`, "Dependabot watches the
actions, not the R packages".

What has to be replaced is small but not nothing: `main`/`xlab` titling, `axes = F, box = F`,
`add = TRUE` overlay onto a `rasterImage()` background, and the **colour-bar legend** that
`raster::plot()` draws by default (the `legend = F` at line 147 is what makes the undecorated
path simple). `graphics::image()` covers everything except the legend, which would have to be
drawn by hand — reaching for `fields::image.plot` would just trade one dependency for another.

**The release gate will not catch a mistake here.** `tools/compare-harness.R` captures
`ci$zmap`, the z-map *matrix*, which is computed before anything is plotted — so a rendering
change is invisible to it and passes green. Verification has to be visual, plus the
uniform-background trick in `CLAUDE.md`'s test conventions. Pixel orientation is the specific
trap: `raster::plot()` and `image()` do not agree on row order, and getting it wrong flips the
z-map vertically over a base face that is drawn separately by `rasterImage()`.

- [ ] Replace the three `raster::plot()` calls with `graphics::image()`, hand-drawing the
      colour bar for the decorated path.
- [ ] Compare rendered PNGs before and after at 128/256/512px, decorated and not, with and
      without a background image.
- [ ] Drop `raster` from `DESCRIPTION` and the two `importFrom(raster, ...)` lines in
      `NAMESPACE`.

The submission on 2026-07-29 lifted the original hold — the tarball is sent, so nothing here
can move it. **But do not land this while 1.2.1 is in the queue.** If CRAN comes back asking
for a change, the answer should be a minimal 1.2.2 addressing exactly what they asked, not one
that also swapped the plotting backend under three rendered outputs. Wait for the verdict,
then give it the same care item 23 got.

---

### 37. Error paths are largely untested **[verified] [own review]**  ✅ **FIXED**

Found in the pre-submission review, 2026-07-29. The suite is strong on *intent* for the
numeric paths — properties rather than stored numbers, grouping actually checked, both
z-map methods, recovery, parallel equivalence, a golden master. Where it is thin is the
**failure** side: 5 `expect_error()` and 4 `expect_warning()` against **27 `stop()` and 6
`warning()` calls** in `R/`.

That matters more than usual here, because the issue tracker is dominated by confusing
failure modes rather than missing features (item 14 is the general form), and because an
untested error path is how items 6, 23 and 28 each stayed broken for years — a guard that
never runs is indistinguishable from one that works.

Untested, roughly in order of how likely a user is to hit them:

- [ ] **`generateCI()` stimuli/responses length mismatch** (`R/generateCI.R:106`). The single
      most likely user error in the package, and nothing asserts the message.
- [ ] **The four `.Rdata` "did not contain X" guards** in `generateCI()`
      (`R/generateCI.R:128–143`: `s`/`p`, `base_faces`, `stimuli_params`, `img_size`). Only
      the fifth — an unknown base image label — is covered. These are exactly the
      returning-user path the CRAN reinstatement is for.
- [ ] **The five equivalents in `computeCumulativeCICorrelation()`**
      (`R/computeCumulativeCICorrelation.R:53–83`) — none covered.
- [ ] **All four mask-import failures in `generateCI()`** (`R/generateCI.R:407–431`:
      non-greyscale PNG, neither string nor matrix, wrong dimensions, values outside
      {0,1}). `plotZmap()` has one of the four; `generateCI()` has none — and its mask code
      is the copy that has always been live.
- [ ] **`generateStimuli2IFC()` on an unreadable base image file** (`:80`) and on a
      **non-square** base image (`:86`). The related *`img_size` mismatch* is tested
      (`test-fixed-bugs.R:93`); these two are different branches.
- [ ] **`computeCumulativeCICorrelation(targetci = ...)` is never exercised.** The only test
      covers the no-`targetci` path, where the assertion is that the last correlation is 1 —
      true by construction. Supplying a target CI is what the function is actually for.

### 38. Two error messages paste the base image matrix instead of its label **[verified] [own review]**

`stop(paste0('No parameters found for base image: ', base))` at
`R/generateCI.R:177` and `R/computeCumulativeCICorrelation.R:83`. In both, `base` is the
base image **matrix** (assigned a few lines earlier from `base_faces[[baseimage]]`); the
intended variable is `baseimage`, the label. `paste0()` recycles over every pixel, so the
message repeats once per pixel.

**Reachable in `computeCumulativeCICorrelation()`**, measured: `stimuli = integer(0)` gives
an **8,190-character** message at a 32×32 base image, and would give roughly 7 MB at a
realistic 512×512. **Latent in `generateCI()`** — `aggregate()` fails first with "no rows to
aggregate", so the branch cannot currently be reached there; fix both anyway, since what
protects it is an unrelated function's behaviour.

Same shape as the `applyMask()` message that named `img_size` outside its own scope
(recorded under "Also found and fixed while working through the P0 items"): **an error
message referring to a variable that is in scope but is not the one meant**.

✅ **FIXED.** Both sites now name `baseimage`. Only the text of an error changed — the
release gate reports the tree still reproduces v1.2.1 (135 checks identical, 0 unexpected),
and the golden master is green. The regression test in `test-fixed-bugs.R` asserts the
message both matches the intended text and stays under 200 characters; the length bound is
the discriminating half, since the buggy message exceeded it by 7,990 characters.

**The sweep this item asked for is done, and found nothing else.** All 32 `stop()`/
`warning()` sites in `R/` were checked, and the ~14 that interpolate a value all interpolate
a scalar that is the right one: lengths (`generateCI.R:107`), dimensions
(`generateStimuli2IFC.R:86`, `generateCI.R:423`), filenames and labels
(`generateStimuli2IFC.R:80`/`:101`, where `base_face` is a name drawn from
`names(base_face_files)`, not a matrix), a CI name (`autoscale.R:43`), a scaling method
(`generateCI.R:478`), and `baseimage` at the two `did not contain any reference to base image
label` guards. **These two were the only instances of the class.** Do not re-run this sweep;
re-check only when a new message interpolates something.

One unrelated oddity noted in passing, not fixed here: `generateStimuli2IFC.R:68-69` writes
its message to `stderr()` and then calls a bare `stop()`, so the condition carries an empty
message. That belongs to item 14, and its error path is one item 37 should cover.


✅ **FIXED** in `tests/testthat/test-error-paths.R`. Suite **291 → 323**, 0 skipped. Every
box above is covered, plus the `base_face_files`-not-a-list guard found during the item 38
sweep.

**How these were verified, since the usual method does not apply.** This item adds no `R/`
change, so `git stash push -- R/` proves nothing. The real risk for a test of an error path
is that it passes on a *different* error than the guard it names. Each of the 13 was
therefore run directly and its `conditionMessage()` read, confirming it fires its own guard
— and every assertion matches a distinctive fragment of that message rather than merely
expecting "an error".

**Two things worth keeping.** `generateStimuli2IFC.R:68-69` writes its explanation to
`stderr()` and then calls a **bare `stop()`**, so the condition carries an empty message and
no regexp can match it; the test asserts only that it errors, and giving it a real message
belongs to item 14. And the assertions deliberately match short fragments (`"same length"`,
`"did not contain"`, `"other than 0 or 1"`) rather than whole sentences, precisely so item 14
can rewrite the wording without rewriting the suite.

---

### 39. Two of the four `load(rdata)` sites did not guard their arguments **[verified] [own review]**  ✅ **FIXED**

Found while sweeping for item 38, 2026-07-29. `load()` assigns straight into the calling
function's frame, so an object in the `.Rdata` file replaces an argument of the same name.
Of the four functions that read a stimulus set:

| function | state before |
|---|---|
| `generateCI()` | guarded — all arguments kept via `mget(names(formals()))` |
| `generateReferenceDistribution2IFC()` | guarded — each argument restored explicitly |
| `computeInfoVal2IFC()` | **partial** — 3 of 5 arguments (`rdata`, `iter`, `response_seed`) |
| `computeCumulativeCICorrelation()` | **none** |

**There was no live bug.** The 15 names `generateStimuli2IFC()` saves do not intersect the
unguarded arguments, checked rather than assumed. The fix is preventive, and the reason it is
worth having is the direction the one real collision came from: item 32 was created by *adding
a field to the file*, not by adding an argument, so an argument that is safe today stops being
safe with nothing in the function changing.

**The case worth naming** is `computeInfoVal2IFC(target_ci)`. It is read at the very end, to
compute the CI norm, and after a **second** `load()` on the cache-writing path — so a file
carrying that name would have scored somebody else's classification image and returned a
plausible InfoVal rather than an error. Both loads now restore.

Tested by planting the colliding names into a fixture `.Rdata` and asserting the caller's
values survive; verified to fail without the fix (3 failures). Note the tests pass against
unguarded code if the collision is *not* planted, which is why each one plants it explicitly.


## P3 — Features requested by users

Long-standing enhancement requests from the issue tracker, roughly by how often they
recur. None are blocking; all are meaningful to researchers.

- [ ] **Non-2IFC designs** — 4AFC and other variants
      ([#6](https://github.com/rdotsch/rcicr/issues/6)); 2IFC with independently drawn
      (non-inverted) pairs ([#54](https://github.com/rdotsch/rcicr/issues/54)).
- [ ] **Better statistics on CIs** — cluster tests / Random Field Theory
      ([#38](https://github.com/rdotsch/rcicr/issues/38)), cluster-based permutation tests
      per Maris & Oostenveld 2007 ([#46](https://github.com/rdotsch/rcicr/issues/46)),
      ROI analysis ([#47](https://github.com/rdotsch/rcicr/issues/47)), comparing two CIs
      with z-maps ([#71](https://github.com/rdotsch/rcicr/issues/71)).
- [ ] **Correlations between CIs** ([#7](https://github.com/rdotsch/rcicr/issues/7)) and
      equal-weight group CIs ([#27](https://github.com/rdotsch/rcicr/issues/27)).
- [ ] **Batch InfoVal** — compute for many CIs at once
      ([#85](https://github.com/rdotsch/rcicr/issues/85)).
- [ ] **Oblong (non-square) base images** ([#76](https://github.com/rdotsch/rcicr/issues/76)) —
      currently hard-rejected; a frequent stumbling block.
- [ ] **Analysis provenance log** ([#9](https://github.com/rdotsch/rcicr/issues/9)) — write
      package version, seed, and parameters alongside outputs. Cheap, and valuable for
      reproducibility.
- [ ] Refactor batch CI generation, deprecating `batchGenerateCI`
      ([#87](https://github.com/rdotsch/rcicr/issues/87)) — must follow the deprecation
      policy above.

---

### 41. `generateStimuli2IFC()` leaves the user's RNG stream where it landed **[verified] [own review]**

**Open, triage. Size S.** Found 2026-08-07 by the package-wide sweep of CRAN's seven
review points, as the one thing adjacent to their point 7 ("do not change the user's
options, `par` or working directory") that is still true.

`generateStimuli2IFC()` calls `set.seed(seed)` at `R/generateStimuli2IFC.R:82` and never
restores the stream. So in a script that generates stimuli and then does anything else
random, the "anything else" depends on whether stimulus generation ran — and on its
`seed` argument.

**This is not the usual version of that bug.** The seed is a *documented, user-supplied
argument*, and reproducible stimuli are the package's central promise, so the user asked
for the state change. Nothing is silent. Risk of a reviewer raising it is low, and it has
been this way since the CRAN-era 0.x versions.

**Why it was not fixed in 1.2.3:** the fix is to capture `.Random.seed` and restore it
under `on.exit()`, which changes **nothing this package computes** — every stimulus,
classification image, z-map and InfoVal is drawn before the restore, so the release gate
would report no deviation. What it changes is what the *user's next* `runif()` returns.
That is a reproducibility change for any analysis script that draws randomly after
generating stimuli, and it belongs in a version that says so in `NEWS.md` under
"Reproducibility impact", not in one whose entire claim is that nothing computed differs.

Do it with the gate run and a `NEWS.md` entry, once CRAN settles.

### 40. Retire `ChangeLog` as a live file **[verified] [own review]**

**Open, triage. Size S.** Raised 2026-08-07: if `NEWS.md` exists, does `ChangeLog` need to?

**It is not mandatory.** R indexes `NEWS.md`, `NEWS` and `inst/NEWS.Rd` for
`news(package = "rcicr")` and `utils::readNEWS()`; it does not parse `ChangeLog` at all,
so nothing a user or CRAN reads comes from it. It is a GNU convention rather than an R
one. It ships in the tarball and is *not* `.Rbuildignore`d, and `R CMD check --as-cran`
has never NOTEd it — confirming it is on R's list of known top-level files, so it is
permitted, just unused.

**But it cannot simply be deleted.** Its 27 entries run back to 2014 and cover **0.2.2
through 1.0.1**, and `NEWS.md` starts at 1.1.0. Deleting the file would destroy the only
record of the package's first seven years, including the entire CRAN-era history up to the
0.3.4.1 that was archived.

Only the 1.1.0-and-later entries are duplication, and they are already thin — each is a
pointer saying "see `NEWS.md`" plus a three-line summary.

**Proposed:** freeze it. Retitle it as the historical changelog for versions up to 1.0.1,
delete the 1.1.0–1.2.3 pointer entries, and stop adding to it. That leaves each file doing
exactly one job — `ChangeLog` the pre-`NEWS.md` archive, `NEWS.md` everything since — and
removes a per-release step that currently has to be remembered. Three places document the
current convention and would need updating with it: `AGENTS.md`, `CONTRIBUTING.md` and
`DECISIONS.md`.

The alternative, migrating all 27 entries into `NEWS.md` and deleting the file outright,
gives users the full history through `news()` but bloats `NEWS.md` with 2014-era detail in
a format that would have to be converted by hand. Not worth it unless someone asks for it.

**Do this after the CRAN submission settles**, not before: it touches a file inside the
tarball, so doing it now would invalidate the external checks for no benefit CRAN can see.

### 42. The superseded Medium link is the only URL that ever 403s a checker **[verified] [own review]**

**Open, triage. Size S.** Raised 2026-08-07, on the question of whether `cran-comments.md`
needed a paragraph explaining the 403.

`README.md:68` links to the 2016 Medium post. Local `R CMD check --as-cran` reports a 403
on it, from inside the CRAN incoming feasibility NOTE: the site refuses datacenter
networks and resolves normally in a browser.

**It is intermittent, not fixed and not gone.** For 1.2.3 it appeared on *no* external
check — absent from both win-builder runs and from all three R-hub platforms. It was also
absent from CRAN's own pretest of 1.2.1. But win-builder *did* flag it on 1.2.1, so it is
a property of which network the checker sits on, not of the version.

The awkward part is that the README already calls the post superseded: the sentence
carrying the link says the vignette "supersedes it and is the version kept current with
the code". So the package ships a link it tells you not to use, which is the only URL in
it capable of failing a check.

**Why it was not done for 1.2.3:** `README.md` ships inside the tarball, so removing the
link invalidates `rcicr_1.2.3.tar.gz` and forces win-builder and R-hub to be re-run — a
real cost for a note no CRAN-side check has ever raised. `cran-comments.md` is
`.Rbuildignore`d and costs nothing to edit, which is why the *paragraph* went and the
*link* stayed.

**Decide after CRAN settles**, and note the decision is not obvious: dropping the link
removes the only URL that can 403, but the post is the version most existing users will
have bookmarked, so a pointer saying it is superseded has some value. Keeping it costs
nothing as long as nobody explains it in a submission again.

## Beyond v1 — changes that need a major version

Nothing here is scheduled, and nothing here should be started as a side effect of other
work. These are collected so the question does not have to be re-derived each time it comes
up.

### 36. Adopt tidyverse style throughout, as a v2 breaking change

**Not scheduled; no short-term intent.** Logged 2026-07-29 so the shape of the decision is
on record.

The idea is to drop the current mixed convention — camelCase functions, snake_case
arguments — for snake_case throughout, enforced by `styler` and gated by `lintr`:
`generateCI()` → `generate_ci()`, `batchGenerateCI2IFC()` → `batch_generate_ci_2ifc()`, and
so on across all 17 exports.

**The scope is narrower than "adopt tidyverse style" sounds.** Arguments and saved list
fields are *already* snake_case (`img_size`, `n_trials`, `noise_type`), so they conform
today. The break is confined to function names. Everything else — assignment operator,
indentation, `TRUE`/`FALSE`, `seq_len()` — is either already conformant or is ordinary
tidying that needs no version bump at all (`CONTRIBUTING.md` → "Code conventions", and item
13 for the specific cleanups).

**Three things that constrain it:**

- **The `.Rdata` field names must not change, in v2 or ever.** `base_faces`,
  `stimuli_params`, `p`, `img_size`, `seed`, `noise_type` are a *file format*, not code
  style. Renaming them breaks every stimulus set generated since 2014, and the contract is
  append-only. `p` in particular is a poor name and is frozen regardless. A style sweep that
  treats these as ordinary variables is the most likely way to get this badly wrong.
- **Deprecated aliases make it survivable, and are themselves permanent.** Keeping every old
  name as a `lifecycle::deprecate_warn()` shim means old scripts keep running, which is the
  only acceptable path here. But then both names exist indefinitely, doubling the surface to
  document and test — so v2 does not actually retire the old names, it just stops
  recommending them.
- **The names are in the literature.** Published methods sections and shared analysis
  scripts say `generateCI()`. Renaming creates two naming worlds in the record permanently;
  aliases soften that for *running* code but not for *reading* it.

**What is genuinely cheap here:** proving the sweep changed no numbers.
`tools/compare-release-output.R` installs and runs the released code and compares every
output, so a pure rename-and-reformat has to come back with zero differences across the
whole battery. Most projects cannot demonstrate that; this one can. The residual cost is
`git blame`, which no gate can address — hence `DECISIONS.md`'s rule that a `styler` sweep
lands as a commit of its own.

**The case still to be made** is why the benefit exceeds all of that. "Consistent with modern
R" is true and thin. The stronger argument is that the current split invites errors at the
boundary — `img_size` versus `imgSize` is a real thing users get wrong — and that new
contributors arrive expecting snake_case. Neither is obviously worth a major version yet.

---

## Already done (do not redo)

- **Test suite + CI/CD** — testthat suite, GitHub Actions (`R CMD check` on ubuntu ×
  release/devel), covr/Codecov, pre-commit. Closes the substance of issue
  [#41](https://github.com/rdotsch/rcicr/issues/41).
- **Vignette** — `vignettes/getting-started.Rmd`. Closes issue
  [#31](https://github.com/rdotsch/rcicr/issues/31).
- **Documentation pass** — `@examples` added to all 12 exported functions that lacked
  them, plus typo and doc-accuracy fixes. Addresses issue
  [#13](https://github.com/rdotsch/rcicr/issues/13).
- **`generateCI()` core count** — issue
  [#90](https://github.com/rdotsch/rcicr/issues/90) asked for `detectCores() - 1`;
  `generateCI()` already defaults to exactly that. **Verify and close.** (Note the
  *separate*, still-open problem in item #11: `generateReferenceDistribution2IFC()` has no
  `ncores` argument at all.)

### <a name="17"></a>17. InfoVal formula — verified correct, now pinned by a test  ✅

Worth recording so it isn't "re-fixed" by someone reading only the comment paper:

Schmitz et al. (2019) published a
[comment](https://link.springer.com/article/10.3758/s13428-019-01295-1) identifying two
discrepancies in the original infoVal metric (Brinkman et al., 2019): use of the **one
norm** instead of the **Euclidean norm**, and omission of the **k constant** in the
denominator. The subsequent
[erratum](https://link.springer.com/article/10.3758/s13428-020-01367-7) established that
the k-constant point was *not* an error — R's `mad()` already applies `constant = 1.4826`
by default — leaving only the norm issue, of minor impact.

**rcicr already implements the corrected version**: `computeInfoVal2IFC()` uses
`norm(matrix(target_ci[["ci"]]), "f")` (Frobenius = Euclidean on a vectorized CI), and
`mad()` supplies k. Git history confirms both were addressed deliberately —
`01e547e "Euclidean norm and scaling factor k (#97)"` and
`ae8fa9c "Removed k constant, since it is already included in the R mad function."`

- [x] Add a regression test pinning infoVal to a hand-computed expected value, so this
      cannot silently regress. **Done** — `test-regression-baseline.R` pins both the
      formula (recomputed independently) and a literal expected value.
