# rcicr modernization backlog

A prioritized backlog for bringing `rcicr` to a modern, maintainable state **without
breaking the API that researchers depend on**.

**Compiled:** 2026-07-26, against `main` @ `b6ab269` (v1.0.1).

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

If only a few things get done, do these — they are small, high-impact, and all but the
first are contained bug fixes:

| # | Item | Why first | Size |
|---|---|---|---|
| 7 | Base image must already equal `img_size` | Blocks new users at step one; explains open issue #124 | S |
| 6 | `mask` broken on R ≥ 4.2 | Documented feature, 100% broken on current R | S |
| 4 | Tibble inputs fail | Tibbles are the modern default; reported twice, 7 years apart | S |
| 2 | `nscales`/`sigma` missing from `.Rdata` | Silently corrupts a **published statistic** | M |
| 3 | `force_gen_ref_dist` ignored | Silently ignores an explicit user instruction | S |
| 1 | CRAN archived | Highest reach, but a process/decision task, not a code fix | M |
| 9 | Drop 13 unused `Imports` | Trivial deletion; removes install-failure risk, unblocks CRAN | S |
| 16 | Pointless 1.5 GB array per worker | One-line fix; likely resolves most of issue #12 | S |

Items 2, 3, 6 and 7 all share a shape worth noting: **the package fails silently or
misleadingly rather than telling the user what went wrong.** Item 14 (better errors)
is the general form of that problem.

---

## P0 — Correctness and availability

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
      old `r.dotsch@uu.nl` address that was presumably the one that bounced. So the
      remaining work is clearing the `R CMD check` NOTEs below and resubmitting; archived
      packages can be reinstated.
- [ ] Either way, add a `NEWS.md` so users can see what changed between 0.3.4.1 and 1.0.1.

### 2. `nscales` and `sigma` are not saved in the `.Rdata` file → silently wrong InfoVal
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

- [ ] Save `nscales` and `sigma` in the `.Rdata` file.
- [ ] Have `generateReferenceDistribution2IFC()` read and forward them.
- [ ] Fall back to the old default *with a loud warning* when loading a pre-fix `.Rdata`
      that lacks the fields (backward compatibility, but not silent).
- [ ] Regression test asserting the regenerated basis matches the original.

### 3. `computeInfoVal2IFC(force_gen_ref_dist = TRUE)` silently does nothing
**Confirmed by reproduction.** Issue [#113](https://github.com/rdotsch/rcicr/issues/113).

With a pre-existing `reference_norms` in the `.Rdata`, setting `force_gen_ref_dist = TRUE`
leaves the distribution completely untouched (verified: median identical before and
after). The flag short-circuits the *lookup-table* branch but never reaches the
regeneration branch, because `reference_norms` still `exists()` after `load()`.

The user gets no error — they simply believe they forced a recomputation when they did not.

- [ ] Make `force_gen_ref_dist = TRUE` actually bypass the cached `reference_norms`.
- [ ] Regression test.

### 4. `generateCI()` / `generateCI2IFC()` break on tibbles
**Confirmed by reproduction.** Issues [#70](https://github.com/rdotsch/rcicr/issues/70)
and [#123](https://github.com/rdotsch/rcicr/issues/123) — the same root cause, reported
seven years apart, which suggests it keeps catching people.

Passing tibble columns yields `arguments must have same length` from
`aggregate.data.frame()`, because `tbl[, "col"]` stays a 1-column tibble where
`df[, "col"]` drops to a vector. Since `readr`/`dplyr` return tibbles by default, this is
now the *normal* path for a modern user, and the error message gives no hint of the cause.

- [ ] Coerce inputs defensively at the top of `generateCI()`
      (e.g. `stimuli <- unlist(stimuli, use.names = FALSE)`).
- [ ] Do the same in `batchGenerateCI()` / `computeCumulativeCICorrelation()`.
- [ ] Tests covering `data.frame`, `tibble`, and bare-vector input.

### 5. `simulateNoiseIntensities()` is dead on arrival
**Confirmed by reproduction.** Errors 100% of the time with
`object of type 'closure' is not subsettable`: it sizes its progress bar with
`data[, by]`, but neither `data` nor `by` is a parameter of the function (copy-pasted
from `batchGenerateCI()`); `data` resolves to `utils::data`, the function object. It also
ignores its own `img_size` argument, hardcoding `generateNoisePattern(img_size = 512)`.

- [ ] Fix both bugs, or formally deprecate the function if it's no longer intended for use.
- [ ] It is currently pinned by a `\dontrun{}` example and a test documenting the failure —
      update both when fixed.

### 6. The `mask` argument is completely unusable on R ≥ 4.2 **[verified] [own review]**
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

- [ ] Replace both `if (!is.na(mask))` guards with a scalar test
      (e.g. `if (!identical(mask, NA))` or `!is.null(mask)`).
- [ ] Replace the hardcoded `512` with `img_size`.
- [ ] Tests covering a matrix mask and a PNG-file mask at a non-512 `img_size`.
- [ ] Consider defaulting `mask` to `NULL` rather than `NA` (a cleaner sentinel), keeping
      `NA` accepted for backward compatibility.

### 7. `generateStimuli2IFC()` requires the base image to already be exactly `img_size` **[verified] [own review]**
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

- [ ] Validate up front that `dim(img) == c(img_size, img_size)`, with an error naming the
      file, its actual dimensions, and the requested `img_size`.
- [ ] Better: restore automatic resizing with a maintained package
      (`magick::image_resize()`, or `raster`/`terra`), so the original intent works again.
      Gate it behind an argument if silent resizing is considered too implicit.
- [ ] Test covering a base image that doesn't match `img_size`.

### 8. Legacy `sinusoids`/`sinIdx` `.Rdata` support is unreachable
`generateNoiseImage()` validates `length(params) != max(p$patchIdx)` **before** the block
that renames pre-0.3.3 `sinusoids`/`sinIdx` to `patches`/`patchIdx`. So an old-style `p`
hits `max(NULL)` → `-Inf` and always errors with "Stimulus generation aborted", despite
the code clearly intending to support it. This defeats the backward compatibility the
package explicitly promises for old `.Rdata` files.

- [ ] Move the rename block above the validation.
- [ ] Test with a genuine pre-0.3.3-shaped `p`.

---

## P1 — Modernize dependencies and toolchain

### 9. Remove 13 unused `Imports` (and shrink a 27-package dependency surface)
`R CMD check` flags: `DBI`, `abind`, `assertthat`, `deldir`, `ggplot2`, `goftest`,
`gridExtra`, `iterators`, `munsell`, `plyr`, `polyclip`, `sp`, `tensor` — declared but
never used. These look like transitive dependencies of old `spatstat` that were pinned
directly and never cleaned up.

Every one is an install-failure risk for users and a CRAN-resubmission blocker. `sp` in
particular is in the retiring `sp`/`rgdal` ecosystem.

- [ ] Delete all 13 from `Imports`. Low risk: nothing references them.
- [ ] Consider `raster` → `terra` (raster is superseded and depends on `sp`).
- [ ] Move genuinely optional deps (`ggplot2`, `gridExtra` if reintroduced) to `Suggests`.

### 10. Replace deprecated function calls
Verified deprecation warnings on current package versions:

| Call | Status | Used in |
|---|---|---|
| `dplyr::progress_estimated()` | deprecated in dplyr 1.0.0 | `batchGenerateCI`, `batchGenerateCI2IFC`, `computeCumulativeCICorrelation`, `generateReferenceDistribution`, `simulateNoiseIntensities` |
| `purrr::rbernoulli()` | deprecated in purrr 1.0.0 | `generateReferenceDistribution` |
| `citEntry()` / `citHeader()` | superseded by `bibentry()` | `inst/CITATION` |

These still work today but emit warnings on every run and will eventually be removed —
at which point five functions break at once.

- [ ] Swap progress bars for `utils::txtProgressBar()` (already used elsewhere in the
      package, so no new dependency) or `cli`.
- [ ] Replace `purrr::rbernoulli(n, p)` with `stats::rbinom(n, 1, p)` — this may drop the
      `purrr` dependency entirely.
- [ ] Modernize `inst/CITATION` to `bibentry()`.

### 11. Make parallelism robust and user-controllable
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

- [ ] Add `ncores` to `generateReferenceDistribution2IFC()` (default preserving current
      behavior).
- [ ] `on.exit()` cluster cleanup everywhere a cluster is created.
- [ ] Skip the cluster entirely when `ncores == 1`.

### 12. Fill out the test suite
A testthat suite now exists (84 tests) covering the pure functions, light I/O paths, and
an end-to-end smoke test. Gaps worth closing:

- [ ] **Lock in the InfoVal formula.** See the note under item #17 — the implementation is
      currently correct per the published erratum, but *nothing tests it*, so a future
      refactor could silently regress a published metric.
- [ ] Cover the `scaling` methods (`none`/`constant`/`matched`/`independent`) — currently
      only `independent` is meaningfully exercised, and scaling choice is the most
      documented user-facing decision in the package.
- [ ] Cover `mask` handling in `generateCI()` — **currently untestable; see item 6**, where
      masks are shown to be broken outright on R ≥ 4.2.
- [ ] Cover the `zmap = TRUE` paths (`quick` and `t.test`).
- [ ] Cover `participants` (per-participant → group averaging).

---

## P2 — Usability and maintainability

### 13. Modernize the R code itself
- [ ] **[own review]** Version metadata is inconsistent and stale. `generateStimuli2IFC()`
      hardcodes `generator_version <- '0.4.0'` (`R/generateStimuli2IFC.R:168`) into every
      `.Rdata` file, while `generateNoisePattern()` writes the real
      `utils::packageVersion('rcicr')`. So a file written today claims to come from 0.4.0
      — useless for the provenance purpose the field exists for, and it undercuts the
      `pre_0.3.0`-style compatibility checks that key off it. Covers issue
      [#29](https://github.com/rdotsch/rcicr/issues/29). Replace the literal with
      `utils::packageVersion('rcicr')`, but keep *reading* tolerant of old values.
- [ ] **[own review]** The `matlab` package exports its own `sum()` with MATLAB semantics
      (column sums for a matrix, not a single total), which masks `base::sum()` wherever
      `@import matlab` is in effect — six files. Package code is currently safe (its only
      two `sum()` calls are on vectors, where the two agree), but this is a silent trap
      for future edits, and it already bit the test suite. Prefer `@importFrom matlab ...`
      over `@import matlab`, importing only what is actually used.
- [ ] Replace bare `T`/`F` with `TRUE`/`FALSE` (11 occurrences across `autoscale.R`,
      `generateCI.R`, `generateStimuli2IFC.R`, `plotZmap.R`). `T`/`F` are ordinary
      rebindable variables, so this is a genuine (if rare) correctness hazard.
- [ ] Adopt a consistent style (`styler` + `lintr`). Deliberately *not* bundled into the
      current pre-commit config, because it would reformat nearly every file in one sweep
      and destroy `git blame`. Do it as one clearly-labelled commit, or not at all.
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
- [ ] Check `length(stimuli) == length(responses)` before `aggregate()`.

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

### 16. Memory ceiling on large stimulus sets — concrete root cause found **[own review]**
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

- [ ] Replace the pre-allocated array with a per-iteration local variable. Small, safe,
      and likely resolves the bulk of #12 on its own.
- [ ] Re-measure the ceiling afterwards and document it, so users can size jobs up front.
- [ ] Only if still needed: chunking, or regenerating noise from `seed` on demand rather
      than storing every parameter vector.

---

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

### <a name="17"></a>17. InfoVal formula — verified correct, but untested

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

- [ ] Add a regression test pinning infoVal to a hand-computed expected value, so this
      cannot silently regress. Currently nothing in the test suite would catch it.
