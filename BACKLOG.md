# rcicr modernization backlog

A prioritized backlog for bringing `rcicr` to a modern, maintainable state **without
breaking the API that researchers depend on**.

**Compiled:** 2026-07-26, against `main` @ `b6ab269` (v1.0.1).
**Last updated:** 2026-07-26 — P0 items 2–8 fixed (v1.0.1.9000); see `NEWS.md`.

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

**All seven P0 code bugs (items 2–8) are fixed.** What remains, in the order I would
take it:

| # | Item | Why | Size |
|---|---|---|---|
| 9 | Drop 13 unused `Imports` | Trivial deletion; removes install-failure risk and clears an `R CMD check` NOTE, which unblocks CRAN | S |
| 16 | Pointless 1.5 GB array exported to every worker | Roughly a one-line fix; likely resolves most of issue #12 | S |
| 10 | Replace deprecated `progress_estimated()` / `rbernoulli()` / `citEntry()` | Still work today, but warn on every run and will eventually break five functions at once | S |
| 1 | CRAN archived | Highest reach of anything here, but a process/decision task rather than a code fix | M |
| ~~11~~ | ~~Cluster cleanup (`on.exit`), serial fallback~~ | **Done.** `on.exit` cleanup in #130; the `ncores == 1` serial fast path landed via `startBackend()`. Test suite 140s → 4s, and serial/parallel output verified bit-identical | M |
| 12 | Widen test coverage (scaling methods, z-maps, `participants`) | The suite covers the fixed bugs well; these paths are still untested | M |
| ~~18~~ | ~~Codecov step fails for want of a token~~ | **Done** — `fail_ci_if_error: false`; a red `main` now means the package is broken | S |
| ~~19~~ | ~~Close the 8 issues already fixed in `main`~~ | **Done** — 7 closed, #12 commented and left open as only partly fixed. ~22 remaining issues still unswept | S |
| 20 | CRAN resubmission checklist | Verified against a real `--as-cran` run; six concrete fixes then submit | M |
| 21 | Announcement post | Drafted in `notes/`; hold until the CRAN outcome is known | S |
| 22 | Move the Medium walkthrough into a vignette | Removes the last real CRAN NOTE, and makes the tutorial execute at build time so it cannot go stale | M |

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

- [x] Add `ncores` to `generateReferenceDistribution2IFC()` (default preserving current
      behavior). **Done** — this also removed the need for the `detectCores()` mock that
      the test suite had been using to stay under CRAN's 2-core check limit.
- [x] `on.exit()` cluster cleanup everywhere a cluster is created. **Done** — a
      `stopClusterSafely()` helper is registered via `on.exit()` at all three cluster
      sites, so workers are released even when an error interrupts the `foreach` loop.
      Verified: a mid-loop failure now leaks zero connections (issue #50).
- [ ] Skip the cluster entirely when `ncores == 1`.

### 12. Fill out the test suite  ✅ **DONE**
A testthat suite now exists (**177 tests**) covering the pure functions, light I/O paths,
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

**Current status — re-run 2026-07-27 against `main` @ `5428b79`, in a check environment
that now genuinely matches CRAN's.** `texlive` and `tidy` are installed, and the run set
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

Of these, only **"Version contains large components"** is an actual blocker, and it is the
one unticked box below. "New submission" and "Package was archived" are expected for a
reinstatement and are what `cran-comments.md` exists to explain. The `codecov.io` URL that
this NOTE previously also flagged is gone.

#### Must do before submitting

- [ ] **Bump the version to `1.1.0`.** *(left deliberately — it is a release decision,
      not a mechanical fix, and it belongs with the CRAN go/no-go in item 1.)* `1.0.1.9000` trips "Version contains large
      components" — the `.9000` development suffix is not acceptable in a submission.
      Date the `NEWS.md` heading at the same time.
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

- [ ] Check on platforms this machine cannot provide: `devtools::check_win_devel()` and
      `rhub::rhub_check()` (macOS, Windows, R-devel). CRAN tests those and we do not.
- [x] **Write `cran-comments.md`. Done** — at the repo root, `.Rbuildignore`d. It states
      plainly that the 2021-06-08 archival was for an undeliverable maintainer address
      rather than any code or policy problem, that the address now works, and explains both
      remaining NOTEs. Two `<!-- TODO -->` markers remain for the win-builder and R-hub
      results, which cannot be filled in until those actually run.
- [ ] Submit at <https://cran.r-project.org/submit.html> (or `devtools::release()`).
      **Ron has to do this himself** — CRAN emails the maintainer address for
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

- [ ] Port the walkthrough into a vignette (alongside `getting-started.Rmd`), with its
      code chunks actually running so they stay honest.
- [ ] **Keep the Medium post published**, updated with a pointer to the package docs. It
      has nine years of inbound links and citations; deleting it would break them. This
      is a *move of the canonical copy*, not a takedown.
- [ ] Once the content lives here, the remaining `medium.com` references in `README.md`
      and `vignettes/getting-started.Rmd` become "further reading" rather than the primary
      source — and if they are dropped entirely at that point, the CRAN URL NOTE goes with
      them.

Worth doing before the CRAN submission if there is appetite, since it removes the one
non-trivial NOTE. Not a blocker: the NOTE is explainable in `cran-comments.md`.

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
