# rcicr (development version)

## Behaviour changes

- **A decorated z-map on a device too small for it now says so, and can be made to fit.**
  `plotZmap(decoration = TRUE)` needs room for margins, labels and the colour scale — about
  160px at the default text size on Linux and macOS, and about 200px on Windows, whose
  graphics device is 96 ppi where theirs is 72 — and below that base R stopped with `figure
  margins too large`, naming neither `rcicr` nor a way forward. It now stops with a message giving the
  size, the minimum, and the three ways out. `plotZmap()` gains **`pointsize`** and
  `generateCI()` gains **`zmappointsize`**: the decoration is measured in lines of text, so a
  smaller size fits it onto a smaller image. This matters most through `generateCI()`, which
  sizes the z-map to `img_size`, so a stimulus set below about 160px could not produce a
  decorated z-map at all — a `128px` set now can, with `zmappointsize = 6`.

  Both arguments default to the graphics device's own `12`, so every existing call renders
  exactly as before; a 512px decorated z-map is byte-identical.

- **`plotZmap()` no longer depends on `raster`, and its `...` arguments now go to
  `graphics::image()` instead of the raster package's plot method.** `col` behaves the same
  way in both (and now actually works — see the bug fixes below); a call passing an argument
  specific to that method will now be rejected as unused. Dropping the dependency also removes
  `terra`, `sp` and `Rcpp`, and with them the GDAL/GEOS/PROJ system libraries that every
  Linux CI job had to install before it could start.

  **The z-map itself renders identically** — the undecorated figures `generateCI()` writes
  are pixel-for-pixel the same as before, within colour quantisation, and a golden reference
  rendered by the old `raster` code is committed as a test fixture to keep it that way. The
  palettes are unchanged, including the quirk that a z-map drawn over a background image uses
  the default palette rather than the viridis one.

- **A `decoration = TRUE` z-map is laid out slightly differently.** The colour bar is now
  drawn by hand rather than by `raster`, and the map is a few pixels wider at 512px. If you
  regenerate a decorated z-map figure, it will not be pixel-identical to one saved with 1.2.3
  — the same is already true of regenerating it on a different operating system, and for the
  same reason: what a graphics device paints is not part of what this package computes. No
  z-score, classification image, scaling result or informational value changes, and
  `generateCI()`'s own z-map figures are undecorated and unaffected. See
  `?plotZmap`, "Reproducibility across platforms".

- **`generateStimuli2IFC()` now checks `base_face_files` before it generates anything,
  and names the entry it cannot use.** Four inputs used to get past the old check and
  fail from inside a parallel worker with `attempt to select less than one element in
  get1index`, which names neither the argument nor the file: a list with no names, a
  list with some names missing, an empty list, and an element that is not a single file
  name. They now stop immediately with a message saying which entry is wrong and why.

  One of them could previously appear to work. **A `base_face_files` with two entries
  under the same name silently dropped all but the first** — the loop looked each name
  up by string, so `list(face = 'a.png', face = 'b.png')` produced one set of stimuli,
  from `a.png`, and nothing at all from `b.png`. That is now an error naming the
  duplicated name. If you have a script that relies on it, the stimuli it produced were
  never what the call asked for; give each base image its own name.

- **The PNG-or-JPEG test now looks at the file extension rather than anywhere in the
  path.** It was `grepl('png|PNG', filename)` against the whole path, so a JPEG stored
  under a directory called `png` was handed to `png::readPNG()` and died with `file is
  not in PNG format`, blaming the file for a choice the package had made. Files are now
  recognised by a `.png`, `.jpg` or `.jpeg` extension, case-insensitively. A base image
  whose *extension* does not say what it is is rejected up front, by name, instead of
  reaching a reader that cannot parse it.

## Bug fixes

- **`plotZmap(col = ...)` works.** Supplying a palette is how `?plotZmap` has always told
  you to change the colours, and doing it stopped the call with `formal argument "col"
  matched by multiple actual arguments` before anything was drawn — the function passed its
  own `col` alongside yours. This affected every released version; the argument is now taken
  as an override, and is used for the colour bar as well as the map.

- **A base image with no contrast no longer becomes an all-`NaN` base image.** Under
  `maximize_baseimage_contrast = TRUE`, the default, `generateStimuli2IFC()` rescales
  with `(img - min(img)) / (max(img) - min(img))` — which is 0/0 when every pixel is the
  same value. The resulting `NaN` base face was written into the `.Rdata` with no error
  and no warning, every classification image computed from that stimulus set inherited
  it, and the stimuli themselves came out uniformly black, `png::writePNG()` clamping
  `NaN` to zero. It now stops with an error naming the file.

  A photograph is never uniform, so this bit synthetic and accidentally-blank base
  images — but it failed silently, and the symptom appeared a long way from the cause.
  The error fires only under `maximize_baseimage_contrast = TRUE`: a flat base image is
  perfectly usable with the rescale switched off, and the message says so.

- **The `base_face_files` type check raises an error you can actually read.** It wrote
  its explanation to `stderr()` and then called `stop()` with no arguments, so the
  condition it raised carried an empty message: `conditionMessage()` returned `""`, and
  anything wrapping the call — a Shiny app, a batch script, knitr — caught an error it
  could not report. Under `sink()` or `capture.output(type = "message")` the explanation
  could be separated from the failure entirely. The text now travels with the condition.

- **A base image that cannot be read is reported with the file it came from.** The
  reader's own complaint (`file is not in PNG format`, and the like) is kept and
  prefixed with the entry's name and path, and a file that does not exist is now
  reported as missing rather than as unopenable.

- **Base images are validated before the noise basis is built and before the output
  directory is created.** A mistyped path used to cost the full noise-pattern
  computation — appreciable at the default 512px — and left an empty stimulus directory
  behind before reporting the problem.

  No numeric output changes. Every check here either replaces an error with a clearer
  error, or rejects input that could not have produced correct stimuli; a call that
  succeeds today produces exactly what it did before.

## Documentation

- **There is now a documentation website: <https://rdotsch.github.io/rcicr/>.** The
  function reference, both vignettes and this changelog are readable without installing
  the package first. It is generated from the same sources — `README.md` is the home page,
  so "How it works" and "Anatomy of the `.Rdata` file" are not maintained twice — and
  rebuilt by GitHub Actions on every push, so it cannot drift from the code.

# rcicr 1.2.3 (2026-08-07)

**Documentation only. Nothing this package computes has changed** — no function,
argument, return value or number differs from 1.2.2, and no analysis script needs
revisiting.

The release exists because the package-level help page, `?rcicr`, had gone stale enough
to contradict the release before it. It was a hand-maintained `.Rd` file that roxygen
never touched, so the sweeps behind 1.2.2 — which all worked from `R/` — went straight
past it. Every code snippet on it was wrong:

* Two of the three example calls, `generateStimuli2IFC(base_face_files, n_trials = 770)`
  and `generateCI2IFC(stimuli, responses, baseimage, rdata)`, would now **error**: they
  omit the write paths that 1.2.2 made required.
* The third, `autoscale(cis, saveasjpegs = TRUE)`, named an argument that has never
  existed under that spelling. It is `save_as_pngs`.
* The page twice promised output "saved as jpegs to a folder called stimuli in your
  current working directory" — the writing-by-default behaviour that 1.2.2 removed, and
  in a format the package does not write.
* Its `\examples` section was a single commented-out line, `#simple examples will be
  added soon.`, left over from 2016.
* It carried a hand-typed `Version: 0.4.0` and `Date: 2017-07-25`, five releases and
  nine years out of date.

**The page is now generated by roxygen from `R/rcicr-package.R`**, so its title,
description, author and URLs come from `DESCRIPTION` and cannot drift from it again. The
version and date table is gone rather than corrected — the way to keep a fact current is
to stop writing it down twice. What remains is a short pointer to
`vignette("reverse-correlation-walkthrough")` and to the three functions a new user
starts with; the walkthrough itself lives in the vignette and `README.md`, which are
tested on every build.

The page also now lists all three key references — Dotsch & Todorov (2012), Brinkman,
Todorov & Dotsch (2017), and Dotsch, Wigboldus, Langner & Van Knippenberg (2008) — each
with a DOI. The 2008 paper had been on the old page and was the only one not carried
anywhere else.

`?rcicr` and `package?rcicr` both still work.

# rcicr 1.2.2 (2026-08-07)

This release exists to answer the changes CRAN asked for when reviewing the previous
submission. Nothing it changes affects a number this package computes: classification
images, scaling, z-maps and informational value are identical to 1.2.1, and the release
gate confirms that against both 1.2.1 and 1.0.1.

## Breaking changes

- **Functions that write files now require you to say where.** `stimulus_path`
  (`generateStimuli2IFC()`), `targetpath` (`generateCI()`, `generateCI2IFC()`,
  `batchGenerateCI()`, `batchGenerateCI2IFC()`, `autoscale()`, `plotZmap()`) and
  `zmaptargetpath` (`generateCI()`) have lost their defaults. They used to be `./stimuli`,
  `./cis` and `./zmaps`, which meant a default call created directories in whatever your
  working directory happened to be — writing to your filespace without being asked, which
  CRAN policy does not permit.

  **What to change in your scripts.** If you relied on the old defaults, name them:

  ``` r
  # before
  generateCI(stimuli, responses, "face", rdata)

  # after
  generateCI(stimuli, responses, "face", rdata, targetpath = "./cis")
  ```

  You will not silently get files somewhere new — a call that would have written to a
  default path now stops with an error naming the argument to supply. If you do not want
  files at all, `save_as_png = FALSE` (or `save_as_pngs = FALSE` for `autoscale()`) needs
  no path.

## Bug fixes

- **`batchGenerateCI()` no longer produces a spurious CI for rows with no group.** Rows
  whose `by` column was `NA` were kept and collapsed into an extra group named after `NA`,
  so a data frame with any missing grouping value returned one more classification image
  than it had groups — computed from whatever rows happened to be missing that value.
  `batchGenerateCI2IFC()` has always dropped those rows; the two now agree.

- **`generateCI(mask = )` accepts a logical matrix.** The matrix branch tested
  `typeof(mask) == 'double'`, so a mask built the obvious way — as `TRUE`/`FALSE` rather
  than `1`/`0` — fell through to `The mask argument is neither a string nor a matrix!`,
  despite the documentation describing exactly that form. It is now tested with
  `is.matrix()`.

- **`generateCI()` and `computeCumulativeCICorrelation()` no longer print the entire base
  image when they cannot find stimulus parameters.** The "No parameters found for base image"
  error named the base image *matrix* where its label was meant. Because `paste0()` is
  vectorized, this did not paste one matrix into one message — it built one complete message
  per pixel, so the error came back as 1,024 concatenated copies at a 32x32 base image (8,190
  characters) and roughly 7 MB at the 512x512 size researchers actually use, with the reason
  for the failure buried inside it. The message now reads, in full, `No parameters found for
  base image: <label>`.

  Only the text of an error changed. No function's return value, arguments or numeric output
  are affected, and the condition that triggers the error is unchanged — if your analysis
  script runs today, it behaves identically.

- **`generateReferenceDistribution2IFC()` no longer leaves a stray `stimuli` directory
  behind.** It re-derives the noise basis by calling `generateStimuli2IFC()` with both save
  options off, purely to work in memory — but the directory was created before either
  option was consulted, so every call to it, and to `computeInfoVal2IFC()` when no
  reference distribution was cached, created an empty `./stimuli` wherever you happened to
  be working.

- **`plotZmap()` restores the graphics parameters it changes.** The undecorated branch set
  `par(mar = ...)` and left it set. It also now closes its PNG device through `on.exit()`,
  so a failure part-way through plotting can no longer leak the device or leave a
  half-written file.

## Documentation

- The `DESCRIPTION` description no longer opens with the redundant "Functions to", and
  cites the two method references: Dotsch and Todorov (2012)
  <doi:10.1177/1948550611430272> and Brinkman, Todorov and Dotsch (2017)
  <doi:10.1080/10463283.2017.1381469>.

- **Every example runs.** The `\donttest{}` wrappers are gone from all eight examples that
  carried them, `simulateNoiseIntensities()`'s `\dontrun{}` is gone, and
  `generateNoiseImage()`'s example is real code rather than three commented-out lines that
  would not have worked (`p` was never defined, and `params` was the wrong length for the
  pattern). The whole example set now runs in about nine seconds.

- `simulateNoiseIntensities()`'s note claiming the function always errors is removed. It
  described two bugs that were fixed in 1.1.0; the note was left behind.

## Internal

- Bare `T` and `F` are replaced by `TRUE` and `FALSE` throughout `R/`. Two of these were
  public API defaults visible in the documentation (`generateCI(zmap =, zmapdecoration =)`
  and `plotZmap(decoration =)`); the values are unchanged.

- The guard that keeps a function's arguments across `load()` is now a shared helper,
  `captureArgs()`. It skips required arguments that were not supplied — necessary once
  paths became required, since `mget()` forces the promise and a wrapper forwarding its own
  missing argument would abort there. Defaulted arguments are still captured: `missing()`
  reports those missing too, and their default is exactly as vulnerable to being replaced
  by a field in the `.Rdata` file as a value passed explicitly.

- **The failure paths are tested.** The suite had 9 assertions covering 33 `stop()` and
  `warning()` calls, so most of the package's error messages had never been run. They now
  are: the stimuli/responses length mismatch, every "this `.Rdata` file did not contain X"
  guard in `generateCI()` and `computeCumulativeCICorrelation()`, all four mask-import
  failures, and base images that are unreadable or not square. No behaviour changed — this
  is coverage of messages that were already there. It matters because an unexercised guard
  is indistinguishable from one that works, which is how three separate bugs in this package
  stayed live for years, the most recent being the one fixed just above.

- **Every function that reads a stimulus set now keeps its arguments across the
  `load()`.** `load()` assigns straight into the calling function's frame, so an object
  stored in an `.Rdata` file silently replaces an argument of the same name. `generateCI()`
  and `generateReferenceDistribution2IFC()` already guarded against this; `computeInfoVal2IFC()`
  guarded three of its five arguments, and `computeCumulativeCICorrelation()` none.

  **No file this package has ever written triggers the problem**, so no result changes and
  no analysis needs revisiting — the guard is preventive. It is worth having because the one
  collision that did occur (the z-map `sigma`, fixed in 1.2.0) was created by *adding a field
  to the file*, not by adding an argument, so an argument that is safe today stops being safe
  without anything in the function changing. The case now closed in `computeInfoVal2IFC()` is
  the one that would have mattered most: `target_ci` is read at the very end to compute the
  CI norm, and after a second `load()`, so a file carrying that name would have scored a
  different classification image and returned a plausible number rather than an error.

# rcicr 1.2.1 (2026-07-28)

**No user-facing changes.** Nothing this package computes differs from 1.2.0 — no function,
argument, return value or number has changed, and no analysis script needs revisiting.

The release exists because the 1.2.0 source tree does not pass `R CMD check` on macOS. The
fault was in the package's own test suite, not in the package: a test asserted properties of
a rendered PNG that belong to the graphics device rather than to what was drawn, and those
properties differ between macOS and Linux. No released function was ever affected. A package
still has to pass its own checks on the platforms CRAN builds for, which is what this
release restores.

## Documentation

- **`?plotZmap` and the README now describe what is and is not reproducible across operating
  systems.** Classification images, scaling, informational value and the z-scores themselves
  are ordinary R arithmetic and do not depend on your platform — as of this release that is
  verified on Linux, macOS and Windows on every change, rather than assumed. The PNG written
  by `plotZmap()` is the one exception: it is drawn through a graphics device, and devices
  differ by platform in colour management and in whether they write an alpha channel, so the
  same z-map produces visibly identical figures whose files are not byte-identical.

  The practical advice, now stated in both places: **compare numbers, not rendered figures**,
  when checking that an analysis reproduces. A z-map image that differs pixel-for-pixel on a
  colleague's machine is not a different result. Every other PNG the package writes —
  stimuli, classification images, autoscaled classification images — is written directly
  from the pixel array and carries no such dependence.

## Internal

- `R CMD check` now runs on macOS and Windows as well as Linux, on every change. It
  previously varied only the R version against a single platform, which is how the macOS
  failure above went unnoticed.
- The test suite pins z-map values as well as classification images, scaling and
  informational value, so cross-platform agreement of the numbers is checked rather than
  assumed.

# rcicr 1.2.0 (2026-07-28)

*Upgrading from the CRAN version?* The last release on CRAN was 0.3.4.1, before the package
was archived in 2021. 1.0.1 and 1.1.0 were GitHub-only releases made in the meantime, so the
1.1.0 section below applies to you too — it is where the bulk of the bug fixes are.

## Reproducibility impact

- **`generateCI(zmap = TRUE)` blurred z-maps with the wrong sigma, for stimulus sets
  generated with 1.1.0 only.** `generateCI()` reads the stimulus set with `load()`, which
  assigns into the function's own frame, and 1.1.0 began storing the *noise* `sigma` there —
  the same name as the z-map blur argument. The saved value replaced the argument, so z-maps
  were blurred with 25 rather than the documented 3, and passing `sigma` did nothing.

  **In practice this affects nobody.** 1.1.0 was a GitHub tag that stood for about a day and
  was never on CRAN, so almost no stimulus set carries the field that triggers it. Files from
  1.0.1 and earlier have no `sigma` at all and were never affected. Only z-maps are involved
  — classification images, scaling, InfoVal and every saved number are untouched. It is
  recorded here because it did change a number, and because if you are the one person who
  generated stimuli that day, regenerating the z-map is a one-line rerun.

  Every argument is now kept across the `load()`, so a field added to the `.Rdata` later
  cannot quietly capture another one. Found by `tools/compare-release-output.R`, the release
  gate introduced in this version — the first bug it caught.

## Behaviour change

- **`plotZmap(mask = ...)` now actually masks the z-map.** The argument has been documented
  since 2016 — "if a cell evaluates to TRUE, the corresponding zmap pixel will be masked" —
  and until now it did nothing at all: the mask was read from its PNG or matrix, checked
  against the z-map's dimensions, validated as binary, converted to boolean, and then
  discarded before plotting. A correct mask produced an unmasked z-map, with no error and no
  warning. Masked cells are now dropped from the z-map exactly as sub-threshold cells are.

  **Who is affected.** Only direct calls to `plotZmap(mask = ...)`. `generateCI()` does not
  pass `mask` to `plotZmap()` — it masks the classification image itself, via a separate and
  working code path — so z-maps produced through the normal pipeline are unchanged, and no
  stored numbers change anywhere. If you have been passing a mask and your z-maps looked
  unmasked, that is why; they will now come out masked, and the earlier images were wrong
  about which regions carry signal.

  A second bug is fixed alongside: the conversion to boolean set every cell to `FALSE`
  whatever you passed, so even once applied the mask would have masked nothing.

- **The `mask` convention is documented correctly for the first time**, in both
  `?plotZmap` and `?generateCI`. Both said a *matrix* masks where the value is `1`/`TRUE`
  while a *PNG* masks where it is black (`0`) — two opposite conventions in one sentence.
  `generateCI()`'s implementation has always masked where the value is `0`, for a matrix and
  a PNG alike, so the matrix half of the documentation was simply wrong. **The code is
  unchanged and the documentation now matches it**, since existing masks were built against
  the behaviour, not the prose. `plotZmap()` follows the same single convention, so one mask
  can be passed to both functions — which is now asserted by a test rather than assumed.

  If you built a mask by reading `?generateCI` rather than by looking at your output, check
  it: `0`/black/`FALSE` is the region that gets masked away.

## New features

- `generateReferenceDistribution2IFC()` and `computeInfoVal2IFC()` gained a
  `response_seed` argument, so the null distribution InfoVal is scored against can be
  varied deliberately. Until now there was no way to draw a second, independent null from
  the same stimuli — which meant you could not check how much Monte Carlo error your choice
  of `iter` was leaving in your InfoVal. `response_seed` seeds the simulated responses only;
  the stimuli, and so the noise basis the null is built on, are untouched.

  **Existing calls are unaffected.** The default (`NULL`) issues no `set.seed()` call at
  all, so the reference distribution is byte-identical to what earlier versions produced.
  Verified against norms generated before the change, not merely assumed.

  In `computeInfoVal2IFC()`, passing `response_seed` forces the reference distribution to be
  regenerated even when the `.Rdata` file already holds one, and the result is deliberately
  *not* written back — a one-off check of the Monte Carlo error cannot silently become the
  number every later analysis of that stimulus set reports.

- `generateReferenceDistribution2IFC()` gained `save_rdata` (default `TRUE`, i.e. unchanged)
  and now returns the reference distribution invisibly instead of returning nothing, so the
  norms are reachable when you ask it not to write them to the `.Rdata` file.

- The `.Rdata` file gained a `reference_norms_seed` field recording the `response_seed` the
  stored `reference_norms` were generated with (`NULL` for the default). Purely additive;
  files written by earlier versions simply lack it. A stimulus set carrying a deliberately
  varied null is no longer indistinguishable from one carrying the default.

## Bug fixes

- `autoscale()` works on masked classification images. `generateCI(mask = ...)` sets masked
  pixels to `NA` by design, and `autoscale()` took a bare `range()` over them, so the scaling
  constant became `NA` and the call died with `missing value where TRUE/FALSE needed`. The
  scaling constant is now computed from the unmasked pixels, exactly as `generateCI()`'s own
  scaling has always done, and masked pixels stay masked in the result. A CI that is
  *entirely* `NA` now raises an error naming the CI instead of failing the same opaque way.

- `generateCI()` accepts a pre-0.3.0 `.Rdata` file when computing a CI from a **single**
  trial. rcicr 0.3.0 stopped drawing four random contrasts per trial that no patch index
  ever referred to (4096 → 4092), and `generateCI()` has truncated older files ever since —
  but the single-trial branch tested for a length of 4092 and then truncated to 4092, a
  no-op that could never fire on the 4096-parameter input it existed for. Such a call failed
  with `Stimulus generation aborted: number of parameters doesn't equal number of patches!`.
  The multi-trial path was always correct and is unchanged.

- `computeInfoVal2IFC()` and `generateReferenceDistribution2IFC()` work on `.Rdata` files
  written before `noise_type` was saved (#94). Such a file failed outright with
  `object 'noise_type' not found`, and the workaround on record was to load the file and
  assign the variable by hand. It now falls back to `sinusoid` with a loud warning, matching
  how `nscales` is handled — a warning rather than a silent default, because guessing wrong
  means the null is built on a different *kind* of noise than participants saw, and the
  resulting InfoVal would be wrong. Files written by 1.1.0 or later already store the field
  and are unaffected.

- `generateStimuli2IFC(return_as_dataframe = TRUE)` shows its progress bar (#82). The
  `return` handing back each trial's noise exits the entire loop body, so it jumped past the
  progress-bar update and the bar sat at zero for the whole run — on the slowest path there
  is, since `generateReferenceDistribution2IFC()` takes it for every InfoVal.

- The `.Rdata` file written by `generateStimuli2IFC()` now records the rcicr version that
  actually wrote it (#29). `generator_version` was a hardcoded `'0.4.0'` string from 2016
  onwards, so every file produced by 0.4.0 through 1.1.0 claims to come from 0.4.0 —
  useless for the provenance the field exists for, and it disagreed with
  `p$generator_version`, which held the real version all along.

  No result changes: nothing in the package has ever read this field. If your own code
  does, note two things. Existing files cannot be trusted to say what wrote them, so treat
  a stored `'0.4.0'` as "unknown, somewhere between 0.4.0 and 1.1.0" rather than as a
  version. And the field is now a `package_version` object rather than a character string,
  so compare with `utils::packageVersion()` or `numeric_version()`, never as text —
  `'0.10.0' < '0.4.0'` is `TRUE` when compared as strings.

## Documentation

- `README.md` now describes the package's architecture and, more usefully, the **anatomy of
  the `.Rdata` file** field by field. That file is the only link between stimulus generation
  and analysis — nothing about a stimulus set is recoverable without it — and until now its
  contents were documented nowhere a user would look.

- A `CONTRIBUTING.md` sets out how to work on the package, leading with the constraint that
  makes it unusual: researchers re-run old analysis scripts years later and publish the
  results, so numeric output must not change silently.

- `?generateStimuli2IFC` documents a restriction on `return_as_dataframe = TRUE`: the frame
  holds one noise image per trial, so it is meaningful only under the default
  `use_same_parameters = TRUE`. With `use_same_parameters = FALSE` and more than one base
  image, only the first base image's noise comes back — the frame's shape cannot represent
  trial × base image. Behaviour is unchanged, and the files written to disk were never
  affected; the restriction simply was not stated.

- `?generateReferenceDistribution2IFC` now documents as a **guarantee** what was previously
  only true by accident: with the default `response_seed`, the reference distribution — and
  therefore InfoVal — is reproducible from the stimulus `.Rdata` file alone, independent of
  the calling session's random number state and of `ncores`. This held before, but nobody
  had chosen it: it is a consequence of `generateStimuli2IFC()`'s internal `set.seed()`
  landing before the simulation draws. That call now carries a comment saying what depends
  on it, so it is not moved casually.

# rcicr 1.1.0 (2026-07-27)

> First release since 1.0.1, and the version submitted to CRAN to reinstate the
> package after its 2021 archival. The minor bump rather than 1.0.2 is deliberate:
> some of the changes below alter behaviour rather than only repairing it. The
> public API is unchanged — no function, argument or argument meaning was removed
> or redefined — so a 2.0.0 is not warranted. Anyone re-running an old analysis
> script should read the "Reproducibility impact" section below.

## Bug fixes

- `generateCI(zmap = TRUE, zmapdecoration = FALSE)` works at all. The undecorated
  branch of `plotZmap()` tested its background image with `if (bgimage != '')`, a
  condition of length `img_size^2`, which R >= 4.2 treats as an error rather than
  silently taking the first element. `generateCI()` always supplies a background
  image, so this path could never run. Same root cause as the `mask` bug above; the
  decorated branch already used `identical()`.
- `plotZmap(decoration = FALSE)` works on small images. `plot.new()` was called
  before the margins were reset to zero, and it rejects a device too small to hold
  the default margins, so any z-map below roughly 100 pixels failed with
  `figure margins too large`. Rendered output at usual sizes is unchanged.
- `generateCI(zmaptargetpath = ...)` is honoured. The argument was documented and
  accepted but never forwarded to `plotZmap()`, so z-maps were always written to
  `./zmaps` relative to the working directory regardless of what was requested.
- `generateStimuli2IFC()` now validates that the base image matches `img_size`, with an
  error naming the file and both sizes, instead of failing inside a parallel worker with
  `non-conformable arrays` (#124).
- The `mask` argument of `generateCI()` works again. It was unusable on R >= 4.2, was
  restricted to 512px stimuli by a hardcoded size check, and rejected greyscale-as-RGB
  PNG masks even when it had successfully converted them.
- `generateCI()`, `batchGenerateCI()` and `computeCumulativeCICorrelation()` accept
  tibble columns, not only `data.frame` columns (#70, #123). `generateCI()` also checks
  that `stimuli` and `responses` are the same length up front.
- `generateStimuli2IFC()` saves `nscales` and `sigma` to the `.Rdata` file, and
  `generateReferenceDistribution2IFC()` uses them (#81). See below.
- `computeInfoVal2IFC(force_gen_ref_dist = TRUE)` regenerates the reference distribution
  instead of silently reusing the cached one (#113).
- `generateNoiseImage()` accepts pre-0.3.3 `.Rdata` noise patterns using the
  `sinusoids`/`sinIdx` layout. This backward-compatibility path had never worked.
- `simulateNoiseIntensities()` runs at all, and honours its `img_size` argument.
- `generateReferenceDistribution2IFC()` gained an `ncores` argument instead of always
  using `detectCores() - 1`.
- `computeInfoVal2IFC()` no longer creates a stray file called `data` in your working
  directory. One `write()` call omitted its destination, and base R's `write()` defaults
  to `file = "data"` — so the "reference distribution has been saved" message was written
  to that file instead of the console.
- `generateCI()` accepts a greyscale-encoded RGB PNG as a `mask`. The conversion worked,
  but the error for genuinely non-greyscale images was raised unconditionally afterwards,
  so every such mask was rejected anyway.
- `generateReferenceDistribution2IFC()` no longer writes its own arguments into the
  `.Rdata` file it updates. `load()` assigns straight into the calling function's frame,
  so a stored `rdata` (or, newly, `ncores`) object overwrote the argument of the same name
  on the next call — meaning a second call ignored the `ncores` you passed and wrote back
  to the path recorded during the first call. `computeInfoVal2IFC()` additionally guards
  against this when reading `.Rdata` files written by older versions, which still contain
  the stale path.

## Performance and dependencies

- `generateStimuli2IFC()` no longer preallocates the full stimulus array before starting
  its parallel workers. At the default 512px / 770 trials that array was ~1.5 GB, and
  `foreach` copied it to every worker, which each wrote one slice into and discarded.
  Addresses issue #12.
- Parallel clusters are now stopped via `on.exit()`, so workers are released even when an
  error interrupts the loop. Fixes the "closing unused connections" warnings (issue #50).
- `generateNoiseImage()` is about **6x faster** — 1.66s to 0.28s per call at the default
  512px with `nscales = 5` — and allocates about 30% less memory. The per-pixel average
  across patch layers is now computed with `rowMeans(..., dims = 2)` instead of
  `apply(..., 1:2, mean)`; that step alone is ~31x faster, and what remains is building
  the weighted patch array, which is unavoidable. Because this function is called for
  every trial during stimulus generation and again for every CI and z-map, the saving
  compounds. Thanks to [@hvalev](https://github.com/hvalev), who diagnosed this and
  benchmarked it in #122. See "Reproducibility impact" below — the result is not
  *bit*-identical to the old one.
- **`ncores = 1` no longer builds a cluster.** The parallel loops in
  `generateStimuli2IFC()` and `generateCI()` previously started a one-worker cluster even
  when asked for a single core — a second R process, loading the package, receiving each
  iteration and sending it back. They now run in the current process. The test suite went
  from 140s to 4s, and `R CMD check` from `[8s/126s]` on tests (eight seconds of CPU
  against 126 elapsed, over 22 cluster spawns) to a fraction of that. Results are
  unaffected and this is verified rather than asserted: `ncores = 1` and `ncores = 2`
  produce bit-identical stimulus parameters, noise bases and classification images, now
  pinned by `tests/testthat/test-parallel-equivalence.R`. Neither loop draws random
  numbers — parameters are drawn under `set.seed()` before the loop — so there is no
  worker RNG stream to diverge. Addresses the remainder of issue #50's neighbourhood.
- `Imports` shrank from 27 packages to 15; none of the removed ones were used.
- Deprecated calls replaced: `dplyr::progress_estimated()`, `purrr::rbernoulli()`, and
  `citEntry()`/`personList()` in `inst/CITATION`. The `rbernoulli()` replacement was
  chosen to preserve the random number stream exactly — see below.

## Reproducibility impact — read this if you have published or in-progress results

This release fixes several long-standing bugs. Most of them turn a **crash into
working behaviour**, and therefore cannot change any result you already have —
if the old code errored, there was no number to change.

Only two fixes can change a value you may already have reported. Both are listed
first, with the exact conditions under which they apply.

A golden-master test (`tests/testthat/test-regression-baseline.R`) pins the
numeric output of the default pipeline — noise basis, classification image,
scaling, and infoVal — to the values produced *before* these fixes. It passes,
which is the evidence for the claim that default-configuration results are
unchanged.

### Results CAN differ

**1. infoVal, if you used a non-default `nscales` (issue #81)**

`generateStimuli2IFC()` did not save `nscales` or `sigma` into the `.Rdata` file.
`generateReferenceDistribution2IFC()` re-generates the stimulus set from that
file to build the infoVal null distribution, so with `nscales` missing it
silently fell back to the default `nscales = 5` — building the reference
distribution on a different noise basis than the stimuli participants actually
saw.

- **Affected:** anyone who called `generateStimuli2IFC(..., nscales = <not 5>)`
  and then reported an infoVal (or `sigma` with `noise_type = "gabor"`).
- **Not affected:** the default `nscales = 5`, where the fallback happened to
  coincide with the real value. This is the common case.
- **After the fix:** infoVal for affected stimulus sets will change, because it is
  now computed against the correct null distribution. The old value was wrong.
- **What to do:** if this applies to you, recompute infoVal. Old `.Rdata` files do
  not contain `nscales`, so the fix warns rather than guessing — regenerate the
  stimulus set, or pass the value explicitly.

**2. infoVal, if you passed `force_gen_ref_dist = TRUE` (issue #113)**

The flag was silently ignored whenever a `reference_norms` vector already existed
in the `.Rdata` file, so the reference distribution was reused rather than
recomputed.

- **Affected:** anyone who passed `force_gen_ref_dist = TRUE` expecting a fresh
  simulation. You received the cached distribution instead.
- **Not affected:** the default `force_gen_ref_dist = FALSE`, which was always
  meant to reuse the cache.
- **After the fix:** the flag regenerates the distribution, so infoVal will differ
  from a cached run. Because the reference distribution is simulated, it also
  varies slightly between regenerations — that is expected.

### Results cannot differ

These fixes only affect code paths that previously errored out, so no previously
obtained result changes:

- **Tibble inputs to `generateCI()`** (issues #70, #123) — previously failed with
  `arguments must have same length`.
- **The `mask` argument** — previously failed on R >= 4.2 with
  `the condition has length > 1`, and was restricted to 512px stimuli by a
  hardcoded size check.
- **Base images not matching `img_size`** (issue #124) — previously failed inside a
  parallel worker with `non-conformable arrays`.
- **Pre-0.3.3 `.Rdata` files** using the `sinusoids`/`sinIdx` layout — the
  backward-compatibility path was unreachable and always errored.
- **`simulateNoiseIntensities()`** — errored on every call.
- **Arguments overwritten by `load()`** — only reachable if you renamed or moved an
  `.Rdata` file between calls, in which case the old code either errored with
  `cannot open the connection` or wrote the reference distribution back to the file's
  former path. The `ncores` half affects speed only. No infoVal changes.

### Changed, but below any scale that can matter: the patch average

`generateNoiseImage()` now averages patch layers with `rowMeans(..., dims = 2)`
rather than `apply(..., 1:2, mean)`. These compute the same quantity but sum in a
different order, so they are **not bit-identical**: they differ by about one unit
in the last place, ~1e-19 in absolute terms on pixel values of order 0.01.

This is a different class of thing from the `rbernoulli` case below. That one would
have changed the random *stream* — a large, systematic divergence. This is
floating-point summation order, and it was checked against an independent oracle
(the average written as an explicit triple loop, using neither `apply()` nor
`rowMeans()`) across noise types, spatial scales and seeds: both the old and new
forms sit ~5.6e-17 from that oracle. Neither is "more correct" than the other.

- **Affected:** nothing you would report. No CI, z-map, infoVal or scaling decision
  changes at any precision a paper prints, and the difference is far below the
  noise you would get from re-running with a different seed.
- **At the configuration the golden master pins**, the results came out bit-identical.
- It is recorded here anyway, because the standard for this package is that any
  numeric change is written down rather than discovered later by someone re-running
  a five-year-old script.
- `tests/testthat/test-generateNoiseImage.R` now pins the new implementation against
  the original `apply()` form across 12 configurations, so it cannot drift further.

### Deliberately unchanged: the random number stream

`purrr::rbernoulli(n, p)` (deprecated) is internally `runif(n) > (1 - p)`. The obvious
modern replacement, `rbinom(n, 1, p)`, draws from the stream differently — verified across
150 seed/probability combinations that the two diverge. Swapping it in would have silently
changed every reference distribution, and therefore every infoVal, computed from a given
seed. The `runif` form is used instead, verified bit-identical to the old behaviour.

### Unchanged and verified

- **The infoVal formula itself.** rcicr already implements the corrected version
  from the erratum to Schmitz et al. (2019): the Euclidean norm, with the *k*
  constant supplied by R's `mad()` (`constant = 1.4826`). This was fixed years
  ago; it is now covered by a regression test so it cannot silently drift.
- **The noise basis** (`generateNoisePattern()`), the stimulus parameter draw, and
  all four scaling methods, for the default configuration.

---

*Older changes are recorded in `ChangeLog`.*
