Warning: truncated output (original token count: 16422)
Total output lines: 994

# rcicr (development version)

## Reproducibility impact

- **`generateStimuli2IFC()` now ignores an image's alpha channel when reading a base face.**
  Greyscale and RGB images, and the usual opaque PNG with constant alpha at the default contrast
  setting, are unaffected. A base face with varying transparency could previously look visibly
  distorted in the stored base, stimuli and rendered CI; the noise parameters and raw
  classification image are unchanged. The problem is apparent on visual inspection, so no
  results are known to be invalidated. If it matters for a past study, regenerate the stimuli
  from the original base image rather than rerunning `generateCI()` from its existing `.Rdata`.

## Documentation

- **`ChangeLog` is frozen at 1.0.1 and is no longer updated.** It keeps the record of
  everything before 1.1.0, including the whole CRAN era up to the archived 0.3.4.1.
  `NEWS.md` is the changelog from 1.1.0 onward and the only one R indexes, so
  `news(package = "rcicr")` is unaffected. Entries for 1.1.0 to 1.3.0, which only
  pointed here, were removed from it.

# rcicr 1.3.0 (2026-08-18)

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

- **`generateCI(mask = matrix(NA, 1, 1))` now reports the malformed mask instead of silently
  ignoring it.** The internal test for "was a mask supplied?" asked only whether the argument
  was a single `NA` — which a one-cell `NA` matrix, a one-element `array(NA)` and a `list(NA)`
  all are. Such a mask was mistaken for the `NA` default, discarded without a word, and an
  entirely **unmasked** classification image came back. A call that passes a mask and gets an
  unmasked CI is the failure worth catching early; it now stops in the same mask validation
  every other malformed mask reaches. The sentinel is now specifically an atomic scalar with no
  dimensions, so `mask = NA` still means "no mask" — it is the default, and every unmasked call
  relies on it — while `matrix(1, 1, 1)`, larger matrices, PNG paths and `NULL` are unchanged.

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

- **The `.Rdata` validation errors now name the version of rcicr that wrote the file.**
  `generateCI()` and `computeCumulativeCICorrelation()` already said which field a file was
  missing; they now also say where the file came from, which is what turns "this file has no
  `stimuli_params`" into "it predates the version that added it — regenerate the stimuli, or
  install that version". The version is read tolerantly, because the field cannot be taken at
  face value: `p$generator_version` is preferred over the top-level `generator_version`, which
  every release from 0.4.0 through 1.1.0 recorded as a hardcoded `0.4.0`, and a file that only
  claims `0.4.0` is reported as *unknown* rather than as 0.4.0. A file with no version field at
  all reports the absence and stops there — an absent field is equally a file older than 0.4.0,
  a truncated one, or one rcicr never wrote.

  `generateReferenceDistribution2IFC()`'s two warnings about a missing `nscales` or
  `noise_type` were reworded for the same reason. They said the file "was written by a version
  of rcicr that did not save" the field, which an absence does not establish; they now report
  what is missing, and keep the version as context rather than as a conclusion. The advice is
  unchanged, and was always right either way: regenerate the stimulus set with this version.

- **A stimulus file with gabor noise and no saved `sigma` now says so.**
  `generateReferenceDistribution2IFC()` assumes the historical default of 25 when a file
  predates 1.1.0 and lacks the field, which it has always done silently — unlike the loud
  warnings for a missing `nscales` or `noise_type`. For gabor noise that silence hid the same
  hazard those warnings exist for: `sigma` is what shapes the Gaussian mask, so guessing it
  wrong rebuilds the null on a different noise basis than participants saw, and the resulting
  InfoVal is wrong. On a 1.0.1 gabor stimulus set the reference norms move from
  0.681/0.689/0.680 at `sigma = 25` to 0.615/0.620/0.626 at `sigma = 10`.

  **Sinusoidal files are unaffected and stay silent**, which is the point of doing this by
  noise type rather than by field: `sigma` reaches the basis through `generateGabor()` alone,
  so for sinusoidal noise the norms are identical whatever it is, and a warning would be pure
  noise on the far more common legacy file. Nothing warns that did not previously produce a
  wrong answer, and no numeric output changes.

- **`plotZmap(mask = ...)` now accepts a mask with an alpha channel.** It previously required
  every colour channel of a multi-channel PNG mask to match exactly, so a greyscale-plus-alpha
  or RGBA mask whose alpha plane happened to differ from its colour planes was rejected — alpha
  carries no colour information and is now always ignored, matching `generateCI(mask = ...)`.
  A rectangular `zmap`/mask pair continues to work as before. `plotZmap(mask = ...)` also now
  rejects a mask that is neither a string nor a matrix with a clear error, instead of failing
  later inside `png::readPNG()` with an unrelated message.

- **`plotZmap(mask = NA)` now means "no mask", as it already does in `generateCI()`.** The two
  functions detected a supplied mask differently — `plotZmap()` asked only whether the argument
  was non-`NULL` — so the same sentinel meant opposite things: `generateCI()` read `NA` as "no
  mask" while `plotZmap()` passed it on and stopped with `The mask argument is neither a string
  nor a matrix!`. `NaN` behaved the same way. Both now render an unmasked z-map. No call that
  worked before changes: the inputs affected all raised that error.

- **`generateStimuli2IFC()` no longer saves `trial` in the `.Rdata` file.** It was the loop
  counter left over from stimulus generation — always equal to `n_trials`, which is already
  saved. Nothing in the package or the documented contract reads it. Existing `.Rdata` files
  that contain `trial` continue to work; the field is simply ignored on load.

## Reproducibility impact

- **Individual-CI PNGs written by `generateCI(save_individual_cis = TRUE)` carried the wrong
  participant's name, and are now named correctly.** The per-participant loop selects each
  participant's trials by *sorted* order and took the output filename from order of
  *appearance*. At each position where those two orders differ, the file in
  `<targetpath>/individual_cis` was given another participant's ID. The pixels were always
  right; only the names were wrong.

  **Affected:** a direct call to
  `generateCI(participants = ..., save_individual_cis = TRUE)` where the `participants`
  vector is not already in sorted order. Such a call produced correct images under incorrect
  filenames, so a figure published as participant `p2` may be someone else's classification
  image.

  **Do not assume tidy data was safe — the common case is affected.** Sorting is *lexical*,
  so text labels like `"p10"` sort before `"p2"`. A study whose participants appear in
  collection order `p1, p2, ..., p10` is therefore affected from the tenth participant
  onward, even though nothing about that data is out of order in any ordinary sense:

  | `participants`, in collection order | affected |
  |---|---|
  | `"p1"` … `"p9"` | no |
  | `"p1"` … `"p10"` or beyond | **yes** |
  | `"p01"` … `"p12"` (zero-padded) | no |
  | `"1"` … `"12"` (character) | **yes** |
  | `1` … `12` (numeric) | no |

  Unpadded text IDs and ten or more participants is the ordinary shape of a reverse
  correlation study, so if that describes yours and you used the call above, check rather
  than assume.

  **Not affected — and this is the documented route, so most analyses are in this group:**

  - **`batchGenerateCI()`, `batchGenerateCI2IFC()` and `generateCI2IFC()`.** None of them
    exposes `participants` or `save_individual_cis` at all; they split the data themselves
    and name each image from the grouping value directly, never reaching the code that was
    wrong. Per-participant classification images computed the way the package documentation
    has recommended since the CRAN era — `batchGenerateCI2IFC()`, and `batchGenerateCI()` in
    the current vignette — were correctly labelled throughout.
  - **Everything `generateCI()` returns.** The group classification image is a mean across
    participants and so does not depend on their order; the per-participant stack, and every
    z-map built from it, were internally consistent throughout. No returned value, and no
    number in any `.Rdata` file, ever changed.
  - **Any call that left `save_individual_cis` at its default `FALSE`.** No individual-CI
    file was written, so there was nothing to mislabel.
  - **Participants already in sorted order**, where the two orderings coincide.
  - **Every CRAN release.** rcicr was on CRAN from July 2014 until it was archived in June
    2021; the last release, 0.3.4.1 (July 2016), predates the `save_individual_cis` option
    entirely. `install.packages('rcicr')` never produced a mislabelled file at any point.

  **If you have published from individual-CI images, check your analysis scripts.** Two steps, neither of which requires re-running anything. The [individual-CI filename advisory](https://rdotsch.github.io/rcicr/articles/rcicr-individual-ci-advisory.html) carries the fuller version of this, including what the mislabelling does to a second-stage analysis and how much of an effect survives it; where the two differ, that page is kept current and this section is not.

  *Step 1 — did you use the relevant call at all?* Look for an `individual_cis/` directory in your output. Nothing else in the package creates it, so finding it means the direct call ran. For output made before 1.3.0, continue to step 2. If your per-participant images came from `batchGenerateCI()` or `batchGenerateCI2IFC()`, those cannot reach the defect.

  Check the output rather than the script. `save_individual_cis` is the sixth formal of
  `generateCI()`, so a call can set it positionally —
  `generateCI(stim, resp, "face", rdata, pids, TRUE, ...)` — and write those files without
  the argument name appearing anywhere; `do.call()` hides it the same way. Searching your
  scripts and finding nothing does not clear an analysis.

  Not finding the directory clears it only if the output still sits where the run left it. Renamed, archived or moved, it will not answer to that name, and the individual files are then recognisable but not conclusive: `ci_<participant>.png` is what the affected call writes, and also what `generateCI(save_as_png = TRUE, filename = "<participant>")` writes for a group CI one directory up.

  *Step 2 — if the directory is there, was your `participants` vector in sorted order?*
  Note "sorted" means lexically sorted for text labels, per the table above — sorting your
  data frame by a `"p1"`-style column does not make it safe once you reach ten participants.
  Where you still have the vector, one expression settles it:

  ```r
  identical(unique(participants), sort(unique(participants)))   # TRUE = names were correct
  ```

  Two conditions on that `TRUE`, both of which change what `sort()` returns: use the original `participants` object rather than a character copy, since a factor carries its own level order and `sort()` follows it; and run it under the collation the original analysis ran under, since identifiers mixing case, accents or punctuation sort differently by locale. Numeric identifiers are not exposed to the second, numeric sorting being locale-independent, and `p1`/`p12` are not realistically exposed either. Lowercase ASCII is no guarantee in general, though, since some locales collate letter pairs as units (Czech `ch`, Danish `aa`).

  If the output is long gone and only the script survives, search it for `individual` rather
  than the full argument name, and read by hand any `generateCI()` call passing six or more
  arguments positionally.

  **Recovering existing output is a rename, not a re-run.** The mapping is exact: the file
  named `unique(participants)[i]` holds the classification image of
  `sort(unique(participants))[i]`.

  ```r
  old <- unique(participants)         # the names the files were given
  new <- sort(unique(participants))   # the participants they actually hold
  ```

  **Every release from 1.0.1 (January 2023) through 1.2.3 (August 2026) carries the
  defect**, as does every install from the development branch from 15 August 2017 onward —
  the date it was introduced, in the commit that first made individual-CI saving work at
  all. There were no releases between 2017 and 2023, so a copy obtained in those years came
  from the branch and is affected too.

  **Getting it from GitHub took a deliberate step, until December 2021.** The default branch
  was `master`, which carried an older repository layout with no `DESCRIPTION` at its root, so
  `install_github('rdotsch/rcicr')` did not install anything at all — it failed. The affected
  code was on the `development` branch, reached only by asking for it:
  `install_github('rdotsch/rcicr', ref = 'development')`, which the README of the time
  labelled "AT YOUR OWN RISK" and specifically advised against for analyses meant for
  publication. `development` became the default branch on 2021-12-28, and from then a plain
  install returns affected code.

  If you need to work out which version a stored analysis actually used, the version number
  alone will not tell you: `0.4.0` sat on the development branch for five years, with the
  defect entering partway through. A note listing every release, its date, whether it came
  from CRAN or GitHub, and what "the latest version" gave you in each window is at
  <https://github.com/rdotsch/rcicr/blob/main/notes/individual-ci-mislabelling.md>.

- **`computeCumulativeCICorrelation()` with a masked `targetci`** now returns numeric
  correlations where it previously returned all-`NA`. No existing analysis could have used
  the old result — it carried no information — but code that checked for `NA` on the returned
  curve will see a different answer. Unmasked targets are bit-identical.

## Bug fixes

- **`generateCI(save_individual_cis = TRUE)` names each PNG for the participant it was
  actually computed from.** The loop selects participants by sorted level order and the
  filename was taken from order of appearance, so the two disagreed for any data not already
  sorted by participant. See "Reproducibility impact" above for how to tell whether files you
  already have are affected, and how to correct them without recomputing anything.

- **`computeCumulativeCICorrelation()` now returns real correlations when `targetci` was
  generated with a mask.** `generateCI()` stores `NA` in every pixel a `mask` excludes,
  and the correlation was taken over all pixels (`use = "everything"`), so a single masked
  pixel made every point on the curve `NA` however strong the signal. Correlations are now
  computed over the unmasked pixels only. A fully masked target (every pixel `NA`) still
  returns an all-`NA` curve — there are no complete pairs to correlate.

- **`plotZmap(col = ...)` works.** Supplying a palette is how `?plotZmap` has always told
  you to change the colours, and doing it stopped the call with `formal argument "col"
  matched by multiple actual arguments` before anything was drawn — the function passed its
  own `col` alongside yours. This affected every released version; the argument is now taken
  as an override, and is used for the colour bar as well as the map.

- **A base image with no contrast no longer becomes an all-`NaN` base image.** Under
  `maximize_baseimage_contrast = TRUE`, the default, `generateStimuli2IFC()` rescales
  with `…6422 tokens truncated…the value is `1`/`TRUE`
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

*Changes before 1.1.0 are recorded in `ChangeLog`, which is frozen at 1.0.1 and covers
the package back to 2014, including its CRAN era.*
