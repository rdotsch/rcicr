# CRAN comments

## Response to your review of 1.2.1

Thank you for the review. This is version 1.2.2, addressing each point.

**"Please omit the redundant 'Functions to' at the start of your description."** Done.

**References in the description field.** Added, in the requested form:

> ... For the method see Dotsch and Todorov (2012) <doi:10.1177/1948550611430272>; for a
> practical primer see Brinkman, Todorov and Dotsch (2017)
> <doi:10.1080/10463283.2017.1381469>.

**"Please write TRUE and FALSE instead of T and F."** Done — all 11 occurrences in `R/`
are gone, including the four that were public argument defaults and therefore visible in
`man/generateCI.Rd` and `man/plotZmap.Rd`, the two files you quoted. The values are
unchanged.

**"Some code lines in examples are commented out."** `generateNoiseImage()`'s example was
three commented-out lines and is now working code. On `generateNoisePattern.Rd` we could
not find a commented example — its `\examples` section contains the single live call
`generateNoisePattern(256)`, in the 1.2.1 tarball as well as in this one. If we have
misread which line you meant, please say and we will fix it.

**"`\dontrun{}` should only be used if the example really cannot be executed... Please
unwrap the examples if they are executable in < 5 sec."** The one `\dontrun{}`
(`simulateNoiseIntensities`) is gone, and so are all eight `\donttest{}` wrappers. **Every
example in the package now runs**, the complete set in about nine seconds. The wrappers
were inherited from a time when the examples ran at full stimulus size; they already ran at
32 pixels over six trials and were simply never revisited.

**"Please ensure that your functions do not write by default... Please omit any default
path in writing functions."** Done, and this is the largest change in the release. Every
argument naming a directory to write to — `stimulus_path`, `targetpath`, `zmaptargetpath` —
has lost its default. They were `./stimuli`, `./cis` and `./zmaps`. A call that would
previously have written to one of those now stops with an error naming the argument to
supply, so no user gets files somewhere new without noticing. This is a breaking change and
`NEWS.md` opens with it and the one-line migration.

Fixing it also removed a related violation you did not see:
`generateReferenceDistribution2IFC()` re-derives its noise basis in memory and wrote nothing,
but created an empty `./stimuli` on every call regardless.

All examples, tests and vignettes write only to `tempdir()`.

**"Please make sure that you do not change the user's options, par or working directory...
-> R/plotZmap.R."** Fixed. `plotZmap()` now captures the previous value at the point of
change and restores it through an immediate `on.exit()`, and closes its PNG device the same
way so a failure mid-plot cannot leak it.

One deviation from the pattern in your mail, deliberately: we restore with
`oldpar <- par(mar = ...)` / `on.exit(par(oldpar))` rather than
`par(no.readonly = TRUE)`. `par()` returns the previous values of exactly the parameters
being set, which is what needs restoring here; `no.readonly = TRUE` additionally captures
derived parameters such as `pin`, which the subsequent `plot.window()` call invalidates, so
restoring them raises `invalid value specified for graphical parameter "pin"`. This is the
second form shown in your mail (`oldpar <- par(mfrow = c(1,2)) ... par(oldpar)`).

**"inst/doc/reverse-correlation-walkthrough.R, please reset the par()."** Done. The
vignette now records `old_par <- par(no.readonly = TRUE)` in its setup chunk and restores
it in a final chunk, and its plotting helper restores `par()` directly instead of through
`on.exit()`.

## Submission type

This is a **resubmission of an archived package**, not a routine update.

`rcicr` was archived on CRAN on **2021-06-08**. The CRAN repository database records the
reason as:

> Archived on 2021-06-08 as email to the maintainer was undeliverable.

The archival was an administrative one — no policy violation and no code defect was
involved. The maintainer address on file had stopped working. The maintainer address is
now `rdotsch@gmail.com`, it is monitored, and it is the address this submission is sent
from.

The package has been actively maintained since. This version fixes a number of bugs
found by a systematic source review, adds a test suite (previously there was none) and
continuous integration, and removes 13 unused dependencies.

**On the version number.** The archived CRAN version is 0.3.4.1 and this submission is
1.2.2. Nothing is missing: 1.0.1, 1.1.0, 1.2.0 and 1.2.1 exist only as GitHub releases,
made while the package was off CRAN, and none of them was published anywhere a user could
install it from with `install.packages()`. `NEWS.md` carries their entries, so the record
between 0.3.4.1 and 1.2.2 is complete.

## Test environments

* local: Ubuntu 24.04, R 4.3.3 — `R CMD check --as-cran`, with
  `_R_CHECK_CRAN_INCOMING_=TRUE` and `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`
* GitHub Actions: ubuntu-latest on R release and R devel, macos-latest on R release,
  windows-latest on R release — all green
* win-builder, R-devel — R Under development (unstable) (2026-08-06 r90366 ucrt): 1 NOTE
  (the incoming feasibility note below). Install 39s, check 408s; examples 24s, tests 80s.
* win-builder, R-release — R 4.6.1 (2026-06-24 ucrt): 1 NOTE (the same one). Install 37s,
  check 391s; examples 24s, tests 75s.
* R-hub, all on R-devel — Linux, Windows, macOS:
  <!-- TODO: fill in from the 1.2.2 run -->. R-hub does not run the CRAN incoming
  feasibility check, which is why the note below does not appear there.

All of the above were run against the tarball being submitted.

## R CMD check results

0 errors | 0 warnings | 2 notes

The second is local to our build machine; on win-builder only the first appeared.

### NOTE 1 — CRAN incoming feasibility

```
New submission
Package was archived on CRAN
```

Expected for a reinstatement; see the explanation above.

```
Found the following (possibly) invalid URLs:
  URL: https://medium.com/@rondotsch/reverse-correlation-image-classification-using-r-a0701648fb0/
    From: README.md
    Status: 403
```

This URL is valid and resolves normally in a browser. `medium.com` returns 403 to
requests originating from datacenter networks — this is a block by network origin rather
than by user agent, so it reproduces from CI machines and not from an ordinary desktop.
Consistent with that, it appears on our build machine but on **neither** win-builder run
for this version, nor in CRAN's own incoming pretest of 1.2.1.

It points to the original tutorial walkthrough of the method. That walkthrough has been
ported into this release as the `reverse-correlation-walkthrough` vignette, so it is no
longer the primary documentation. The single remaining reference, in `README.md`, is a
"further reading" pointer to the historical copy, which has nine years of inbound links
from published papers. We are happy to remove the link if you prefer — it is now a
one-line change, since the content itself ships with the package.

win-builder additionally reports, under the same note:

```
Possibly misspelled words in DESCRIPTION:
  Brinkman (13:5)
  Dotsch (11:68, 13:27)
  Todorov (12:5, 13:15)
  psychophysical (10:33)
```

All four are correct. `Brinkman`, `Dotsch` and `Todorov` are the surnames of the authors of
the two references you asked us to add — `Dotsch` is also the maintainer's own name. They
appear only in the citations.

`psychophysical` is the standard adjective for psychophysics, the field this package's
method comes from, and is used in the sense of "psychophysical task" — the 2-image
forced-choice procedure the package generates stimuli for.

### NOTE 2 — future file timestamps

```
checking for future file timestamps ... NOTE
  unable to verify current time
```

This is our build machine being unable to reach the network time service used for the
check. It does not reproduce elsewhere.

## Downstream dependencies

There are none. The package has been off CRAN since 2021, so nothing currently depends
on it.

## Notes for the reviewer

* No example is wrapped in `\dontrun{}` or `\donttest{}` any more; all of them run. The
  ones that generate stimuli do so at 32x32 pixels over six trials, against a synthetic
  base image built in the example itself, so the whole set takes about nine seconds. A
  realistic run for a researcher is 512x512 pixels over several hundred trials, which is
  why the reduced parameters are there.
* `simulateNoiseIntensities()`'s example draws a boxplot, so running the examples
  non-interactively produces an `Rplots.pdf` in the check directory. That is R's default
  device rather than anything the package writes; no function in the package opens a file
  the caller has not named.
* The package uses `parallel`/`doParallel`. All defaults respect
  `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2 when that variable is set and
  `parallel::detectCores() - 1` otherwise, so no example, test or vignette uses more than
  two cores under check.
* Four test files call `skip_on_cran()`. They are development guards — a golden-master
  regression baseline, an end-to-end pipeline smoke test, a statistical signal recovery
  test, and a check that serial and parallel runs agree bit for bit — that protect against
  changes to numeric output during development. They are not needed to validate an
  installation, and they continue to run on every push in GitHub Actions.
