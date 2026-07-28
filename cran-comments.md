# CRAN comments

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
1.2.1. Nothing is missing: 1.0.1, 1.1.0 and 1.2.0 exist only as GitHub releases, made while
the package was off CRAN, and none of them was published anywhere a user could install it
from with `install.packages()`. `NEWS.md` carries their entries, so the record between
0.3.4.1 and 1.2.1 is complete.

1.2.1 follows 1.2.0 by a short interval and contains no user-facing change. 1.2.0 did not
pass `R CMD check` on macOS: a test in the package's own suite asserted properties of a
rendered PNG that belong to the graphics device rather than to what was drawn, and those
differ between macOS and Linux. No released function was affected and nothing the package
computes changed, but the tree could not pass its own checks on a platform CRAN builds for,
so it is not the tree we are asking you to accept.

## Test environments

* local: Ubuntu 24.04, R 4.3.3 — `R CMD check --as-cran`, with
  `_R_CHECK_CRAN_INCOMING_=TRUE` and `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`
* GitHub Actions: ubuntu-latest on R release and R devel, macos-latest on R release,
  windows-latest on R release — all green
* win-builder, R-devel  <!-- TODO: re-run against the 1.2.1 tarball and record result -->
* win-builder, R-release  <!-- TODO: re-run against the 1.2.1 tarball and record result -->
* R-hub: Linux, Windows, macOS  <!-- TODO: run against 1.2.1 and record result -->

The win-builder results previously recorded here were for the **1.2.0** tarball (R-devel
2026-07-26 r90304 ucrt and R 4.6.1, 1 NOTE each — the incoming feasibility note below).
They are not carried forward: the submitted tarball has changed, so they are re-run rather
than assumed still to hold.

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

It points to the original tutorial walkthrough of the method. That walkthrough has been
ported into this release as the `reverse-correlation-walkthrough` vignette, so it is no
longer the primary documentation. The single remaining reference, in `README.md`, is a
"further reading" pointer to the historical copy, which has nine years of inbound links
from published papers. We are happy to remove the link if you prefer — it is now a
one-line change, since the content itself ships with the package.

win-builder additionally reports, under the same note:

```
Possibly misspelled words in DESCRIPTION:
  psychophysical (10:51)
```

`psychophysical` is spelled correctly. It is the standard adjective for psychophysics, the
field this package's method comes from, and it is used in the sense of "psychophysical
task" — the 2-image forced-choice procedure the package generates stimuli for.

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

* Examples that generate stimuli are wrapped in `\donttest{}` and run at a much reduced
  image size and trial count. A realistic run for a researcher is 512x512 pixels over
  several hundred trials, which is far too slow for a check.
* The package uses `parallel`/`doParallel`. All defaults respect
  `_R_CHECK_LIMIT_CORES_`: `default_ncores()` returns 2 when that variable is set and
  `parallel::detectCores() - 1` otherwise, so no example, test or vignette uses more than
  two cores under check.
* Four test files call `skip_on_cran()`. They are development guards — a golden-master
  regression baseline, an end-to-end pipeline smoke test, a statistical signal recovery
  test, and a check that serial and parallel runs agree bit for bit — that protect against
  changes to numeric output during development. They are not needed to validate an
  installation, and they continue to run on every push in GitHub Actions.
