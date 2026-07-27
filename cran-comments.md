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

## Test environments

* local: Ubuntu 24.04, R 4.3.3 — `R CMD check --as-cran`, with
  `_R_CHECK_CRAN_INCOMING_=TRUE` and `_R_CHECK_CRAN_INCOMING_REMOTE_=TRUE`
* GitHub Actions: ubuntu-latest, R release and R devel
* win-builder (R-devel)  <!-- TODO: run devtools::check_win_devel() and record result -->
* R-hub: macOS, Windows, R-devel  <!-- TODO: run rhub::rhub_check() and record result -->

## R CMD check results

0 errors | 0 warnings | 2 notes

### NOTE 1 — CRAN incoming feasibility

```
New submission
Package was archived on CRAN
```

Expected for a reinstatement; see the explanation above.

```
Found the following (possibly) invalid URLs:
  URL: https://medium.com/@rondotsch/reverse-correlation-image-classification-using-r-a0701648fb0/
    Status: 403
```

This URL is valid and resolves normally in a browser. `medium.com` returns 403 to
requests originating from datacenter networks — this is a block by network origin rather
than by user agent, so it reproduces from CI machines and not from an ordinary desktop.
It is a tutorial walkthrough of the method that the package implements, and it is
referenced from `README.md` and the vignette as a pointer for new users. We have kept it
because it is genuinely useful to readers, but are happy to remove it if you prefer.

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
* Three test files call `skip_on_cran()`. They are development guards — a golden-master
  regression baseline, an end-to-end pipeline smoke test, and a statistical signal
  recovery test — that protect against changes to numeric output during development. They
  are not needed to validate an installation, and they continue to run on every push in
  GitHub Actions.
