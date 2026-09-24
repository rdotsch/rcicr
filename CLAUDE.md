<!--
This file exists only so that Claude Code loads this repository's conventions.
Claude Code reads CLAUDE.md and does not read AGENTS.md:
https://code.claude.com/docs/en/memory#agents-md

AGENTS.md remains the single source of truth for the repository's conventions.
Put those there, not here. What belongs here is the little that is true only of
the Claude Code container and would be noise in a file other agents read --
currently the bootstrap below, and nothing else.

An @-import is used rather than `ln -s AGENTS.md CLAUDE.md` because symlinks
require Administrator privileges or Developer Mode on Windows, and this package
has Windows contributors and a Windows CI runner.

Block-level HTML comments are stripped before the file is loaded, so this note
costs no context. To confirm the file is loading, run /context in a session and
look for CLAUDE.md under "Memory files".
-->

@AGENTS.md

## Getting R in a Claude Code container

A fresh container has **no R at all**, and the egress proxy refuses every CRAN mirror and every GitHub tarball host: `403` on `CONNECT` to `cloud.r-project.org`, `cran.r-project.org`, `cran.rstudio.com`, `packagemanager.posit.co`, `github.com` and `codeload.github.com`. So `install.packages()` and `remotes::install_github()` are both dead ends. That is an organization egress policy, not a fault: report it, do not route around it, and never disable TLS verification or unset `HTTPS_PROXY`.

Ubuntu's `universe` pocket *is* reachable. It has a Debian build of every dependency this package and its pinned v1.0.1 gate reference need, except one:

```sh
bash tools/setup-container-r.sh      # toolchain + rcicr; safe to re-run
```

After that the ordinary workflow applies, including the release gate. Measured here: `testthat::test_local()` gives 1141 passing and 1 skip.

Three things to know before trusting or editing that script:

- **The spatstat packages keep their dots**: `r-cran-spatstat.explore`, not `r-cran-spatstat-explore`. The usual dash convention reports MISSING and suggests, wrongly, that CRAN is needed after all.
- **`yesno` is the one package without a Debian build**, so the script installs a stub. `computeInfoVal2IFC()` only reaches the real one from a single `if (interactive())` branch, behind a `ref_lookup` whose rows have been commented out since 2018, so a batch run can never reach the prompt. The stub raises an error instead of answering, so nothing can quietly come to depend on a stubbed reply.
- **The v1.0.1 reference's own imports are installed too** (`raster`, `sp`, `ggplot2`, `plyr` and the others this package has since dropped), so the gate never needs `--install-deps` and never reaches for CRAN. `raster` also pulls in the GDAL/GEOS/PROJ system libraries it links against.

Two differences from CI will mislead you if you forget them. Sessions run as **root**, so a file made read-only is still writable, and the one test that depends on that skips instead of running. And the full release gate needs about 20 minutes per reference plus about 1.5 GB of RAM at 512px, so prefer running it in CI: a `workflow_dispatch` of `reproducibility.yaml` runs the **full** battery against both references, because the `--quick` choice depends on a `pull_request` event.

### Isolating one change from what its base branch already carries

`--ref` takes any git revision, not just a tag, and `EXPECTED` entries are filtered to those naming the resolved reference. A branch or SHA therefore matches **no** entry, and every difference is reported:

```sh
Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"
```

That answers what the release runs cannot. An entry without a `check` predicate excuses *any* deviation in its output. So where `main` already deviates from a tag for a documented reason, a new regression in the same output would be absorbed by that entry instead of reported. Against `main` itself nothing is excused, which turns "the gate is green" into "this change moves no number". Expect `0 expected deviations`; anything else is caused by the change.
