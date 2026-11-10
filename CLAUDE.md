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

A fresh container has **no R at all**. The egress proxy allows or refuses each host separately, and the list changes, so **test a host before relying on it** (`curl -sS -o /dev/null -w '%{http_code}\n' https://<host>/`). On 2026-11-10, `cloud.r-project.org`, `cran.r-project.org` and `github.com` answered and `download.packages()` from `cloud.r-project.org` worked; `cran.rstudio.com`, `packagemanager.posit.co` and `bioconductor.org` returned `403`. A `403` is an organization egress policy, not a fault: report it, do not route around it, and never disable TLS verification or unset `HTTPS_PROXY`.

The setup script needs none of those hosts. Ubuntu's `universe` pocket has a Debian build of every dependency this package and its pinned v1.0.1 gate reference need, except one:

```sh
sudo bash tools/setup-container-r.sh   # toolchain + rcicr; safe to re-run; apt needs root
```

After that the ordinary workflow applies, including the release gate and `R CMD check --as-cran` with the PDF manual. With R 4.5.2 on Ubuntu 26.04 that check leaves one NOTE only this sandbox gives, `unable to verify current time`. R 4.3.3 also NOTEd `Uses the superseded package: 'doSNOW'`, from a hard-coded list that win-builder's R did not flag for 1.5.0. Because `bioconductor.org` is refused, the CRAN incoming feasibility check does not run (`NB: need Internet access to use CRAN incoming checks`); win-builder covers it. Measured here: `testthat::test_local()` gives 2069 passing and 1 skip, as uid 1000. The skip is the test that needs a case-insensitive file system.

Three things to know before trusting or editing that script:

- **The spatstat packages keep their dots**: `r-cran-spatstat.explore`, not `r-cran-spatstat-explore`. The usual dash convention reports MISSING and suggests, wrongly, that CRAN is needed after all.
- **`yesno` is the one package without a Debian build**, so the script installs a stub. rcicr itself no longer imports it (#392); the gate's reference versions, v1.0.1 through 1.5.0, still do, and reach it only from a prompt their empty lookup table can never trigger. The stub raises an error instead of answering, so nothing can quietly come to depend on a stubbed reply.
- **The v1.0.1 reference's own imports are installed too** (`raster`, `sp`, `ggplot2`, `plyr` and the others this package has since dropped), so the gate never needs `--install-deps` and never reaches for CRAN. `raster` also pulls in the GDAL/GEOS/PROJ system libraries it links against.

Two differences from CI will mislead you if you forget them. A session may run as **root** (check with `id -u`). Then a file made read-only is still writable, and the one test that depends on that skips instead of running. And the full release gate needs about 20 minutes per reference plus about 1.5 GB of RAM at 512px, so prefer running it in CI: a `workflow_dispatch` of `reproducibility.yaml` runs the **full** battery against both references, because the `--quick` choice depends on a `pull_request` event.

### Isolating one change from what its base branch already carries

`--ref` takes any git revision, not just a tag, and `EXPECTED` entries are filtered to those naming the resolved reference. A branch or SHA therefore matches **no** entry, and every difference is reported:

```sh
Rscript tools/compare-release-output.R --ref="$(git rev-parse origin/main)"
```

That answers what the release runs cannot. An entry without a `check` predicate excuses *any* deviation in its output. So where `main` already deviates from a tag for a documented reason, a new regression in the same output would be absorbed by that entry instead of reported. Against `main` itself nothing is excused, which turns "the gate is green" into "this change moves no number". Expect `0 expected deviations`; anything else is caused by the change.
