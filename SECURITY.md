# Security policy

**Keep this file under 600 words.**

## Reporting a vulnerability

Email the maintainer at the address in `DESCRIPTION` rather than opening
a public issue. Include a reproducible example and the rcicr version.

`rcicr` generates images and reads `.png`, `.jpg` and `.Rdata` files
supplied by the user. It opens no network connections and runs no code
from the files it reads. The realistic attack surface is malformed input
reaching an image decoder, or an untrusted `.Rdata` file:
[`load()`](https://rdrr.io/r/base/load.html) can execute code on read, a
property of R rather than of this package. **Only load `.Rdata` files
that you or your collaborators generated.**

## Supported versions

Only the most recent release. Fixes are not backported; upgrade to the
current version.

## Dependency posture

`.github/dependabot.yml` covers GitHub Actions and nothing else, on
purpose.

**The R dependencies cannot be watched, and pretending otherwise is
worse than the gap.** Dependabot has no CRAN ecosystem, and CRAN
publishes no security advisory database. Sonatype’s OSS Index via
`oysteR` was considered and rejected: its CRAN coverage is thin, so it
reports “clean” because its database is empty, not because the
dependencies are safe. The realistic dependency failure, an import being
archived, is caught by `R CMD check` on four platforms for every pull
request, every push to `main`, and weekly.

The real CVE surface is not in R code. It is in the C libraries that
compiled dependencies bind to (libpng via `png`, libjpeg via `jpeg`),
which the user’s operating system patches, not CRAN and not this
package.

**The actions do need watching.** Twelve of the thirteen were pinned to
floating major tags, so whoever controls those repositories could change
what runs in CI at any time. That is the one supply-chain surface here
with a real incident history, and the one Dependabot supports. Updates
arrive as a single grouped PR: `.github/workflows/` is not on the
reproducibility gate’s inert allowlist, so each separate bump would cost
a full CI round.
