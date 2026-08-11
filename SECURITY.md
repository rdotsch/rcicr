# Security policy

## Reporting a vulnerability

Email the maintainer at the address in `DESCRIPTION` rather than opening a public issue. A
report is most useful with a reproducible example and the rcicr version.

`rcicr` generates images and reads `.png`, `.jpg` and `.Rdata` files supplied by the user. It
opens no network connections and runs no code from the files it reads, so the realistic surface
is malformed input reaching an image decoder, or an untrusted `.Rdata` file — `load()` can
execute code on read, which is a property of R, not of this package. **Only load `.Rdata` files
you or your collaborators generated.**

## Supported versions

The most recent release. Fixes are not backported; upgrade to the current version.

## Dependency posture

`.github/dependabot.yml` covers `github-actions` and nothing else. The asymmetry is deliberate.

**The R dependencies cannot be watched, and pretending otherwise is worse than the gap.**
Dependabot has no CRAN ecosystem, and CRAN publishes no security advisory database. Sonatype's
OSS Index via `oysteR` was weighed and rejected: a scanner whose CRAN coverage is thin reports
"clean" because its database is empty, not because the tree is safe, and a green badge that
means nothing is worse than a known gap. The realistic dependency failure — an import getting
archived — is caught by `R CMD check` on four platforms on every push.

The real CVE surface is not in R code at all. It is in the C libraries the compiled
dependencies bind to — libpng via `png`, libjpeg via `jpeg` — and those are patched by the
user's operating system, not by CRAN, and not by anything this package does.

**The actions genuinely needed watching**: twelve of the thirteen in use were on floating major
tags, so whoever controls those repositories could change what runs in CI at any time. That is
the one supply-chain surface here with real incident history, and the one Dependabot supports.
Updates are grouped into a single PR, because `.github/workflows/` is deliberately not on the
reproducibility gate's inert allowlist and ungrouped bumps would each spend a full CI round.
