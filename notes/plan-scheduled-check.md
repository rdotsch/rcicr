# Plan: give `R CMD check` a schedule

## The gap

No workflow in `.github/workflows/` carries a `schedule:` trigger. All six fire on push to
`main`, pull requests to `main`, version tags, or manual dispatch:

```
$ grep -rn "schedule\|cron" .github/
.github/dependabot.yml:10:    schedule:
```

The only cron in the repository belongs to dependabot. So CI runs exactly as often as someone
pushes, and stops when they stop. A break caused from outside the repository — an R-devel
change, a dependency's API, a platform's graphics backend — is then invisible until CRAN's own
checks catch it, and CRAN's remedy is an archival deadline. This package has been archived
once: `ChangeLog` records the CRAN era up to the archived 0.3.4.1, and 1.3.0 was the
reinstatement.

## The change

Add to `.github/workflows/R-CMD-check.yaml`:

```yaml
on:
  push:
    branches: [main]
  pull_request:
    branches: [main]
  schedule:
    - cron: '0 5 * * 1'
  workflow_dispatch:
```

Weekly, Mondays 05:00 UTC. `schedule` always runs on the default branch, so no branch filter is
needed. `workflow_dispatch` comes with it so the same matrix can be run on demand — the check a
returning maintainer wants first, and currently only reachable by pushing.

**This workflow rather than a new one.** `ubuntu-latest (release)` and `ubuntu-latest (devel)`
are required checks matched by name and enforced by a ruleset; the file's own header warns that
renaming a job makes every PR unmergeable. Adding a trigger renames nothing and reuses the
four-platform matrix. A separate canary workflow would duplicate the matrix and drift from it.

**Weekly rather than monthly.** CRAN's notice period is measured in weeks, so a monthly cadence
can put first warning after the deadline that matters. Public-repository Actions minutes are
not metered, and GitHub notifies on failure only, so a green week costs nothing to read.

## What this rests on, and what it does not fix

Two documented GitHub behaviours, neither verifiable from here, both to confirm once a run has
fired:

- **Scheduled workflows in public repositories are disabled after 60 days with no repository
  activity**, with an email announcing it. dependabot's weekly branch pushes are activity of
  that kind, but `open-pull-requests-limit: 1` means an unmerged PR stops new ones — the
  quiet-repository case is exactly when the keepalive lapses. The honest description is a
  warning system whose own death is announced by mail, not a guarantee.
- **Failure notifications go to whoever last edited the cron syntax**, not the repository
  owner, and re-enabling a disabled workflow transfers them to whoever re-enabled it. The
  maintainer must be the committer of this line for the mail to reach him.

Both go in a comment beside the `schedule:` block rather than in `MAINTENANCE.md`, which is at
1795 words against its stated 1800 budget. The workflow file is where someone editing the cron
will read them, and this repository already keeps its load-bearing CI explanations there — the
warning about renaming a required job sits at the top of this same file.

## Most likely to fail

A scheduled run executes the two `if:`-gated steps (`Documentation is up to date`,
`CITATION.cff matches its sources`) on the same conditions as a push, since they key off the
matrix row rather than the event. Both regenerate from pinned generators (`roxygen2@7.3.1`,
`cffr@1.4.1`), so they are deterministic and should stay green — but they are the steps that
would turn a canary into a recurring false alarm, and the first scheduled run is what confirms
otherwise.

## Verification

`workflow_dispatch` is unavailable until the trigger is on the default branch, so before merge
the check is a YAML parse and a diff review; after merge, one manual dispatch confirms the
matrix still runs green and that the two gated steps behave on a non-push event.
