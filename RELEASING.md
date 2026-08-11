# Releasing rcicr

The checklist, and why it is ordered the way it is. Maintainer-facing: contributors want
`CONTRIBUTING.md`. `NEWS.md` records what each release changed, `DECISIONS.md` why the package
behaves as it does.

---

## The checklist

Releases are squash-merged onto `main` and marked with a tag; there is no `develop` branch.
The version in `DESCRIPTION` carries a `.9000` suffix between releases, and the release
commit is the one that drops it.

**Step 1 is the reproducibility gate, and a release does not go out while it is red.**

```sh
Rscript tools/compare-release-output.R                 # vs v1.0.1 -- the published baseline
Rscript tools/compare-release-output.R --ref=v1.1.0    # vs the last release
```

Both runs are required, and they answer different questions. The first asks whether this
tree still produces the numbers that are *in the literature* — it stays pinned at v1.0.1,
the last version on CRAN, and never advances, because a reference that moves forward with
each release lets a tree drift away from those numbers one tolerated epsilon at a time. The
second asks whether anything broke since the last release, and covers more ground: calls
that used to crash (masks, small z-maps, undecorated z-maps) have a comparable value in
v1.1.0 and none at all in v1.0.1.

Exit status is `0` for a clean run, `1` for a difference that is not accounted for, and `2`
if the comparison could not be run at all. On a difference there are exactly two honest
outcomes:

- **it was not intended** — fix it; or
- **it was intended** — write it up in `NEWS.md` under "Reproducibility impact", saying who
  is affected and what to do, and add an entry to `EXPECTED` in the script naming the
  reference, the affected output, the reason, and that `NEWS.md` heading. The script fails
  if an `EXPECTED` entry stops firing (so the list cannot rot) and if it fires while
  `NEWS.md` says nothing (so a deviation cannot be quiet).

Silently widening a tolerance, deleting a configuration from the battery, or skipping the
run is none of those things. The whole point is that a change to researchers' numbers can
only get out deliberately and loudly.

### 1. Prepare the release PR

Branch from `main` (`release-X.Y.Z` by convention) and make these five edits together:

```sh
Rscript tools/compare-release-output.R              # both gate runs green
Rscript tools/compare-release-output.R --ref=vA.B.C # the previous release
R CMD build . && R CMD check --as-cran rcicr_*.tar.gz
```

- **`NEWS.md`**: rename `# rcicr (development version)` to `# rcicr X.Y.Z (YYYY-MM-DD)`.
  Until you do, none of this release's entries are in the news database at all — R indexes
  only sections under a heading it can read a version from. Keep version numbers *out* of
  `##` section headings, or the whole file stops parsing and `R CMD check` NOTEs.
- **`DESCRIPTION`**: drop the `.9000`. This is also what switches CI from the everyday
  `--quick` gate to the full battery, so it has to happen **before** the PR is opened, not
  after review.
- **`ChangeLog`**: add a dated pointer entry deferring to `NEWS.md`, never a duplicate — two
  changelogs both claiming authority is worse than one that defers.
- **`CITATION.cff`**: regenerate it, after the `DESCRIPTION` edit above, with
  `cffr::cff_write("DESCRIPTION", dependencies = FALSE, gh_keywords = FALSE)` — it carries the
  version. `ubuntu-latest (release)` regenerates and compares, so it fails until you do.
- **`cran-comments.md`**: re-read every claim in it rather than trusting it. It has gone
  stale *within a single day* before, by describing a URL reference that another PR had
  just removed, one step short of a CRAN reviewer reading it. Leave the test-environment
  results for step 2; everything else is written now.

Open the PR. The full battery runs on it, takes ~20 minutes, and is a required check.

### 2. Run the external checks — on the release branch, before merging

win-builder and R-hub run against the release branch, and their results go into
`cran-comments.md` **in the same PR**, so the release commit carries its own evidence.

```sh
R CMD build .                                # at the repo root, never in a worktree
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-devel/
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-release/
```

Then trigger R-hub from the Actions tab against the release branch. Two NOTEs are expected
locally and explained in `cran-comments.md` — CRAN incoming feasibility (new submission of
an archived package) and future file timestamps (no network clock in a sandbox).

This ordering works because **`cran-comments.md` is `.Rbuildignore`d and never reaches the
tarball**, so writing the results into it cannot invalidate the checks those results
describe — the tarball is unchanged. And `DESCRIPTION` on the release branch already carries
the clean version, so the checks see the exact version string that will be submitted.

Neither check is run casually: win-builder mails the maintainer, and both put the tarball in
front of a third party.

- **win-builder** — `devtools::check_win_devel()` and `devtools::check_win_release()` do the
  same thing as the `curl` lines above. Each FTPs the tarball to `win-builder.r-project.org`
  and mails a check log to the maintainer address within the hour; `226` means the transfer
  completed. The server **will not overwrite an existing upload**: a second attempt at the
  same filename fails with a bare `Failed FTP upload: 550`, which means "already queued", not
  "rejected". Confirm with `curl --list-only ftp://win-builder.r-project.org/R-devel/` before
  re-uploading, or you will conclude the upload failed when it succeeded. The mail carries a
  result URL; `curl -s https://win-builder.r-project.org/<key>/00check.log` fetches the full
  log, which beats transcribing the one-line summary. The pages expire after ~72 hours, so
  `cran-comments.md` is the durable copy.
- **R-hub** — `.github/workflows/rhub.yaml` is the stock R-hub v2 workflow, so the check runs
  on this repository's own Actions rather than uploading anywhere. Trigger it from the
  Actions tab ("R-hub" → Run workflow), which needs nothing but repository access.
  `rhub::rhub_check()` does the same over the API but requires a GitHub PAT with the `repo`
  scope, stored via `gitcreds::gitcreds_set()` — `rhub::rhub_doctor()` reports whether yours
  has it. Because the workflow is `workflow_dispatch`-only, it must be merged to the
  **default branch** before it can be triggered at all, and R-hub wants the file on the
  default branch *and* on whichever branch is being checked.

  To read the results, **download the artifacts — do not use `gh run view --log`**, which
  truncates long jobs: at 1.2.1 the 38-minute macOS job returned 1548 lines ending four
  minutes in, with no `Status:` line anywhere. `gh run download <run-id> -D <dir>` and read
  `*/rcicr.Rcheck/00check.log`. Expect `Status: OK` where win-builder reports 1 NOTE — R-hub
  does not run the CRAN incoming feasibility check, so that is agreement, not a discrepancy.

  **`RHUB_TOKEN` is deliberately unset, and nothing is missing.** The stock workflow passes
  `secrets.RHUB_TOKEN` to four actions, which makes it look like a prerequisite. It is an
  optional slot for your own PAT to reach private repositories — not a credential R-hub
  issues — and as of `r-hub/actions@v1` none of those four actions references the input in
  any step. This package is public; an unset secret expands to an empty string and the jobs
  run normally.

### 3. Merge and tag

```sh
gh pr merge <n> --squash --delete-branch
git checkout main && git pull && git fetch --prune
git tag -a vX.Y.Z -m "rcicr X.Y.Z" && git push origin vX.Y.Z
gh release create vX.Y.Z --title "rcicr X.Y.Z" --notes-file <(...)   # NEWS section
```

The tag push re-runs the full gate against the released tree, which is the copy that
matters. Squash-merge, always: one commit per PR on `main`, and the squash message is where
measurements and rejected alternatives have to live, because the branch commits disappear.

**Confirm the squash did not change the tree** the external checks ran on:

```sh
git diff <pr-head-sha> vX.Y.Z --stat    # must be empty
```

Empty output is what lets step 2's results stand for the tagged tree.

### 4. Submit to CRAN

```sh
git switch --detach vX.Y.Z      # never build from main HEAD, never in a worktree
R CMD build .
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
git switch -                    # back to where you were
```

Building from the tag is what keeps the development suffix out of the submitted version:
`Version contains large components` is only a CRAN blocker if the *tarball* carries it.

Submit at <https://cran.r-project.org/submit.html>, pasting the body of `cran-comments.md`
into the "Optional comment" field — that file is `.Rbuildignore`d and reaches CRAN only this
way, never inside the tarball. `devtools::submit_cran()` does the same thing.

**Paste the tag's copy, not `main`'s** — `git show vX.Y.Z:cran-comments.md`. Submission can
trail tagging by weeks, and `main` moves in between; its copy would describe check results
from a tree that is not the tarball, and nothing in the submission would contradict it.

**Ron submits to CRAN personally.** CRAN emails the maintainer address to confirm, and for a
package archived over an undeliverable address, that confirmation *is* the point of the
resubmission.

### 5. Reopen development

Bump `DESCRIPTION` to `X.Y.Z.9000` and start a fresh `# rcicr (development version)` heading
in `NEWS.md`. `main` is protected, so this goes through a small PR like anything else.

`.github/workflows/reproducibility.yaml` runs the gate on every PR to `main` (`--quick`,
which skips the 512px configuration) and in full on release PRs, version tags and manual
dispatch. The `compare` job is a **required status check** on `main`, alongside
`R CMD check` on release and devel, so a red gate blocks the merge rather than merely
reporting. That is a repository ruleset, not something the workflow can enforce on itself.


---

## Why the procedure is shaped this way


The conventions and their rationale are in `AGENTS.md` — trunk-based with tags and no
`develop` branch, `.9000` on `main` between releases, tarballs built from the tag, squash
merges, delete merged branches, `NEWS.md` ordered largest-impact first. Two things learned
the hard way that those entries do not cover:

- **Squash merging is what makes a bad commit cheap.** One branch here carried ~75,000 lines
  of committed `R CMD check` artifacts across three commits; `main` only ever saw the final
  tree, so no history rewrite was needed. (Both patterns are now in `.gitignore` and
  `.Rbuildignore`, and the `git diff --stat main...HEAD` habit is in `CONTRIBUTING.md`.)
- `git push origin --delete a b` is **atomic**: one stale name fails the whole command and
  takes the other branch down with it. (`gh pr merge --delete-branch` has already removed the
  remote branch anyway.)

### External checks run on the release branch; the tag and `.9000` still precede acceptance
The procedure is in `CONTRIBUTING.md` → "Releasing", step 2. What is not obvious is why it
can work at all.

**The circularity that appears to force checks after the tag is not real.** It looks as
though results cannot be recorded in the commit they describe, because writing them changes
the tree and so invalidates the tarball that was checked. But `cran-comments.md` is
`.Rbuildignore`d and verifiably absent from the built tarball, so editing it produces the
same tarball byte for byte. It is not a package artifact at all — it is the text pasted into
the "Optional comment" field of the CRAN submission form. At 1.2.1 the checks ran *after* the
tag and the `.9000` reopen, so the results landed on a development tree and the tag recorded
no evidence that the tree it names had passed anything.

**Deliberately *not* copied from `usethis:::release_checklist()`: tagging and reopening
`.9000` still happen before CRAN accepts.** The `.9000` suffix is what selects `--quick` over
the full ~20-minute battery, so holding `main` at a clean version across the weeks CRAN can
take would run the full gate on every unrelated PR in that window. And a rejection means
shipping X.Y.Z+1 anyway, leaving the tag naming a tree that was never on CRAN — already true
of `v1.2.0` and `v1.2.2`, so the repository tolerates that cheaply.


### GitHub releases carry notes only — the built tarball is not attached
A release page already offers "Source code (tar.gz)", which GitHub generates from the tag, and
that is **not** the artifact `R CMD build` produces: the git archive has no built vignettes in
`inst/doc/` and does carry every `.Rbuildignore`d development file. Attaching the real tarball
would put two different "source" downloads on one page and a second artifact that can drift from
the tag it claims to be. Follows r-lib and tidyverse practice: CRAN hosts the tarball,
`remotes::install_github('rdotsch/rcicr@vX.Y.Z')` covers everyone else. Being off CRAN is the
strongest counter-argument and was considered; it is temporary.

Accepted cost: the tarball is **not byte-reproducible from the tag**, because `R CMD build`
stamps `Packaged: <timestamp>; <user>`. What has to be reproducible is the *tree*, which the tag
pins exactly and which the gate and `R CMD check` actually operate on.

### Answer CRAN from the review text, and never send them a question a sweep can answer
CRAN's review of 1.2.1 listed two `.Rd` files under "some code lines in examples are
commented out". Working from a summary that had kept one of them, 1.2.2 replied that we could
not find the commented line and **asked the reviewer which line she meant** — while the fix
for the file she had named correctly sat in the same commit, made under a different point.
Two further drafts repeated it.

The logging and re-verification rules that came out of it are in `AGENTS.md`. The one that is
not: **a question to the reviewer costs a round trip measured in weeks, and a sweep costs
minutes.** Sweep everything the point could apply to, then report what was found.

The same instinct applies in reverse — **do not explain a note the reviewer does not have.**
`cran-comments.md` carried a paragraph answering a `medium.com` 403 that appears only on our
local check, on no external check for 1.2.3, and not in CRAN's own pretest of 1.2.1. It was
removed; the link stays in `README.md`, because `README.md` ships and removing it would
invalidate the built tarball and force all five external checks to re-run for something no
CRAN-side check has ever raised (issue #192).

