# Releasing rcicr

Maintainer-facing; contributors want `CONTRIBUTING.md`. Releases are
squash-merged onto `main` and marked with a tag — no `develop` branch.
`DESCRIPTION` carries a `.9000` suffix between releases and the release
commit drops it.

**Keep this file under 1600 words.**

------------------------------------------------------------------------

## 0. The reproducibility gate

**A release does not go out while this is red.** Both runs are required;
they answer different questions.

``` sh
Rscript tools/compare-release-output.R               # vs v1.0.1 -- the published baseline
Rscript tools/compare-release-output.R --ref=v1.1.0  # vs the previous release
```

The v1.0.1 run asks whether this tree still produces the numbers that
are *in the literature*. It never advances, because a reference that
moves forward with each release lets a tree drift away from those
numbers one tolerated epsilon at a time. The second asks whether
anything broke since the last release, and reaches further: calls that
used to crash — masks, small z-maps, undecorated z-maps — have a
comparable value in v1.1.0 and none at all in v1.0.1.

Exit `0` clean, `1` an unaccounted difference, `2` could not run. On a
difference there are two honest outcomes: **it was not intended**, so
fix it; or **it was intended**, so write it up in `NEWS.md` under
“Reproducibility impact” — who is affected, what to do — and add an
`EXPECTED` entry naming the reference, the output, the reason and that
heading. The script fails if an `EXPECTED` entry stops firing (the list
cannot rot) and if one fires while `NEWS.md` is silent (a deviation
cannot be quiet).

Widening a tolerance, dropping a configuration, or skipping the run is
none of those.

## 1. Prepare the release PR

Branch from `main` (`release-X.Y.Z`) and make five edits together:

- **`NEWS.md`** — rename `# rcicr (development version)` to
  `# rcicr X.Y.Z (YYYY-MM-DD)`. Until you do, none of this release’s
  entries are in the news database: R indexes only sections under a
  heading it can read a version from. Keep version numbers *out* of `##`
  headings or the file stops parsing and `R CMD check` NOTEs.
- **`DESCRIPTION`** — drop the `.9000`. This also switches CI from the
  `--quick` gate to the full battery, so it must happen **before** the
  PR is opened.
- **`ChangeLog`** — a dated pointer entry deferring to `NEWS.md`, never
  a duplicate.
- **`CITATION.cff`** — regenerate *after* the `DESCRIPTION` edit, with
  `cffr::cff_write("DESCRIPTION", dependencies = FALSE, gh_keywords = FALSE)`.
  CI compares it and fails until you do.
- **`cran-comments.md`** — re-read every claim rather than trusting it.
  It has gone stale *within a single day*, describing a URL reference
  another PR had just removed, one step short of a CRAN reviewer reading
  it. Test-environment results come in step 2.

Open the PR; the full battery runs on it (~20 min) and is a required
check.

## 2. External checks, on the release branch, before merging

win-builder and R-hub run against the release branch and their results
go into `cran-comments.md` **in the same PR**, so the release commit
carries its own evidence.

``` sh
R CMD build .                                # at the repo root, never in a worktree
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-devel/
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-release/
```

Then trigger R-hub from the Actions tab against the release branch. Two
NOTEs are expected locally and explained in `cran-comments.md`: CRAN
incoming feasibility (new submission of an archived package) and future
file timestamps (no network clock in a sandbox).

**Why results can be written into the commit they describe:**
`cran-comments.md` is `.Rbuildignore`d and verifiably absent from the
tarball, so editing it produces the same tarball byte for byte. It is
not a package artifact at all — it is the text pasted into the
submission form. At 1.2.1 these checks ran *after* the tag and the
`.9000` reopen, so the results landed on a development tree and the tag
recorded no evidence that the tree it names had passed anything.

Neither check is run casually: win-builder mails the maintainer, and
both put the tarball in front of a third party.

- **win-builder** — `226` means the transfer completed. The server
  **will not overwrite an existing upload**: a second attempt at the
  same filename fails with a bare `550`, meaning “already queued”, not
  “rejected”. Confirm with
  `curl --list-only ftp://win-builder.r-project.org/R-devel/` before
  re-uploading. The mail carries a result URL;
  `curl -s https://win-builder.r-project.org/<key>/00check.log` fetches
  the full log. Pages expire after ~72 hours, so `cran-comments.md` is
  the durable copy.
- **R-hub** — the stock v2 workflow, run from this repository’s Actions
  tab. Because it is `workflow_dispatch`-only it must reach the
  **default branch** before it can be triggered at all. **Download the
  artifacts; do not use `gh run view --log`**, which truncates: at 1.2.1
  the 38-minute macOS job returned 1548 lines ending four minutes in,
  with no `Status:` line. `gh run download <run-id> -D <dir>`, then read
  `*/rcicr.Rcheck/00check.log`. Expect `Status: OK` where win-builder
  reports 1 NOTE — R-hub does not run the incoming feasibility check, so
  that is agreement, not a discrepancy.
- **`RHUB_TOKEN` is deliberately unset and nothing is missing.** The
  stock workflow passes it to four actions, which makes it look
  required. It is an optional slot for your own PAT to reach *private*
  repositories, not a credential R-hub issues, and no step in
  `r-hub/actions@v1` references it. This package is public; an unset
  secret expands to empty and the jobs run.

## 3. Merge and tag

``` sh
gh pr merge <n> --squash --delete-branch
git checkout main && git pull && git fetch --prune
git tag -a vX.Y.Z -m "rcicr X.Y.Z" && git push origin vX.Y.Z
gh release create vX.Y.Z --title "rcicr X.Y.Z" --notes-file <(...)   # NEWS section
git diff <pr-head-sha> vX.Y.Z --stat                                 # must be empty
```

The tag push re-runs the full gate against the released tree. The empty
diff is what lets step 2’s results stand for the tagged tree.

**The release page carries notes only — the tarball is not attached.**
GitHub already offers “Source code (tar.gz)” generated from the tag, and
that is *not* what `R CMD build` produces: no built vignettes in
`inst/doc/`, and every `.Rbuildignore`d development file present. Two
different “source” downloads on one page is a support question waiting
to happen. CRAN hosts the real tarball and
`remotes::install_github('rdotsch/rcicr@vX.Y.Z')` covers everyone else.
Accepted cost: the tarball is not byte-reproducible from the tag,
because `R CMD build` stamps `Packaged: <timestamp>; <user>`. What must
be reproducible is the *tree*, which the tag pins.

## 4. Submit to CRAN

``` sh
git switch --detach vX.Y.Z      # never from main HEAD, never in a worktree
R CMD build .
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
git switch -
```

Building from the tag keeps the development suffix out of the submitted
version: `Version contains large components` only blocks if the
*tarball* carries it.

Submit at <https://cran.r-project.org/submit.html>, pasting the body of
`cran-comments.md` into the “Optional comment” field — `.Rbuildignore`d,
so it reaches CRAN only this way. **Paste the tag’s copy**
(`git show vX.Y.Z:cran-comments.md`): submission can trail tagging by
weeks while `main` moves, and `main`’s copy would describe check results
from a tree that is not the tarball.

**Ron submits personally.** CRAN emails the maintainer address to
confirm, and for a package archived over an undeliverable address, that
confirmation *is* the point of the resubmission.

### Answering a review

**Never send CRAN a question a sweep can answer.** A round trip costs
weeks; a sweep costs minutes. The review of 1.2.1 named two `.Rd` files
with commented-out example lines; working from a summary that had kept
only one, 1.2.2 asked the reviewer which line she meant — while the fix
for the file she had named sat in the same commit under a different
point. Two further drafts repeated it. Sweep everything the point could
apply to, then report what was found.

**And do not explain a note the reviewer does not have.**
`cran-comments.md` once answered a `medium.com` 403 that appears only on
our local check — not on any external check, not in CRAN’s own pretest.
Removed; the link stays in `README.md` (issue \#192).

## 5. Reopen development

Bump `DESCRIPTION` to `X.Y.Z.9000` and start a fresh
`# rcicr (development version)` heading in `NEWS.md`, through a small PR
like anything else.

**Tagging and reopening `.9000` deliberately precede CRAN acceptance**,
unlike `usethis:::release_checklist()`, which tags only once CRAN has
accepted. Two reasons. The `.9000` suffix is what selects `--quick` over
the full ~20-minute battery, so holding `main` at a clean version for
the weeks CRAN can take would run the full gate on every unrelated PR in
that window. And **a tag here marks a release, not an acceptance**:
GitHub is a real distribution channel for this package — 1.0.1 through
1.2.3 were released there and nowhere else — and while a submission is
under review it is the only way anyone can install the new version.
Tagging on acceptance would point
`remotes::install_github('rdotsch/rcicr@*release')` at the previous
release for as long as CRAN takes.

So a tag naming a tree CRAN never took is expected here rather than a
defect — already true of `v1.2.0` and `v1.2.2` — and answering a review
means shipping X.Y.Z+1, an ordinary release whose `cran-comments.md`
happens to be a point-by-point reply. **Which versions CRAN accepted is
recorded rather than inferred from the tags**: in
`notes/cran-review-<version>.md` and in the tag’s GitHub release notes.

The alternative is to stop making GitHub releases once the package is
back on CRAN and adopt the `usethis` order wholesale. That trade is only
worth revisiting if CRAN becomes the channel people actually install
from.
