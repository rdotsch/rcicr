# Releasing rcicr

For the maintainer; contributors want `CONTRIBUTING.md`. Releases are
squash-merged onto `main` and marked with a tag; there is no `develop`
branch. `DESCRIPTION` carries a `.9000` suffix between releases, and the
release commit drops it.

**Keep this file under 1600 words.**

------------------------------------------------------------------------

## 0. The reproducibility gate

**A release does not go out while this is red.** Both runs are required,
because they answer different questions.

``` sh
Rscript tools/compare-release-output.R  # vs v1.0.1 -- the published baseline
Rscript tools/compare-release-output.R --ref="$(git tag --list 'v*' --sort=-v:refname | head -1)"
```

The second reference is derived rather than named, so it always is the
previous release; CI picks it the same way.

The v1.0.1 run asks whether this tree still produces the numbers that
are *in the literature*. That reference never advances: one that moved
with each release would let the tree drift away from those numbers one
tolerated epsilon at a time. The second run asks whether anything broke
since the last release. It also reaches further, because calls that used
to crash (masks, small z-maps, undecorated z-maps) have a value to
compare against from v1.1.0 on, and none in v1.0.1.

The script exits `0` when clean, `1` on an unexplained difference, and
`2` if it could not run. A difference has two honest outcomes. **It was
not intended:** fix it. **It was intended:** describe it in `NEWS.md`
under “Reproducibility impact” (who is affected, what they should do)
and add an `EXPECTED` entry naming the reference, the output, the reason
and that heading. The script fails when an `EXPECTED` entry stops
firing, so the list cannot rot, and when one fires while `NEWS.md` says
nothing, so a deviation cannot go unannounced.

Widening a tolerance, dropping a configuration or skipping the run is
neither outcome.

## 1. Prepare the release PR

Branch from `main` (`release-X.Y.Z`) and make four edits together:

- **`NEWS.md`**: rename `# rcicr (development version)` to
  `# rcicr X.Y.Z (YYYY-MM-DD)`. R indexes only sections under a heading
  it can read a version from, so until then none of this release’s
  entries are in the news database. Keep version numbers *out* of `##`
  headings, or the file stops parsing and `R CMD check` NOTEs.
- **`DESCRIPTION`**: drop the `.9000`. This also switches CI from the
  `--quick` gate to the full battery, so do it **before** opening the
  PR.
- **`CITATION.cff`**: regenerate it *after* the `DESCRIPTION` edit, with
  `cffr::cff_write("DESCRIPTION", dependencies = FALSE, gh_keywords = FALSE)`.
  CI compares it and fails until you do.
- **`cran-comments.md`**: re-check every claim. It has gone stale within
  a single day, describing a URL another PR had just removed.
  Test-environment results come in step 2.

Open the PR. The full battery (about 20 minutes) runs on it as a
required check.

## 2. External checks, on the release branch, before merging

Run win-builder and R-hub against the release branch and put their
results in `cran-comments.md` **in the same PR**, so the release commit
carries its own evidence.

``` sh
R CMD build .                                # at the repo root, never in a worktree
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-devel/
curl -T rcicr_X.Y.Z.tar.gz ftp://win-builder.r-project.org/R-release/
```

Then trigger R-hub from the Actions tab against the release branch.
Record each check’s actual errors, warnings and NOTEs, with its R
version and platform. Do not carry the archived-package reinstatement
NOTE into a routine update, and report sandbox-only NOTEs as local
results.

**Why the results can go into the commit they describe:**
`cran-comments.md` is `.Rbuildignore`d and absent from the tarball, so
editing it changes no package source. It is not a package file at all,
but the text pasted into the submission form. (Two fresh builds are not
byte-identical anyway: packaging timestamps and metadata differ.)
Running these checks after tagging instead would leave the tag with no
evidence that the tree it names passed anything.

Neither check is run casually: win-builder emails the maintainer, and
both put the tarball in front of a third party.

- **win-builder**: `226` means the transfer completed. The server **will
  not overwrite an existing upload**, so a second attempt at the same
  filename fails with a bare `550`, meaning “already queued”, not
  “rejected”. Confirm with
  `curl --list-only ftp://win-builder.r-project.org/R-devel/` before
  uploading again. The email carries a result URL;
  `curl -s https://win-builder.r-project.org/<key>/00check.log` fetches
  the full log. Pages expire after about 72 hours, so `cran-comments.md`
  is the lasting copy.
- **R-hub**: the stock v2 workflow, run from this repository’s Actions
  tab. It is `workflow_dispatch`-only, so it must be on the **default
  branch** before it can be triggered at all. **Download the artifacts;
  do not use `gh run view --log`**, which truncates: for a 38-minute
  macOS job it returned a log ending four minutes in, with no `Status:`
  line. Run `gh run download <run-id> -D <dir>`, then read
  `*/rcicr.Rcheck/00check.log`. R-hub skips the incoming feasibility
  checks and uses `--no-manual`, so a green run establishes neither; get
  a complete manual check before submitting.
- **`RHUB_TOKEN` is unset on purpose, and nothing is missing.** The
  stock workflow passes it to four actions, which makes it look
  required. It is an optional slot for your own token to reach *private*
  repositories, not a credential R-hub issues, and no step in
  `r-hub/actions@v1` uses it. This package is public; an unset secret
  expands to empty and the jobs run.

## 3. Merge and tag

``` sh
gh pr merge <n> --squash --delete-branch
git checkout main && git pull && git fetch --prune
git tag -a vX.Y.Z -m "rcicr X.Y.Z" && git push origin vX.Y.Z
gh release create vX.Y.Z --title "rcicr X.Y.Z" --notes-file <(...)   # NEWS section
git diff <pr-head-sha> vX.Y.Z --stat                                 # must be empty
```

The tag push re-runs the full gate on the released tree. The empty diff
is what lets step 2’s results stand for the tagged tree.

**The release page carries notes only; do not attach the tarball.**
GitHub already offers “Source code (tar.gz)” generated from the tag, and
that is *not* what `R CMD build` produces: it lacks the built vignettes
in `inst/doc/` and includes every `.Rbuildignore`d development file. Two
different “source” downloads on one page invite support questions. CRAN
hosts the real tarball, and
`remotes::install_github('rdotsch/rcicr@vX.Y.Z')` covers everyone else.
The accepted cost: `R CMD build` stamps `Packaged: <timestamp>; <user>`,
so the tarball is not byte-reproducible from the tag. What must be
reproducible is the *tree*, and the tag pins it.

## 4. Submit to CRAN

``` sh
git switch --detach vX.Y.Z      # never from main HEAD, never in a worktree
R CMD build .
R CMD check --as-cran rcicr_X.Y.Z.tar.gz
git switch -
```

Building from the tag keeps the development suffix out of the submitted
version: `Version contains large components` only blocks when the
*tarball* carries it.

Submit at <https://cran.r-project.org/submit.html>, pasting the body of
`cran-comments.md` into the “Optional comment” field. The file is
`.Rbuildignore`d, so this is the only way it reaches CRAN. **Paste the
tag’s copy** (`git show vX.Y.Z:cran-comments.md`): submission can trail
tagging by weeks while `main` moves on, and `main`’s copy would then
describe checks on a tree that is not the tarball.

**An agent may upload; only Ron confirms.** Nothing reaches CRAN until
the maintainer clicks the link CRAN emails him.

### Answering a review

**Never send CRAN a question a sweep can answer.** A round trip costs
weeks, a sweep minutes. The review of 1.2.1 named two `.Rd` files with
commented-out example lines. Working from a summary that kept only one,
the reply asked the reviewer which line she meant, while the fix for the
file she had named sat in the same commit. Sweep everything the point
could apply to, then report what you found.

**Do not explain a NOTE the reviewer did not get.** A local-only 403 is
a sandbox result, not something to raise in `cran-comments.md`.

## 5. Reopen development

Bump `DESCRIPTION` to `X.Y.Z.9000` and start a fresh
`# rcicr (development version)` heading in `NEWS.md`, in a small PR like
any other.

**Tagging and reopening at `.9000` come before CRAN acceptance**, unlike
`usethis:::release_checklist()`, which tags only after acceptance. There
are two reasons:

- The `.9000` suffix is what selects the `--quick` gate. Holding `main`
  at a clean version for the weeks CRAN can take would run the full
  20-minute battery on every unrelated PR in that window.
- **A tag here marks a release, not an acceptance.** GitHub is a real
  distribution channel for this package: 1.0.1 through 1.2.3 were
  released there and nowhere else. While a submission is under review,
  GitHub is the only way to install the new version. Tagging on
  acceptance would point
  `remotes::install_github('rdotsch/rcicr@*release')` at the previous
  release for as long as CRAN takes.

So a tag naming a tree CRAN never accepted is expected here, not a
defect; `v1.2.0` and `v1.2.2` are examples. **Which versions CRAN
accepted is recorded, not inferred from the tags**: in
`notes/cran-review-<version>.md` and in the tag’s GitHub release notes.

The `usethis` order was reconsidered once the package was back on CRAN,
and rejected. It would delay GitHub releases until acceptance and drop
them entirely for a declined version. Its one gain, never tagging a
declined tree, does not survive CRAN’s rule that a resubmission carries
a new version number: answering a review means shipping X.Y.Z+1 anyway,
and the earlier tag marks a tree that was built and sent.
