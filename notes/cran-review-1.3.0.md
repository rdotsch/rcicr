# CRAN auto-check of the 1.3.0 submission

**Received** 2026-08-22, from the CRAN teams' auto-check service. Automated, not a team
member's review — the manual inspection it announces had not arrived as of 2026-08-27.

**Submission** `rcicr_1.3.0.tar.gz`, the reinstatement submission after the 2021-06-08
archival, carrying the answers to Konstanze Lauseker's review of 1.2.1
(`notes/cran-review-1.2.1.md`). The point-by-point response is the body of
`cran-comments.md` as it stood at the `v1.3.0` tag.

**Outcome: pretests passed.** One NOTE on each of the two flavours, and it is the incoming
feasibility NOTE that a reinstatement is expected to produce. Nothing in it is a defect and
nothing in it needs a fix — `cran-comments.md` already explains all four parts to the
reviewer. What remains is the manual inspection, quoted below as "within the next 10 working
days", which from 2026-08-22 lands around **2026-09-04**.

**Why this file exists before there is a review to log.** The `incoming_pretest` log
directory is deleted after roughly seven days, so the two `00check.log` URLs below stop
resolving around 2026-08-29. The NOTE text that mattered is reproduced in the mail itself
and therefore survives here; the full logs do not, unless someone saved them.

**Two observations that are easy to misread, recorded so they are not re-derived:**

- On 2026-08-27 `rcicr` appeared in none of the folders under
  <https://cran.r-project.org/incoming/>. That is not evidence of a problem. Absence there
  means the tarball is not sitting in a publicly listed working directory at that moment;
  the auto-check mail is proof it was accepted into the pipeline, and an item awaiting
  manual inspection need not be visible in that listing.
- <https://cran.r-project.org/package=rcicr> not resolving, and the check results still
  showing the 2021 archival, are both expected until a team member publishes the package.

---

## The mail, verbatim

Dear maintainer,

package rcicr_1.3.0.tar.gz has been auto-processed and is pending a manual inspection of this new CRAN submission. A CRAN team member will typically respond to you within the next 10 working days. For technical reasons you may receive a second copy of this message when a team member triggers a new check.

Log dir: <https://win-builder.r-project.org/incoming_pretest/rcicr_1.3.0_20260822_194825/>
The files will be removed after roughly 7 days.
Installation time in seconds: 9
Check time in seconds: 129
R Under development (unstable) (2026-08-21 r90440 ucrt)

Pretests results:
Windows: <https://win-builder.r-project.org/incoming_pretest/rcicr_1.3.0_20260822_194825/Windows/00check.log>
Status: 1 NOTE
Debian: <https://win-builder.r-project.org/incoming_pretest/rcicr_1.3.0_20260822_194825/Debian/00check.log>
Status: 1 NOTE

Last released version's CRAN status: OK: 13
See: <https://cran-archive.r-project.org/web/checks/2021/2021-06-08_check_results_rcicr.html>

No strong reverse dependencies to be checked.

Best regards,
CRAN teams' auto-check service
Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
Check: CRAN incoming feasibility, Result: NOTE
  Maintainer: 'Ron Dotsch <rdotsch@gmail.com>'

  New submission


  Package was archived on CRAN

  Possibly misspelled words in DESCRIPTION:
    Brinkman (13:5)
    Dotsch (11:68, 13:27)
    Todorov (12:5, 13:15)
    psychophysical (10:33)

  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2021-06-08 as email to the maintainer was
      undeliverable.

---

## What the NOTE contains, and why none of it is actionable

Each part was anticipated in `cran-comments.md` before submission:

- **`New submission` and `Package was archived on CRAN`, with the db override** — both are
  what a reinstatement produces by definition. The override records the 2021-06-08
  archival, whose stated reason was an undeliverable maintainer address; `Authors@R` now
  carries `rdotsch@gmail.com`, which is the address this mail was delivered to.
- **Four possibly misspelled words** — `Brinkman`, `Dotsch` and `Todorov` are author
  surnames from the references CRAN asked for in the 1.2.1 review, and `psychophysical` is
  the standard adjective for the field the method comes from.

The pretest agrees with what we ran before submitting: win-builder release and devel each
reported the same single NOTE, and R-hub reported `Status: OK` on Linux, macOS and Windows.
R-hub does not run the incoming feasibility check, so that is agreement rather than a
discrepancy.

**Do not chase this.** `RELEASING.md` → "Answering a review" is written around the cost of a
round trip; there is nothing here to answer, and a status enquiry inside the stated window
would be one. The next action is triggered by a team member's mail, or by that window
passing without one.
