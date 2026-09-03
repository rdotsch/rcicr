# CRAN review of the 1.3.0 submission

**Pretest received** 2026-08-22, from the CRAN teams' auto-check service. This was the automated step before a team member's manual inspection.

**Submission** `rcicr_1.3.0.tar.gz`, the reinstatement submission after the 2021-06-08
archival, carrying the answers to Konstanze Lauseker's review of 1.2.1
(`notes/cran-review-1.2.1.md`). The point-by-point response is the body of
`cran-comments.md` as it stood at the `v1.3.0` tag.

**Outcome: accepted.** CRAN published rcicr 1.3.0 on 2026-09-02. Leonore Hochhauser's confirmation is quoted below. The preceding pretests had one NOTE on each of the two flavours: the incoming feasibility NOTE expected for a reinstatement. Nothing in it was a defect, and `cran-comments.md` already explained all four parts to the reviewer.

**Why this file exists before there is a review to log.** The `incoming_pretest` log
directory is win-builder scratch space on its own cleanup schedule, and it went early: the
mail said "roughly 7 days", but by **2026-08-27**, day five, none of the three URLs below
resolved for the maintainer any more — not the log dir, and not either of the two
`00check.log` files under it, which can outlive an unbrowsable directory but here did not.
The full logs are therefore gone unless someone saved them, and the NOTE text survives only
because the mail reproduces it — which is the whole of what was ever actionable.

**What survives, and what does not.** The NOTE itself survives because the mail quotes it
in full, and `cran-comments.md` holds our own pre-submission evidence: the same single NOTE
from win-builder on release and devel, and `Status: OK` from R-hub on Linux, macOS and
Windows. The pretest run itself is unrecoverable — it was CRAN's own incoming check under
R Under development r90440 of 2026-08-21, on a Debian and a Windows flavour. A later
`devtools::check_win_devel()` sends the same tree to whatever R-devel, dependency versions
and check infrastructure exist at that point, and win-builder has no Debian flavour at all.
It can show that the `v1.3.0` tree still passes; it cannot reconstitute this evidence, and a
fresh log should never be filed as if it were this one.

The directory's removal is a cleanup job, not a status signal: it says nothing about where
the submission sits in CRAN's queue.

**Two observations that are easy to misread, recorded so they are not re-derived:**

- On 2026-08-27 `rcicr` appeared in none of the folders under
  <https://cran.r-project.org/incoming/>. That is not evidence of a problem. Absence there
  means the tarball is not sitting in a publicly listed working directory at that moment;
  the auto-check mail is proof it was accepted into the pipeline, and an item awaiting
  manual inspection need not be visible in that listing.
- On 2026-08-27 <https://cran.r-project.org/package=rcicr> did not yet resolve to a current package, and the check results still showed the 2021 archival. Both changed when CRAN published 1.3.0 on 2026-09-02.

---

## The acceptance mail, verbatim

Thanks,
on its way to CRAN.

Best,
Leonore Hochhauser (she/her)

---

## The auto-check mail, verbatim

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

No reply to the auto-check was needed. The acceptance mail above closed the manual review, and the package's CRAN page records the published outcome.
