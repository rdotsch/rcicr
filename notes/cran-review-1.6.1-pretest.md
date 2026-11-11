# CRAN auto-check of the 1.6.1 submission

**Submitted** 2026-11-10: tarball built from the `v1.6.1` tag, uploaded through CRAN's submission form and confirmed by the maintainer by email; CRAN's "package was uploaded" confirmation followed. The submitter's comment was the tag's `cran-comments.md`.

**Pretest** run 2026-11-10 (log directory `rcicr_1.6.1_20261110_193509`) by the CRAN teams' auto-check service.

**Outcome:** Windows and Debian pretests passed with Status: OK. The auto-check again flagged the MKL additional issue of the *last released* version (1.5.0), so CRAN will check 1.6.1 manually. 1.6.1 is the version that addresses that issue (see its `cran-comments.md`), and the mail asks for a reply only if the submitted version will not pass, so it was not answered.

---

## The mail, verbatim

Dear maintainer,

package rcicr_1.6.1.tar.gz has been auto-processed.
The auto-check found additional issue(s) for the *last* version released on CRAN:
  MKL <https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/rcicr.out>
CRAN will have to check the additional issue(s) manually.
IMMEDIATELY reply-all in case you submitted a version that will not pass these checks.
Please do not respond if you expect the additional checks will pass.

Log dir: <https://win-builder.r-project.org/incoming_pretest/rcicr_1.6.1_20261110_193509/>
The files will be removed after roughly 7 days.
Installation time in seconds: 9
Check time in seconds: 208
R Under development (unstable) (2026-11-07 r90850 ucrt)

Pretests results:
Windows: <https://win-builder.r-project.org/incoming_pretest/rcicr_1.6.1_20261110_193509/Windows/00check.log>
Status: OK
Debian: <https://win-builder.r-project.org/incoming_pretest/rcicr_1.6.1_20261110_193509/Debian/00check.log>
Status: OK

Last released version's CRAN status: OK: 13
See: <https://CRAN.R-project.org/web/checks/check_results_rcicr.html>

Last released version's additional issues:
  MKL <https://www.stats.ox.ac.uk/pub/bdr/Rblas/MKL/rcicr.out>

CRAN Web: <https://cran.r-project.org/package=rcicr>

No strong reverse dependencies to be checked.

Best regards,
CRAN teams' auto-check service
Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
Check: *, Result: OK
