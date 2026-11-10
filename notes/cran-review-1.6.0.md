# CRAN reply to the 1.6.0 submission

**Received** 2026-11-10 from Uwe Ligges, answering the maintainer's reply-all
to the auto-check (`notes/cran-review-1.6.0-pretest.md`).

**Outcome:** 1.6.0 cancelled at the maintainer's request; 1.6.1 follows with
the MKL-sensitive test comparisons made to a tolerance.

---

## The thread, verbatim

Thanks for letting us know, submission cancelled.

Best,
Uwe Ligges


On 10.11.2026 11:45, Ron Dotsch wrote:
Dear CRAN team,
Thank you. rcicr 1.6.0 will not pass the MKL additional check either:
the tests that fail for 1.5.0 under MKL are still present in 1.6.0.
They compare reference norms computed twice in the same session with
expect_identical(), and MKL's results differ between runs in the last
bit (for example 62.502882847115842 against 62.502882847115885). The
package's results are unaffected; the tests are too strict.
I will submit 1.6.1 with those comparisons made to a numerical
tolerance, where bit-identical results are not guaranteed. Please feel
free to cancel 1.6.0.
Best regards,
Ron Dotsch
Op di 10 nov 2026 om 11:29 schreef <ligges@statistik.tu-dortmund.de>

The quoted auto-check mail that follows is identical to the one in
`notes/cran-review-1.6.0-pretest.md`.
