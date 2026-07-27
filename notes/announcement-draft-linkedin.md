# LinkedIn teaser — draft

*Links to the Medium long-form piece. Keep it short; the caveats live in the article and
in `NEWS.md`, which is where anyone with results in progress should be sent.*

---

Nine years ago I wrote **rcicr**, an R package for reverse correlation image classification — the technique where you show people pairs of noise-masked faces, ask which looks more *trustworthy*, and reconstruct the mental image behind their choices from their responses alone.

In 2021 CRAN archived it. The reason wasn't a bug: email to the maintainer was undeliverable. My old university address had stopped working. The package was fine; the mailbox wasn't. But `install.packages('rcicr')` stopped working for everyone, and the README kept telling people to run it anyway.

I've spent the last stretch fixing that properly, and the hardest part wasn't the code.

Research software gets re-run years later — for a revision, a reproduction, a student picking up an old project. So the rule I set before touching anything was: **never silently change the number a function returns.** A "fix" that quietly alters a published result isn't an improvement.

What that meant in practice:

→ Before fixing anything, I pinned the existing numeric output in a golden-master test. It runs on every change. If a commit moves any of those numbers, it goes red.

→ Seven real correctness bugs fixed — and notably, four of them **failed silently rather than crashing.** The dangerous bugs in research software aren't the ones that error out. Crashes are honest. The dangerous ones hand you a number.

→ Two fixes genuinely *do* change a result — both in the "informational value" metric, both because the old number was computed against the wrong null distribution. Those are documented explicitly, with the exact conditions, rather than buried.

→ 27 dependencies down to 15. A 1.5 GB array no longer copied to every parallel worker. Core computation ~6x faster.

I also used an AI coding assistant throughout, and wrote up honestly where it helped and where it got things wrong — including a benchmark it reported as 29x that was actually 6x, because it had timed only the step it changed rather than the whole function. The mistakes weren't hallucinated APIs. They were confident, plausible, specific and wrong — errors of *verification*, which is exactly the part you can't delegate.

If you have an rcicr analysis in a drawer, the thing to read is `NEWS.md` — specifically the section called "Reproducibility impact."

Full write-up: [link to Medium post]
Package: https://github.com/rdotsch/rcicr

#RStats #OpenScience #ReproducibleResearch #Psychology #Psychophysics
