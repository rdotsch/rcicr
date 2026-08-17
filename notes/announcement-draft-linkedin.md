# LinkedIn teaser — draft

*Links to the Medium long-form piece. Keep it short; the caveats live in the article and
in `NEWS.md`, which is where anyone with results in progress should be sent.*

---

Twelve years ago I wrote **rcicr**, an R package for reverse correlation image classification — the technique where you show people pairs of noise-masked faces, ask which looks more *trustworthy*, and reconstruct the mental image behind their choices from their responses alone.

Then I left academia. That was 2017. I hoped someone would take the package up, and in a limited sense people did — a few good fixes arrived over the years, and I'm grateful for each. But a contributor isn't a maintainer, and nobody ever took it on. There's no mechanism for handing this kind of thing over: you publish a tool, people cite it, and then it's yours forever or it's nobody's.

In 2021 CRAN archived it. The reason wasn't a bug: email to the maintainer was undeliverable — a university address I'd stopped having four years earlier. The package was fine; the mailbox wasn't. But `install.packages('rcicr')` stopped working for everyone, and the README kept telling people to run it anyway.

I've finally fixed that properly. Not because I went back — I didn't. Because with an AI assistant the work fit into evenings, and modernizing a twelve-year-old package is mostly unglamorous archaeology: reading your own old code, reproducing bug reports, working out why a branch was never reachable. That's the part that needs a maintainer's full attention, which is exactly what I no longer have.

The hardest part still wasn't the code.

Research software gets re-run years later — for a revision, a reproduction, a student picking up an old project. So the rule I set before touching anything was: **never silently change the number a function returns.** A "fix" that quietly alters a published result isn't an improvement.

What that meant in practice:

→ Before fixing anything, I pinned the existing numeric output in a golden-master test. It runs on every change. If a commit moves any of those numbers, it goes red.

→ Seven real correctness bugs fixed — and notably, four of them **failed silently rather than crashing.** The dangerous bugs in research software aren't the ones that error out. Crashes are honest. The dangerous ones hand you a number.

→ Two fixes genuinely *do* change a result — both in the "informational value" metric, both because the old number was computed against the wrong null distribution. Those are documented explicitly, with the exact conditions, rather than buried.

→ 27 dependencies down to 15. A 1.5 GB array no longer copied to every parallel worker. Core computation ~6x faster.

→ And one correction I'd rather not be writing: for nine years, one optional way of saving per-participant images wrote them under the **wrong participant's filename**. The images were right; the names weren't. It never reached CRAN and never touched the documented route for per-participant images — but if you used that option, the article shows how to check in one line and fix it by renaming rather than re-running.

The timing of that last one is the part I keep thinking about. The feature was a contribution from someone else, and I merged it on 15 August 2017 — which makes it mine; reviewing is the job. My last real commits are from that September. It landed as I was walking out, and then there was nobody there to read it — for nine years. **The AI found it**, not by hunting for bugs but by starting on an unrelated change, reading the line that picks each participant's trials and the line that names their file, and noticing they disagree. An unmaintained package doesn't usually fail from a lack of judgement. It fails from a lack of anyone reading it.

If you have an rcicr analysis in a drawer, the thing to read is `NEWS.md` — specifically the section called "Reproducibility impact."

Full write-up: [link to Medium post]
Package: https://github.com/rdotsch/rcicr

#RStats #OpenScience #ReproducibleResearch #Psychology #Psychophysics
