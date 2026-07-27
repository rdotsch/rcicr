# Bringing rcicr back: modernizing a research package without changing anyone's results

*A note on repairing an R package that people have published with — and on what an AI coding assistant was and wasn't good for.*

---

Nine years ago I wrote [rcicr](https://github.com/rdotsch/rcicr), an R package for reverse correlation image classification — the psychophysics technique where you show people pairs of noise-masked faces, ask which one looks more *trustworthy* (or *criminal*, or *competent*), and reconstruct the mental image driving their choices from their responses alone. People have used it in published work. Then, as happens, I moved on to other things.

In June 2021 CRAN archived it. The stated reason was that email to the maintainer was undeliverable — my old university address had stopped working. The package didn't break; the mailbox did. But the effect was that `install.packages('rcicr')` stopped working for everyone, and the README kept cheerfully telling people to run it.

Last week I finally sat down to fix that. This post is about what I found, what I changed, and — the part that matters most for anyone with an rcicr analysis in a drawer — **what I deliberately did not change**.

## The constraint that shaped everything

Research software has an unusual property: people re-run it years later. A reviewer asks for a revision, a student picks up an old project, someone tries to reproduce a figure. If a "fix" silently changes the number a function returns, you haven't improved the package — you've quietly invalidated a result somebody published.

So I set one rule before touching anything:

> Never change existing call syntax, argument meanings, or the numeric output of a function silently. Deprecate rather than delete. Treat the saved `.Rdata` file as append-only.

Every decision below follows from that.

## What was actually broken

I went through the open issues, the published literature, and the source itself, and reproduced each problem rather than inferring it from reading. Seven correctness bugs came out of that:

- **`generateStimuli2IFC()` didn't save `nscales` or `sigma`** into the `.Rdata` file. That sounds cosmetic. It isn't: the function that builds the null distribution for the "informational value" metric re-generates the stimulus set from that file, and with `nscales` missing it silently fell back to the default of 5. Anyone who used a non-default `nscales` and reported an infoVal got a number computed against the wrong null distribution.
- **`force_gen_ref_dist = TRUE` did nothing.** It was supposed to force recomputation of the reference distribution. It short-circuited one branch but never reached the one that regenerates, so you got the cached distribution and no warning.
- **The `mask` argument was unusable on R ≥ 4.2**, failed a hardcoded size check on anything but 512px stimuli, and rejected greyscale-encoded PNG masks even after successfully converting them.
- **Base images that weren't already exactly `img_size`** failed inside a parallel worker with `non-conformable arrays` — an error that tells you nothing about the actual problem.
- **`generateCI()` broke on tibbles**, which is what you get from any modern data import.
- **`simulateNoiseIntensities()` errored on every call.**
- **The backward-compatibility path for pre-0.3.3 files had never worked** — it was unreachable code that always errored.

There's a pattern in that list worth naming. Four of the seven **failed silently or misleadingly rather than telling the user what went wrong.** The worst bugs in research software aren't crashes. Crashes are honest. The dangerous ones hand you a number.

## How I know your results are safe

This is the part I'd want to read if I had an rcicr analysis in review, so let me be precise rather than reassuring.

**Before fixing anything, I pinned the existing behaviour.** There's now a golden-master test that records the numeric output of the default pipeline — the noise basis, the classification image, all four scaling methods, the infoVal — as it was *before* any fix. It's checked into the repo and runs on every change. If a commit moves any of those numbers, the test goes red and the change doesn't merge without an explicit, documented justification.

It has stayed green through every fix. **That is the evidence for the claim that default-configuration results are unchanged** — not my word, a test anyone can run.

**But "unchanged" isn't the whole story, and I'm not going to pretend it is.** Two fixes genuinely do change a number you might already have:

1. **infoVal, if you used a non-default `nscales`.** The old value was computed against the wrong null distribution. It was wrong. The new one is right. If this applies to you, recompute — and note that old `.Rdata` files don't contain `nscales`, so the fix warns rather than guessing.
2. **infoVal, if you passed `force_gen_ref_dist = TRUE`.** You were getting the cached distribution. Now you get a fresh one, so the value will differ.

Everything else in the fix list only affects code paths that previously *errored out*. If the old code crashed, there was no number to change.

All of this is written up in `NEWS.md` under a section called "Reproducibility impact — read this if you have published or in-progress results," with the exact conditions for each. That file is the thing to read before you re-run anything.

**One more, for completeness.** A performance change (below) made `generateNoiseImage()` compute its average in a different summation order. The results differ by about one unit in the last place — roughly 1e-19 on pixel values of order 0.01. Nothing that reaches a printed number, and at the golden master's own configuration the output is bit-identical. I'm mentioning a difference eight orders of magnitude below anything observable because the standard I set was that *any* numeric change gets written down, rather than discovered later by someone re-running a five-year-old script.

## The other half: making it maintainable

**A test suite where there was none.** 56 tests across 20 files: unit tests for every pure function, targeted tests for the I/O-heavy ones, the golden master, a regression test for each fixed bug, and an end-to-end test that simulates an observer with a known mental template and checks that the pipeline actually recovers it. That last one is my favourite, because nothing else in the suite tested the thing the package exists to do.

**Continuous integration.** `R CMD check` on two R versions on every push and pull request.

**Dependencies: 27 → 15.** None of the twelve removed were used at all. Every unused dependency is a way for someone's install to fail for reasons unrelated to your code.

**A ~1.5 GB array stopped being copied to every parallel worker.** At the default settings, `generateStimuli2IFC()` preallocated the full stimulus array before starting workers, and `foreach` shipped a copy to each one — where each wrote a single slice and discarded the rest. That was the memory ceiling people had been hitting on large stimulus sets.

**`generateNoiseImage()` got about 6× faster**, which matters because it runs once per trial during generation and again for every classification image and z-map.

That last one has a story attached.

## A community contribution, and a lesson about benchmarks

In 2023, [@hvalev](https://github.com/hvalev) opened a pull request replacing an `apply()` call with a vectorized `rowMeans()`, with benchmarks showing a large speedup. I reviewed it, found the numbers didn't reproduce, and we established that the submitted version was wrong: `rowMeans()` on a 3-D array defaults to collapsing two dimensions instead of one, so the result was silently recycled into the wrong shape. No error — just a wrong image.

They posted a corrected version in the thread. And then, because life happens, it sat there for nearly three years.

Coming back to it, I adopted the idea in a slightly different form that avoids materializing an intermediate copy, verified it against an independent implementation of the mathematics — the average written as an explicit triple loop, using neither `apply()` nor `rowMeans()`, and tested on a deliberately *non-square* array so that a transposition couldn't hide — and landed it with credit.

Here's the part worth telling. My first benchmark said **29× faster**. It was wrong. I had precomputed the expensive input *outside* the timed section, so I was measuring only the step I'd changed and not the function anyone actually calls. Measured end to end, it's **6×** — which is exactly what @hvalev had reported in 2023.

Amdahl's law is not optional, and the contributor's original number was better than mine.

## Where AI fit, honestly

I did this work with an AI coding assistant (Claude Code), and since that's now a normal thing to do, it seems worth being specific about what it was and wasn't good for.

**What worked well.** It was genuinely strong at the archaeology: reading nine-year-old code, reproducing bugs from issue reports, tracing why a backward-compatibility branch was unreachable. It wrote the bulk of the test suite. It was good at the kind of care this project needed — when replacing a deprecated random-number call, it checked that the replacement drew from the random stream *identically* rather than just from the same distribution, across 150 seed and probability combinations. The obvious modern replacement would have diverged and silently changed every infoVal computed from a given seed. I would not have thought to check that.

It also found a genuinely subtle bug I'd never have gone looking for: R's `load()` assigns objects straight into the calling function's frame, so an object stored in a saved file can silently overwrite a *function argument* of the same name.

**What didn't.** The 29× benchmark above was its number, and its error. Nothing would have caught it except re-measuring — the code was correct, the speedup was real, the figure was just measuring the wrong thing. It also initially misread the state of that 2023 pull request by reading the diff and not the review conversation underneath it, and told me the contributor's fix was broken when they'd already posted a working version.

Both mistakes have the same shape: **confident, plausible, specific, and wrong in a way that only checking against reality catches.** Neither was a hallucinated API or a syntax error — the failure modes people usually warn about. They were errors of *verification*, which is precisely the thing you cannot delegate.

**So the honest summary:** it made me substantially faster at work I know how to evaluate. Every number in this post exists because something was measured, and the measurements that mattered got redone when they looked too good. The golden-master test isn't there because an AI suggested it; it's there because I wasn't willing to accept "the tests pass" as evidence that published results were safe.

If you take one thing from this section: the value wasn't the code generation. It was having something that would patiently write the test that proves you didn't break anything — the task you always mean to do and never quite do.

## Where it stands

The package is on GitHub, checked on two R versions, and the numbers are pinned. The reason CRAN archived it — an undeliverable maintainer address — no longer applies, and the dependency issues that would have blocked a resubmission are cleared, so returning to CRAN is now a realistic prospect rather than an aspiration.

If you have an rcicr analysis in progress, the thing to read is `NEWS.md`, specifically the "Reproducibility impact" section.

And if you hit a bug: the issue tracker is being read again.

---

*rcicr: https://github.com/rdotsch/rcicr*
*The original introduction to the method: [Reverse correlation image classification using R](https://medium.com/@rondotsch/reverse-correlation-image-classification-using-r-a0701648fb0/)*
