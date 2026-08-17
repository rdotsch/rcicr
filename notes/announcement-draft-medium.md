# Bringing rcicr back: modernizing a research package without changing anyone's results

*A note on repairing an R package that people have published with — and on what an AI coding assistant was and wasn't good for.*

---

Twelve years ago I wrote [rcicr](https://github.com/rdotsch/rcicr), an R package for reverse correlation image classification — the psychophysics technique where you show people pairs of noise-masked faces, ask which one looks more *trustworthy* (or *criminal*, or *competent*), and reconstruct the mental image driving their choices from their responses alone. The first commit is dated 28 July 2014 and its message reads, in full, "V0.1 first commit". It was on CRAN three days later, and people have published with it since.

Then I left academia. That was 2017; my last real commits are from that September.

I did what I think most people in that position do: I hoped someone would take it up. In a limited sense, people did — a handful of good fixes arrived over the years, one in 2018, one in 2019, a couple in 2020, one in 2023. I'm grateful for every one. But a contributor is not a maintainer, and nobody ever took the package on. There is no mechanism in academic research software for handing something over. You publish the tool, people cite it, and then it is yours forever or it is nobody's.

So it sat. In June 2021 CRAN archived it, and the stated reason wasn't a bug in the code — email to the maintainer was undeliverable. It was a university address I'd stopped having four years earlier. The package didn't break; the mailbox did. But the effect was that `install.packages('rcicr')` stopped working for everyone, and the README kept cheerfully telling people to run it.

What changed this year isn't that I went back — I didn't, and I'm not going to. What changed is that the work became possible in the margins. Modernizing a twelve-year-old package is weeks of unglamorous archaeology: reading your own old code, reproducing bug reports, working out why a backward-compatibility branch was never reachable. That is precisely the kind of work an AI assistant is good at, and with one I could do in evenings what would otherwise have needed a maintainer's full attention — which is the thing I no longer have to give.

This post is about what I found, what I changed, and — the part that matters most for anyone with an rcicr analysis in a drawer — **what I deliberately did not change**. It also includes a correction I'd rather not be writing, which the same assistant is the reason I can write at all.

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

## A correction: some individual-CI images were saved under the wrong participant's name

This one is not a rounding difference, and it deserves its own heading rather than a line in a fix list.

If you called `generateCI()` with a `participants` vector **and** `save_individual_cis = TRUE`, the per-participant images written to `individual_cis/` could be saved under the wrong participant's filename. The loop picked each participant's trials in *sorted* order but took the filename from order of *appearance*. When those two orders differ — any data frame not sorted by participant, `"p2"` filed before `"p10"` because sorting is lexical, IDs collected out of order — every file in that folder got somebody else's ID.

The images themselves were always computed correctly. It is purely a naming error. But a naming error on a per-participant figure is not a small thing: if you published one of those images as participant `p2`, it may be a different participant's classification image.

**Here is who this does not touch, and it is most people.** The documented way to compute per-participant classification images has always been the batch functions — `batchGenerateCI2IFC()` in the old CRAN documentation and in the official demo script, `batchGenerateCI()` in the current vignette. None of those, and not `generateCI2IFC()` either, even exposes the two arguments involved. They split the data themselves and name each image from the grouping value directly, so they never reach the code that was wrong. Nothing `generateCI()` *returns* was affected either — the group classification image is an average across participants, so it doesn't depend on their order, and no z-map or stored number ever changed.

**And it never reached CRAN.** rcicr was on CRAN from July 2014 until it was archived in June 2021, and the last release there, 0.3.4.1 from July 2016, predates the `save_individual_cis` option by thirteen months. Nobody who ran `install.packages('rcicr')` ever got a mislabelled file. The defect entered the development branch on 15 August 2017 — in, of all things, the commit that first made saving individual images work at all — and every GitHub release since carries it.

**You can check your own output in one line, without re-running anything:**

```r
identical(unique(participants), sort(unique(participants)))   # TRUE means your files were fine
```

And if it comes back `FALSE`, the fix is a rename rather than a recomputation, because the mapping is exact — the file named `unique(participants)[i]` holds the image of `sort(unique(participants))[i]`.

The timeline is uncomfortably exact, and I'd rather set it out than have someone find it in the log.

The feature that saves per-participant images was contributed by someone else, as a pull request, and merged on **15 August 2017**. I merged it, which makes it mine — reviewing is the job, and I did not catch this. My last substantive commits to the package are from **25 September 2017**, about six weeks later. So it went in essentially as I was walking out of the building, and then there was nobody there to find it. It is the clearest illustration I have of what an unmaintained package actually costs: not that it stops working, but that a small wrong thing inside it gets nine years to sit there looking fine.

None of the work described above found it either. The test suite didn't, because the test covering this function asserted the *set of filenames* it produced and never which participant was inside each one — and it used IDs already in sorted order, where the bug cannot appear. The release gate didn't, because it never exercised that path at all.

**It was the AI that found it.** Not by scanning for bugs — I hadn't asked it to look. It was starting on an unrelated feature in that same function, read the line that picks each participant's trials and the line that names their output file, and pointed out that the two use different orderings. Then it demonstrated it: two participants with opposite responses, presented out of order, each written image compared against that participant's own data. `ci_a.png` held b's classification image and `ci_b.png` held a's.

I want to be careful about what that does and doesn't say. It is not that the machine is a better reviewer than a person. It is that nobody had read those two lines together in nine years, because there was nobody whose job it was, and reading old code carefully is exactly the labour that a maintainer stops having time for first. What the assistant supplied wasn't insight so much as attention — which turns out to be the scarce resource in unmaintained research software.

Both gaps are now closed: the regression test compares image contents against each participant's own data, and the release gate exercises the path with deliberately unsorted IDs and fails if the pixels ever move.

The honest lesson is the one in the section below about benchmarks. A green test suite is evidence about what you thought to check, and a test whose title says "writes one PNG per participant" while its assertions only count filenames is exactly the sort of thing that reads as covered and isn't.

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

**What worked well.** It was genuinely strong at the archaeology: reading twelve-year-old code, reproducing bugs from issue reports, tracing why a backward-compatibility branch was unreachable. It wrote the bulk of the test suite. It was good at the kind of care this project needed — when replacing a deprecated random-number call, it checked that the replacement drew from the random stream *identically* rather than just from the same distribution, across 150 seed and probability combinations. The obvious modern replacement would have diverged and silently changed every infoVal computed from a given seed. I would not have thought to check that.

It also found genuinely subtle bugs I'd never have gone looking for. One: R's `load()` assigns objects straight into the calling function's frame, so an object stored in a saved file can silently overwrite a *function argument* of the same name. The other is the mislabelling correction above — nine years old, invisible to a green test suite, surfaced because it happened to read two adjacent lines while working on something else.

That second one is the honest case for doing this at all. I am not going to claim the assistant has judgement. But an unmaintained package doesn't usually fail from a lack of judgement; it fails from a lack of anyone reading it. Attention is the thing that ran out in 2017, and attention is the thing I could get back.

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
