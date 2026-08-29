# Decisions: the latent-space module

Why the experimental latent-space reverse correlation behaves as it does.

Separate from [`DECISIONS.md`](DECISIONS.md), which was within 80 words of its budget when this module was written. Fold this in and drop the file when the module leaves the experimental gate behind.

## The module is exported and gated rather than unexported

Every function here refuses to run unless `options(rcicr.experimental = TRUE)` is set.

The gate is not to hide the code but to make a user state that they accept a weaker guarantee: the module's numeric output is not covered by the promise that a stored script reproduces its results, so a classification image computed today may not reproduce under the next development build.

Exported rather than internal, because an unexported function gets no manual page, no pkgdown entry, and no `R CMD check` of its examples, and the checking is the point of having continuous integration on it at all.

The gate comes off once a golden master pins the numbers and a release ships them.

## The generator is a contract, not a dependency on torch

R cannot run a StyleGAN. A hard dependency on torch was rejected: it downloads libtorch on first use, would be an Import for a feature most users never touch, and still would not run a StyleGAN without a TorchScript export the user makes themselves.

A generator is instead a list satisfying a documented contract, which `latentGeneratorPCA()` does in base R. That is what lets the module be run, tested and checked with no GPU, no Python, no network and no additional package. An eigenface model is far weaker than a generative adversarial network; it is not there to be a good one, but so the arithmetic is verifiable.

This is the first `class()` in `R/`, where every other object is an unclassed list. A generator is the one thing a user supplies and the package checks, so it needs a type. The classification image is still a plain list.

## `latentGeneratorTorch()` was planned and not written

A TorchScript backend through the torch package was planned and is not here (#291). Nothing available can install torch, so it would ship with no test at any level, and `latentGeneratorCommand()` already reaches such a generator through a short script while being tested end to end against a real subprocess.

## The stimulus file stores a fingerprint, not the model

A multi-gigabyte generator cannot go in an `.Rdata`, so the file records a string identifying the renderer and `generateLatentCI()` refuses one that does not match. Rendering through the wrong model gives a face that looks plausible and belongs to nobody, which is why it errors rather than warns.

`latentGeneratorPCA()` is small enough that its components travel inside the file, so it rebuilds from the file alone months later, as `base_faces` makes the pixel-noise pipeline self-contained.

It is summary statistics at full precision, not a digest: base R has no hashing function and none of the Imports provides one, and a determined caller can construct a collision. For an external generator, passing the weights file gives a fingerprint that follows the model, via `tools::md5sum()`.

## The recovered direction is left weighted by the sampling covariance

Under Gaussian perturbations of covariance S, the response-weighted mean estimates the participant's preference direction `u` multiplied by S, not `u` itself. Dividing S back out was considered and rejected: the whitened direction fits the responses equally well and renders as a face nobody has, because it travels furthest along the dimensions on which real faces vary least. Leaving the weighting in is what keeps the answer on the face manifold. Two participants measured through one generator are weighted identically, so their results stay comparable to each other.

## A 2IFC response carries no distance, so agreement between generations supplies it

A sign says which of two faces is preferred and not by how much. The consequence is exact rather than approximate: the response-weighted mean converges to `sqrt(2/pi) * S %*% u / sqrt(t(u) %*% S %*% u)`, whose length in standard deviations is the same whatever the length of `u`. A single round of trials therefore cannot say how far along the estimated direction to go, which is why the scaling constant in `generateLatentCI()` is the researcher's choice and not something the data can settle.

Two rules for `searchLatent2IFC()` were tried and rejected before the one that is there.

A step proportional to the estimate does not converge: it moved about one standard deviation every generation however close the centre already was, circling the answer while the distance bounced between 0.3 and 1.0.

A fixed schedule of `alpha / generation^step_decay` converges only when `alpha` happens to suit the distance, which is exactly what is unknown. Measured against simulated observers, `alpha = 1` ended 0.12 standard deviations from a target 1.6 away and 12.23 from one 15.6 away, while `alpha = 3` reversed the ordering. No default can be right for both.

The distance is in the sequence of estimates rather than any one of them: generations that agree mean the centre is still short, generations that disagree mean it went past. The step grows while they agree and shrinks when they do not, which settles because `step_grow * step_shrink` is below 1, so an oscillation cannot sustain itself. Against the same observers this reached 0.75 and 3.48 where the fixed schedule reached 2.92 and 12.23, at no cost near the target.

The same sequence is what the `cosine_with_previous` diagnostic reports, so the number that drives the step is also the number a researcher reads to see what the search did.

## Whitening the search direction was tried and is worse

For a direction of travel the covariance weighting looked like a liability: on a generator whose components differ in scale by 65 to 1, the search barely moves along the narrow ones. Dividing it back out measures worse at every distance tried, ending 1.91, 6.32 and 15.61 standard deviations out where leaving it in ended 0.77, 3.89 and 12.93.

An anisotropic face space is genuinely harder to search, and the answer is a better set of images rather than a correction afterwards. `latent_sd` is on the generator so this is visible before a study runs: a first component dwarfing the rest means most variation lies along one direction.

## Two meanings of "so many standard deviations", and where each belongs

`applyLatentScaling()` measures a displacement as a root mean square across dimensions, so `scaling_constant` means the same on generators of different dimensionality. `searchLatent2IFC()` uses a Euclidean length, because a search direction is usually concentrated in a few dimensions and dividing by their number would inflate the step by up to `sqrt(latent_dim)`.

Writing both as a root mean square made a search step overshoot by exactly that factor: `cosine_with_previous` was -1 in every generation and the search ended further from its target than it began. `latentNorm()` follows `applyLatentScaling()`; informational value is unaffected either way, being a ratio in which a constant factor cancels.

## `sigma_decay` defaults to 1, against expectation

Shrinking the perturbations as the search proceeds sounds like it should sharpen the answer. Measured against observers with internal noise it is worse at every level, and worst for the noisiest: it shrinks the difference between the two faces while the participant's uncertainty stays put, so responses get noisier exactly as the search needs them finer.

The argument stays, because a real participant's discriminability is not the simulation's. The default does not.

A noiseless simulated observer cannot inform this at all: its responses depend only on the sign of a projection, which is unchanged by scaling the perturbations, so `sigma_decay` has mathematically no effect on such a run. Any tuning of it needs an observer with internal noise, and anyone re-measuring it should check that first, because a sweep against a noiseless observer returns identical numbers for every value and reads as a bug in the sweep.

## Responses are checked, and the module is stricter than the pixel pipeline

`generateCINoise()` documents that a response "can be changed into a scale", and the pixel pipeline accepts whatever it is given. The latent module does not: 1 and -1 or an error.

The case that decides it is a 0/1 coding, which several experiment programs write. It raises no error anywhere downstream. It turns the response-weighted mean into a quantity with no meaning and still renders a perfectly plausible face, so nothing about the result looks wrong. Nothing in this module depends on the older latitude yet, so a scale can be added later as a decision rather than inherited as an accident.

One check, `checkResponseCoding()`, at every entry point that takes responses. Callback mode had it and the resumable path did not, which is the wrong way round: the resumable path is the one fed by another program.

## The null shuffles the observed responses rather than drawing fresh ones

Written first as `((runif(n) > 0.5) * 2) - 1`, copied from `generateReferenceDistribution2IFC()` where it is correct: that function simulates an observer from scratch, so balanced coin flips are exactly what it wants. Here the documentation promised a permutation null and the code did not deliver one.

The difference is not cosmetic. A participant who pressed one key more often than the other would be compared against balanced null vectors they could never have produced, and the imbalance would be reported as latent signal. In the limit, an all-one response vector has a single possible permutation and must sit at the centre of its own null; against fresh coin flips it scores as though it carried information.

Holding the response multiset fixed and breaking only its pairing with the stimuli is the question being asked, so `computeLatentInfoVal2IFC()` takes the responses as an argument.

Within each participant, when there are participants. A global shuffle moves answers between people, so one who always pressed one key and one who always pressed the other would be compared against a null full of mixed answers neither could have given, and the difference between their key biases would read as signal.

That case also makes the null degenerate: one arrangement means no spread, and the standardised score is 0/0. The limit is reported rather than `NaN`, so the headline number stays readable for exactly the case the permutation null exists to handle.

## One direction computation, shared by the estimate and its null

The null was built by its own copy of the arithmetic, and the copy drifted. It pooled where the classification image had averaged per participant, letting the participant with more trials decide the score. And it permuted responses that pooling had already collapsed over repeated stimuli, shuffling means rather than answers: `c(1, 1, 2)` answered `c(1, 1, -1)` collapses to `c(1, -1)`, whose two arrangements are not the three arrangements of the answers, collapsed.

Both are one mistake with one fix: `latentDirection()` holds the whole path from trials to a direction, and the null calls it per iteration on a trial-level permutation, so the two cannot diverge again without the estimate changing too.

Each iteration re-aggregates, so the null scales with trial count as well as `iter`: a few seconds for 300 trials at the default, far below the pixel-noise equivalent.

## Informational value is bound to the analysis, not the stimulus set

The first version compared generator fingerprints. Two experiments run through one generator have identical fingerprints by design, which is the whole point of a fingerprint that identifies a renderer, so that check accepted a classification image from one experiment scored against another's perturbations.

Widening it to the stimulus file was not enough either: one file's trials can be analysed many ways. The fingerprint now covers the file, the trials analysed, the responses and the participant structure, so it identifies the analysis rather than the material. It is derived on demand, so the `.Rdata` contract is unchanged.

## A resumed search restores its design and its random stream

Resumable mode took its settings from whichever call supplied them, so a search begun at `latent_sigma = 0.5, sigma_decay = 0.8` and resumed by the documented call, which names only the state and the responses, continued at sigma 1 rather than 0.4. `set.seed(seed)` likewise ran only when the first state was written, so a search resumed in another session drew from that session's unrelated stream.

Both now travel in the state, along with the seed, which names the files. Every setting is restored, not only those the current generation reads, because they are written back into the next state and one left behind reverts a resume later. The tests that missed this passed every argument on every call, which no documented usage does.

## Informational value uses a permutation null, not the erratum formula

`computeInfoVal2IFC()` compares a classification image against a simulated reference distribution of noise images, and its definition is fixed by the published erratum. `computeLatentInfoVal2IFC()` permutes the observed responses over the same trials instead.

That null needs no rendering, which is why nothing is cached in the stimulus file and there is no interactive prompt, where the pixel-noise equivalent is expensive enough to need both. Length is measured in standard deviations of the generator's training faces so a generator whose components have very different scales does not have its widest dimension decide the answer alone.

The two numbers are not comparable with each other, and the documentation says so rather than leaving a reader to assume a shared scale.

## Rendering does not go through the foreach backend

The rest of the package parallelises with `parallel` and `foreach`. This module renders in the calling process.

A generator is a single external resource: a GPU, a Python interpreter, a loaded model. Fanning it across workers either duplicates the model in memory or contends for the device, and neither is an improvement. The speed is in batching instead, which is why `render()` takes a matrix of latents rather than one at a time, and why `generateStimuliLatent2IFC()` has a `batch_size` rather than an `ncores`.

## The antithetic pair was kept

Each trial renders the base latent plus a perturbation and the base latent minus it, rather than two independent samples.

It keeps the response coding of the existing pipeline exactly, so a task script producing 1 for the original image and -1 for the inverted one needs no reshaping. The pair is also antithetic, which halves the variance of the estimator built from it.
