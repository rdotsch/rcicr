# Plan: latent-space reverse correlation, on an experimental branch

Tracking issue: #289.

## What this adds, and why it is a second pipeline rather than a change to the first

Reverse correlation in this package works in one basis: a fixed stack of sinusoid or Gabor patches from `generateNoisePattern()`, with the classification image computed by `generateCINoise()` as a response-weighted mean of per-trial contrast weights. That basis is unconstrained, so most of its parameter space is not a face, and the classification image comes out as a greyscale field rather than as a person.

Albohn, Uddenberg and Todorov (2022, Behavior Research Methods) replace the basis with the latent space of a generative face model. Perturb a base latent, render two faces per trial, average the perturbations the participant chose. The arithmetic is the same as `generateCINoise()`. Only the basis and the renderer change. Every point in the space is a face, so the estimate stays on the face manifold, and because the space is low dimensional it can be searched across generations instead of sampled once.

This lands as a parallel pipeline. No existing function, argument, saved field or number moves, with one exception described under "The one risky commit" below.

## Where it runs

Everything is on `claude/reverse-correlation-gan-faces-he3i4f`, with a draft pull request into `main` that is held as a draft rather than merged. The five slices below are commits on that one branch.

No workflow change is needed and none is made. `R-CMD-check.yaml`, `lint.yaml`, `reproducibility.yaml`, `test-coverage.yaml` and `pkgdown.yaml` all trigger on `pull_request: branches: [main]`, and `grep -rn "draft\|types:" .github/workflows/` returns only `pkgdown.yaml:9`, which is the `release: types: [published]` trigger. So a draft pull request from this branch gets the full required suite on every push. The `AGENTS.md` warning about stacked pull requests getting pre-commit and nothing else applies to a pull request whose base is not `main`, which is not the shape used here.

## The experimental gate

`R/experimental.R` provides `requireExperimental()`, called as the first statement of every exported function in this module. It errors unless `options(rcicr.experimental = TRUE)` is set.

Exported and gated rather than unexported, because unexported code gets no manual page, no pkgdown entry, and no `R CMD check` of its examples, and the checking is the reason for having continuous integration on this branch at all.

The gate carries a carve-out that has to be stated rather than assumed. While this module is experimental its numeric output is not covered by the package's guiding constraint, so a latent classification image computed today may not reproduce under the next development build. Setting the option is how a user consents to that. The gate is removed only once a golden master pins the numbers and a release has shipped them.

## The generator contract

R cannot run StyleGAN. The method is separated from the renderer so the method stays testable.

`R/latent-generator.R` defines an `rcicr_generator`: a named list with `kind`, `latent_dim`, `img_size`, `space`, `latent_mean`, `latent_sd`, `render`, `fingerprint` and `state`, carrying a class attribute. `render()` takes a matrix of latents, one row per stimulus, and returns greyscale pixel values in `[0, 1]`, batched because an external renderer pays a large per-call startup cost.

This is the first `class()` and the first `inherits()` in `R/`, which is a `DECISIONS.md` entry rather than a detail.

`fingerprint` is the reproducibility hinge. A multi-gigabyte generator cannot go in the `.Rdata`, so the file stores this string instead and `generateLatentCI()` errors when handed a generator that does not match the one that made the stimuli. `latentGeneratorPCA()` is small enough to embed its components, and does, so that pipeline is self-contained in the way `base_faces` makes the existing one self-contained.

`latentGeneratorPCA()` builds an eigenface generator in base R from a set of aligned face images. It exists so the whole feature runs, tests and checks with no GPU, no Python and no network. `torch` and `reticulate` are `Suggests` reached through `requireNamespace()`, which is the package's first optional dependency and so is also a `DECISIONS.md` entry: there is currently not one `requireNamespace()` call in `R/`.

## The `.Rdata` contract

A new file kind, `<label>_latent_seed_<seed>_time_<ts>.Rdata`, not an extension of the existing one. Fields: `latent_params`, `base_latent`, `latent_sigma`, `generator_spec`, `n_trials`, `img_size`, `seed`, `label`, `stimulus_path`, `generator_version`.

Three existing rules govern it. It is append-only. `load()` assigns into the calling frame, so every new field name is checked against every argument of every function that loads the file and every new argument against every field, with argument copies kept across the `load()` through `captureArgs()` as `computeInfoVal2IFC.R:88-90` does. Every new loaded name goes into `R/zzz.R`.

Reading goes through a new `loadLatentStimulusParams()` in `R/rdata.R`, beside the existing `loadStimulusParams()`.

## Slices

1. `R/experimental.R`, `R/latent-generator.R`, `R/latentGeneratorPCA.R`, and the extraction described below.
2. `generateStimuliLatent2IFC()`, the latent `.Rdata`, `loadLatentStimulusParams()`, the `zzz.R` entries.
3. `generateLatentCI()`, `applyLatentScaling()`, `computeLatentInfoVal2IFC()`, and the recovery test.
4. `searchLatent2IFC()`, in callback and resumable modes, with convergence diagnostics.
5. `latentGeneratorCommand()`, `latentGeneratorTorch()`, a helper script under `inst/python/`, and an article under `vignettes/articles/`.

## The one risky commit

Slice 1 extracts the base-image reading, greyscale conversion, size check and contrast maximization out of `generateStimuli2IFC.R:62-128` into `R/base-images.R`, and calls the extracted helper from both pipelines.

That is the only edit to existing code in the whole plan, and it is the step most likely to fail, because it is the only one that can move a number. Two instruments decide whether it did: `test-regression-baseline.R`, which pins the default pipeline's output, and `tools/compare-release-output.R --quick`, which runs the released v1.0.1 code against this tree. If either goes red the extraction is wrong and gets reverted. Neither is adjusted, and no `EXPECTED` entry is added, because a correct extraction produces no difference to record.

The extracted helper takes the calling function's name so the two advice strings that name `generateStimuli2IFC()` stay accurate for the existing path. `test-error-paths.R:464-491` matches on substrings rather than whole messages, so it would not catch a wording change on its own.

## What decides whether the adaptive mode ships

`searchLatent2IFC()` claims that re-centring across generations reaches a given representation in fewer trials than sampling once around a fixed centre. That is a claim about behaviour, so it is measured rather than asserted.

`tests/testthat/test-latent-recovery.R` builds a PCA generator from synthetic images, fixes a hidden target latent `t`, and answers trials with a simulated observer choosing whichever of the pair is closer to `t`. It asserts that the recovered direction correlates with `t - w0`, that the correlation falls to chance under permuted responses, and that the adaptive search reaches a stated cosine to `t` in fewer total trials than the one-shot method. If that last comparison does not hold, slice 4 does not ship and the module is one-shot only.

## Verification

- `devtools::install()` before anything that renders, because the render loop makes workers call `library(rcicr)`.
- `devtools::test()` green, including the new contract, gate, load-frame and recovery tests.
- `test-regression-baseline.R` green and `tools/compare-release-output.R --quick` green with no `EXPECTED` entry. This is the check that decides whether the branch has touched production behaviour.
- A latent classification image rendered by hand from 100 simulated trials, confirmed to be a face rather than a smear.

## Records

`NEWS.md` leads with the module being experimental and its numbers not yet stable, and carries no "Reproducibility impact" section, because nothing existing changes. `DECISIONS.md` gains entries for the generator contract, the fingerprint, the permutation null, the antithetic pair, the first S3 class, and the gate. `DECISIONS.md` has a 5200 word budget and `AGENTS.md` 2800, checked with `wc -w`, so an entry that pushes either over means something comes out rather than the budget going up.
