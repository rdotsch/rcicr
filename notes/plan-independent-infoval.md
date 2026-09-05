# Correct independent-base InfoVal references (#299)

## Verified evidence and scope

PR #304 ran the audit against production commit a9c696e971f15bd6c8fe29d1965304d525c2334b on Linux release/devel, macOS and Windows. With two independent saved parameter matrices, the current reference followed base1; a base2 synthetic CI scored 2.853360613 rather than 2.920465250 under the base2 reference. Shared-parameter and serial/parallel controls passed. The formula itself is correct and is not changed.

This PR fixes #299. It uses saved noise for independent bases and therefore removes the original-image dependency for that path. The general shared/legacy reconstruction change in #301 is separate: retaining the shared route avoids mixing its older missing-metadata fallbacks and RNG behavior into this correction.

## API and compatibility decision

Append baseimage = NULL to generateReferenceDistribution2IFC and computeInfoVal2IFC, preserving every existing positional argument. A single base or identical saved parameter matrices needs no selector and uses the existing shared reference route. If matrices differ, require one valid base label; omission produces an actionable error listing labels instead of silently selecting a reference. This is an intentional tightening of ambiguous calls, documented under Reproducibility impact. Do not infer the base from pixel equality: distinct bases can use the same image.

Use the saved parameter matrices to detect differing noise, not just use_same_parameters or unreliable generator_version metadata. Explicit invalid selectors fail early. Keep generateCI's return object and all existing functions' numeric formulae unchanged.

## Implementation

1. Read reference selection into an isolated environment before either function's legacy load. Return only the selected independent context when needed, so shared calls keep the existing cache/reconstruction path. Protect the appended argument across legacy loads and exclude new scratch variables from saves.
2. Build the selected independent base's pixel-by-trial noise from its saved p and stimuli_params. Honor ncores using the existing foreach/backend pattern. Normalize the legacy s/sinusoids layout and 4096-column parameters through the existing compatibility semantics.
3. Simulate the same runif response coding and compute norms with the existing matrix multiplication/norm expression. With response_seed supplied, seed responses only. With NULL, reproduce the historical shared-reconstruction response stream: set.seed(saved seed), consume one matrix's trial-by-trial uniform draws at saved nscales (legacy default 5), then simulate responses. Do not consume extra draws for additional bases. Verify this against the old public generator, not merely a duplicate draw-count expression.
4. Append reference_norms_by_base, a named list of entries containing norms and response_seed. Independent calls never reuse the old unscoped reference_norms and never overwrite that legacy vector. Updating one base preserves other bases and all original stimulus fields. Explicit response_seed in computeInfoVal2IFC still forces a fresh, uncached draw; force_gen_ref_dist regenerates only the requested base. Leave the shared path's legacy reference_norms behavior intact.
5. Document selector usage, cache ownership and affected calls in roxygen/man and NEWS. Correct DECISIONS.md's first-base-only claim without broadening the word budget. Delete this plan from the branch.

## Validation

Add focused regression tests before production changes and execute them against unchanged R code in the PR to demonstrate failure. Then implement on the same branch and verify they pass. Tests cover both independent bases against saved-noise oracles under identical draws; positive MAD and a deliberately different wrong-base reference; explicit/default response seeds; serial/parallel equivalence; shared default equality to the old generator; selector omission/invalid names; per-base cache hits and force regeneration; stale unscoped cache rejection; no-write seeded calls; argument collisions; moved sources for the independent path; and legacy basis/parameter layouts. Do not alter golden values or released-signature fixtures.

Run the existing four-platform R checks, lint, documentation regeneration and both release comparisons in GitHub Actions. Local R is unavailable. Generated man changes must match pinned roxygen2 7.3.1; use the remote regeneration output if needed. Record measured output changes without adding EXPECTED entries unless the release comparison actually reaches an affected configuration.

## Most likely failure

Reproducing the historical response RNG offset while replacing stimulus reconstruction. Shared results and default first-base independent results must agree with the old pipeline within the existing floating-point tolerance. A test compares against generateStimuli2IFC's actual resulting RNG state. Cache selection is the other silent-failure boundary: a base2 call must not accept base1 or unscoped norms, and selected-base provenance must survive reloading and saving.

No merge or release is part of this implementation step.
