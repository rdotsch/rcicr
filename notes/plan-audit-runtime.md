# Runtime verification of the rcicr audit

## Scope and evidence

Execute the proposed reproductions for #294 and #299–#303 against unchanged production R code at a9c696e971f15bd6c8fe29d1965304d525c2334b. This PR adds a diagnostic script and a step in the existing R check job. It does not fix or close those issues; their runtime evidence determines the subsequent fixes.

Source inspection verified that reference generation calls generateStimuli2IFC without use_same_parameters, reads original image paths, and resets the response RNG through stimulus reconstruction. DECISIONS.md's first-base equality observation does not establish the reference for later independent bases. The stimulus return is inside the per-base loop. Trial selection and participant grouping lack the validation described by #300 and #294. No R execution has yet occurred in this environment.

## Execution

1. After this plan is reviewed, add tools/audit-runtime.R. Use temporary directories and synthetic PNGs only. Install the checkout before running it so parallel workers load the same package.
2. Run the script as a final step in the existing R-CMD-check job, preserving all job names, triggers, checks and release gates. Run on the existing Linux, macOS and Windows matrix, with explicit core limits. Print a machine-readable result line per measurement and R/session details; persist a Markdown report in the runner temporary directory and upload it even if a diagnostic case fails.
3. Each case records observed values, warnings and errors. Unexpected setup failures are distinguished from reproduced defects and make the diagnostic step fail. Assertions establish the necessary fixture conditions and independent comparisons; diagnostic observations of bugs are not permanent regression tests asserting buggy behavior.
4. Delete this plan in the same branch. Keep the PR a draft while collecting evidence. Report measured results on the tracking issues and in the PR; do not merge or close issues as if they were fixed.

## Measurements

- #299: two bases, independently drawn saved parameter matrices, fixed stimulus and response seeds. Construct each base's pixel-by-trial noise from saved p and stimuli_params, simulate the same response draws, and compare the current reference vector against both oracles. Report norm differences, reference median/MAD, and InfoVal differences for an explicit second-base synthetic CI. Check shared-parameter and first-base controls, default response stream, and ncores 1/2. Treat numerical differences as fixture-specific, not estimates of prevalence or general statistical impact.
- #301: retain a valid stimulus file but remove the source PNG; separately use a uniform grey source generated with contrast maximization disabled. Record reference-generation errors and show the saved basis/parameters still produce finite noise.
- #300: fractional and zero stimulus IDs with 12 columns. Compare against the requested valid integer control and capture warnings. Demonstrate row/response misalignment, not merely a successful return.
- #294: short participant IDs dividing the trial count; compare against explicit recycled IDs and the intended correctly sized grouping. Record differences and warnings with serial and parallel generation.
- #302: two bases and two trials, PNG saving plus returned data frame. Count eight expected versus actual files, enumerate missing base labels, and compare the returned first-base noise against saved parameters. Cover shared/independent parameters and ncores 1/2.
- #303: opposing repeated responses producing exact zero raw CI; record nonfinite scaled/combined pixels. Cover autoscale on all-zero and mixed CIs, plus a masked zero CI without confusing intended masked NA with unmasked NaN.

## Main risk

An oracle that mirrors the defective first-base selection could falsely clear #299. Select the second saved matrix explicitly, assert its noise differs from the first, and use identical response draws for both. Small simulations can have zero MAD; use a sufficiently varied deterministic fixture, assert positive MAD, and report the actual iteration count. No fixture, numeric baseline, API or production R file is changed.
