# Plan: give the release gate the entry points added since 1.5.0

## Why

The gate compares this tree's numbers with a released version's, which is the one thing the test suite cannot do: a test pins values this repository computed for itself. Its battery only calls what the reference version can run, so everything added since 1.5.0 is compared with nothing today:

- `reference_stimuli` and `reference_method = "images"` on `computeInfoVal2IFC()` (#349, #354);
- `batchComputeInfoVal2IFC()` (#85);
- replaying a reference under the file's recorded `rng_kind` from a session of another kind (#315);
- `participants` on the batch CI functions (#87);
- the `scaling` attribute's constants (#9). The gate compares list fields, never attributes, so a regression in a recorded constant passes it today.

One older path is missing as well: InfoVal with a `response_seed`, which every version since 1.2.0 accepts (`git show v1.2.0:R/computeInfoVal2IFC.R`).

## The change

New extras in `tools/compare-harness.R`, each with a floor in its `SINCE` table so a reference version that cannot run them skips them, as `mask` and `zmap_plain` already do:

| extra | what it records | floor |
|---|---|---|
| `infoval_subset` | `computeInfoVal2IFC(reference_stimuli = <every other stimulus>)` | 1.6.0 |
| `infoval_images` | `computeInfoVal2IFC(reference_method = "images", force_gen_ref_dist = TRUE)` | 1.6.0 |
| `reference_norms` | the full vector `generateReferenceDistribution2IFC(save_rdata = FALSE)` returns, for the default and for a subset. An InfoVal reduces the reference to its median and MAD, so a change that alters the norms without moving those two would pass every InfoVal extra. | 1.6.0 |
| `reference_images` | the full reference vector computed the way rcicr 1.5.0 and earlier computed it, compared **exactly**: its key goes into `EXACT_RE` in `tools/compare-release-output.R` (line 346), not the eight-ULP tolerance numeric keys get. Side-specific: a version whose `generateReferenceDistribution2IFC()` has a `reference_method` formal is called with `"images"`; an older one is called with its default, which *is* that calculation. Both sides use `save_rdata = FALSE` and `response_seed = NULL`. This is what makes the promise in `?generateReferenceDistribution2IFC` (bit-for-bit with 1.5.0) a gate check against 1.5.0 itself, before 1.6.0 becomes the baseline. | 1.4.0, the first version whose default reference is built from the saved noise (#301), to be confirmed by running it against v1.4.0, v1.4.1 and v1.5.0 |
| `infoval_batch` | `batchComputeInfoVal2IFC()` over the full-set CI and the subset CI, per-CI `reference_stimuli` | 1.6.0 |
| `infoval_cross_kind` | under `RNGkind("L'Ecuyer-CMRG")`, for a file written under Mersenne-Twister: the full vector `generateReferenceDistribution2IFC(save_rdata = FALSE)` returns, and the InfoVal of the default CI with `force_gen_ref_dist = TRUE`, since by then the stored reference would otherwise be reused and nothing simulated under the other kind. The session's kind is restored afterwards. | 1.6.0 |
| `infoval_seeded` | `computeInfoVal2IFC(response_seed = 7)`, and the full vector `generateReferenceDistribution2IFC(response_seed = 7, save_rdata = FALSE)` returns, since the other vectors are all drawn without a seed | 1.2.0, to be confirmed by running it |
| `batch_participants` | `batchGenerateCI(participants = )` and `batchGenerateCI2IFC(participants = )`, every field: two exported wrappers that forward the argument separately. Its data give participants unequal trial counts within a group, so averaging per participant differs from pooling; the harness's existing `participants` and `batch` data are balanced, and a wrapper that stopped forwarding the argument would change nothing on them. The implementation checks that the two differ. | 1.6.0 |
| `scaling_record` | the **whole** `scaling` attribute, for each scaling method, for a `participants` CI, and after `autoscale()`. It is flattened by name (`unlist()` over the nested list, so `method`, `constant`, `individual.method`, `individual.constant.<id>`, `combined.method`, `combined.constant` and any field added later all appear). The flattened names and every character field (the method labels) are exact keys; the numeric fields go through the numeric comparison. Recording the attribute whole, rather than field by field, means a field that is dropped, renamed, relabelled or added changes a key. The input to `autoscale()` is a list of `'independent'` CIs, so the preserved constant is a number, not `NA`. | 1.6.0 |

The InfoVal extras go on the four existing `infoval` configs, where a reference costs seconds. One of them, `sinusoid-64-twobase-indep-infoval`, has bases with different parameter matrices, where every reference and InfoVal call stops without a `baseimage`. So every new call, the direct reference vectors and the batch call included, forwards the base the existing InfoVal step selects (`key`), under the same condition it uses: only where the function has a `baseimage` formal. `batch_participants` and `scaling_record` go where `batch`, `participants` or `individual_cis` already run.

**How they are compared.** The InfoVal-valued keys (`infoval_subset`, `infoval_images`, `infoval_batch`, `infoval_seeded`, and the cross-kind InfoVal) join the existing InfoVal branch of `compare_one()` (`tools/compare-release-output.R:562-565`). It applies the absolute `INFOVAL_TOL` of `1e-9` to existing InfoVals, not the eight-ULP relative tolerance other numbers get, which a harmless rounding difference between two O(1) z-scores can exceed. That branch takes `abs(a - b)` of a scalar; `infoval_batch` is a vector, so it becomes `max(abs(a - b))`, and its names are recorded as a separate exact key (`infoval_batch_names`, in `EXACT_RE`), since the result promises names and the numeric comparison ignores them. The reference vectors (`reference_norms`, the cross-kind vector) stay on the generic numeric comparison, and `reference_images` is exact, as above.

**When they start comparing.** CI's `reproducibility.yaml` compares against v1.0.1 and against the newest release tag (`.github/workflows/reproducibility.yaml:129-150`); a local `--ref` can name any commit. A floor of 1.6.0 means the extras run once that reference is 1.6.0 or later: CI's previous-release run from the first PR after 1.6.0 is tagged, and a local `--ref` at any `main` commit after it. Until then they are skipped, as the floors intend. The battery is chosen by the reference version, so nothing here can crash a run against v1.0.1 or v1.5.0.

**Not floored at a development version.** `1.5.0.9000` would make the extras run against `main` now, but `main` has carried that version both before and after these functions existed, so the floor could not tell those commits apart. A reference that cannot run an extra aborts the whole gate. Floors are release versions only.

## Verification

- `--quick` against v1.0.1: unchanged check count, since every new extra has a floor above 1.0.1.
- `--quick` against v1.5.0: the 1.6.0-floored extras are skipped, but `infoval_seeded` (floor 1.2.0) and `reference_images` (floor 1.4.0) run on each of the four InfoVal configs, so the count grows by those checks and no others. `reference_images` must come back identical, bit for bit; if it does not, that is a defect in `"images"`, and the plan stops there rather than adding an `EXPECTED` entry. The measured count before and after, and the new keys, are quoted in the PR.
- The driver cannot force the floor: `compare-release-output.R` reads the reference version from the reference tree's `DESCRIPTION` and passes it to both harness runs, overriding the environment (`run_side()`). So `tools/compare-harness.R` is run directly with `RCICR_COMPARE_REF_VERSION=1.6.0`, twice into separate output directories. Then every new extra's keys must be present in the saved output, with their counts listed in the PR, and the two runs must be identical. That proves each extra executes, records something, and is deterministic, which is what a comparison needs. It is a check of the harness only, never how the gate is run.
- `infoval_seeded` against v1.2.0, v1.4.0 and v1.5.0. Where a version's documented reference change (`NEWS.md`, the saved-noise rebuild of #301, the Gram route of #354) moves the value, it gets an `EXPECTED` entry with a `check` predicate bounding it, not a bare entry that would excuse any value. If a version cannot run it, the floor rises to the first that can.

## The step most likely to fail

`infoval_cross_kind`. It changes the session's RNG kind inside the harness, and the harness relies on one random stream for everything after it, which is why `infoval_oracle()` restores `.Random.seed`. The extra saves and restores the full `.Random.seed` around its call (which restores the kind too, as #315 established) and runs last in its config, and it records whether the caller's complete `.Random.seed` (kind and position) is identical immediately after the package calls, before the harness restores anything. That is an exact key, so a release that kept the kind but moved the caller's stream fails the gate.
