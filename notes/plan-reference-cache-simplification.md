# Reference cache simplification

The current shared and independent paths repeat stale-cache detection, iteration inheritance, warning handling, read-only fallback, and RNG preservation. R/reference-base.R fingerprints only three selected norms; R/computeInfoVal2IFC.R separately implements the same migration policy.

## Changes

- Introduce one internal reference-cache resolver accepting the cache values and provenance plus the existing public generation arguments. Both scoring paths delegate generation and migration to it; storage remains in the existing generators.
- Preserve cached iteration counts and ambient RNG state for automatic refreshes, explicit force/response_seed semantics, seeded-cache reuse, read-only refreshes, and numeric scoring. Leave the dormant lookup machinery intact.
- Hash the complete vector with digest::digest using an explicit serialization version and algorithm. Declare digest in Imports. Existing development fingerprints fail validation and refresh once under the same migration policy.
- Use a message consistently when a changed reference supersedes a cache; remove warning-to-error interception. Update notification tests and NEWS to state this policy.
- Trim comments to the compatibility constraints and maintain the decision rationale in DECISIONS.md.

## Verification

The riskiest step is extracting the shared path around load(), which must not overwrite caller arguments or reuse saved scratch fields. Keep the argument guards and test both storage formats. Retain the existing cache, RNG, seeded-null, read-only, and release-comparison regressions. Add a fingerprint test changing a norm outside the formerly sampled positions, retaining the format-option checks. Run tests and lint, then the existing CI/reproducibility gate on the committed branch. Local R is not currently installed; no local runtime results have been claimed.

Remove this file with the implementation so no plan remains in the squash.
