# Plan: let researchers match the InfoVal reference to the stimuli they scored (#349)

## Position

Matching the reference to the target CI's design is the **researcher's responsibility**. The package
makes the paper-covered case easy to ask for and changes no default. Brinkman et al. (2019, <https://doi.org/10.3758/s13428-019-01232-2>, Part I,
after Eq. 1) require the reference to use the identical stimulus set, including its number of
stimuli. `analyses/infoval-design-mismatch.md` measures what ignoring that costs: at 512px, the median pure-noise CI's InfoVal was 0.09 to 0.11 with 1% of trials missing and 0.51 to
0.58 with 5%. Each range spans the four stimulus sets measured, two with 770 trials and two with 300.

## Scope

**In:** an optional `reference_stimuli` argument naming the distinct saved stimuli a CI was built
from, with one response each. This covers missing trials and partial trial sets. It deliberately
does not reuse the name `stimuli`: in `generateCI()` that argument holds one stimulus number per
response and "Numbers may repeat" (`R/generateCI.R:52`), whereas this one is a set of distinct
stimuli.

**In:** a message guard. It never warns or errors, and changes no number (see "Message guard").

**Out:**
- **Repeats and participant averages.** The paper gives no reference for them, and defining one is
  a methodological decision this change does not make. The help page points to per-participant
  InfoVal for full-set participants, which the paper does cover.
- **Any change to a default.** Every existing call returns the same number, bit for bit.
- **#354 (Gram-matrix reference).** It is not needed: a subset reference renders only the subset's
  images, so it costs at most what a full reference costs today.

## Changes

1. **`generateReferenceDistribution2IFC(..., reference_stimuli = NULL)` and
   `computeInfoVal2IFC(..., reference_stimuli = NULL)`**, added last, so positional calls are unaffected.
   - The value must be distinct whole numbers in `1..n_trials`, validated as `generateCI()`
     validates its `stimuli` (`validateStimulusIds()`). Duplicates are an error that says repeats
     are outside what the reference covers.
   - The validated value is then made canonical, `sort(as.integer(x))`, before any comparison. IDs
     often arrive as doubles (`unique()` of a numeric column), and `c(1, 2, 3)` is not `identical()`
     to `seq_len(3)`. The canonical form serves the full-set check, cache matching and the guard's
     attribute alike. The order does not matter, since the null does not depend on stimulus order.
   - `NULL`, or a canonical value equal to `seq_len(n_trials)`, takes today's path unchanged.
   - The reference is built exactly as today (`referenceNoise()` then `seedResponseStream()`), but
     from the selected columns only. Each iteration draws `length(reference_stimuli)` responses. Without a
     `response_seed`, the stream replays the generator's parameter draws for all `n_trials` first,
     as it does now, so a subset reference is reproducible from the file alone.
   - Independent bases (`computeBaseInfoVal()` / `generateBaseReference()`) take the same argument.
     `savedReferenceParams()` gains the selection.
2. **Cache.** A new, append-only `.Rdata` field `reference_norms_by_stimuli`: a list of entries
   `list(reference_stimuli, baseimage, norms, response_seed, source, fingerprint)`, matched with
   `identical()` on `reference_stimuli` and `baseimage`. It uses the same trust policy as `resolveReferenceNorms()`. The
   full-set `reference_norms` and `reference_norms_by_base` are never read or written on the subset
   path. The empty `ref_lookup` table, keyed on `n_trials`, is skipped for subsets.
3. **Docs.** `?computeInfoVal2IFC` and `?generateReferenceDistribution2IFC` state the matching
   requirement with the citation, give the measured sizes from the analysis as it states them (medians, as ranges over the stimulus sets measured), and show the one-line
   use:
   `computeInfoVal2IFC(ci, rdata, reference_stimuli = unique(my_stimuli))`.
   The README's InfoVal section gets one paragraph.
4. **`NEWS.md`**: a new-feature entry. There is no "Reproducibility impact" entry, because no
   existing call changes. **`DECISIONS.md`**: why matching is left to the researcher, and why
   repeats and participant averages get no argument.

## Message guard

`generateCI()` records its trial design as **an attribute** on its return value,
`attr(ci, "trial_design")`, not as a list field. The list's names are pinned
(`test-generateCI.R:29`, `test-batchGenerateCI.R:20`), and scripts may iterate its fields expecting
matrices. `generateCI2IFC()` and `batchGenerateCI()` return `generateCI()`'s value, so they carry
the attribute too. The attribute holds three things:

- `stimuli`: the sorted distinct saved stimuli used, in the same canonical integer form;
- `repeated`: whether any stimulus was presented more than once, pooled or within a participant;
- `n_participants`: 1 when `participants` is not given.

Distinct IDs alone would not be enough. A repeats or participant-average design would look like a
subset, and the guard would recommend `reference_stimuli` for a design no reference covers.

`computeInfoVal2IFC()` prints one `message()`, in one of two forms:

- **One responder, no repeats, and a stimulus set that differs from the reference's.** The
  reference's set is every saved stimulus when `reference_stimuli` is omitted, otherwise the ones
  passed. The message names the number of stimuli in each and shows the one-line use.
- **Repeats, or more than one participant.** Whatever `reference_stimuli` says, the message states
  that every reference here assumes one response per stimulus from one responder, and that the
  package defines none for this design. It points to per-participant InfoVal and never recommends
  `reference_stimuli`.

It never warns or errors and changes no number, so a researcher who knows their design can ignore
it or silence it with `suppressMessages()`. A CI without the attribute (from an older version, or
built by hand) gets no message, because there is nothing to compare.

## Risk: the step most likely to fail

`generateReferenceDistribution2IFC()` re-saves **its whole frame** into the `.Rdata`, excluding a
hand-kept `internals` list, and `load()` assigns into the same frame. The name
`reference_stimuli` removes the collision with the existing local `stimuli` (the noise matrix), but
not the other two ways the argument can go wrong:

- a same-named object in a file could override it;
- it could be written back into the user's file as if it were stimulus metadata.

Excluding the name from the save would fix the second, but it would silently delete any
`reference_stimuli` object a file already holds, which breaks the append-only `.Rdata` contract.
Mitigation:
- read the argument only from `.args`, which is captured before `load()`;
- decide from `load()`'s return value, the names it loaded (`loadRdata()` passes it through):
  - if the file supplied `reference_stimuli`, the frame keeps the file's value and the save writes
    it back unchanged;
  - otherwise the formal is removed from the frame before the save;
- test both cases:
  - a file carrying a `reference_stimuli` object: the argument is unaffected and the object is
    saved back `identical()`;
  - a file without one: the save adds only the cache field.

The independent-bases path (`generateBaseReference()`) already loads into and saves from a separate
environment, `selection$source` (`R/reference-base.R:202`), so it needs only the cache field.

## Verification

- **Default unchanged:** with `reference_stimuli` omitted, and with `reference_stimuli = seq_len(n_trials)`, norms are
  `identical()` to `main`'s. The golden master stays green. The release gate against `origin/main`
  reports `0 expected deviations`, using the `--ref` method in `CLAUDE.md`.
- **Subset correct:** a subset reference equals a direct random-response simulation over those
  stimuli, with the same stream, to 1e-12. This is the Gram check from the analysis, run at 64px.
- **Paper case end to end:** a CI from `generateCI()` on a subset, scored with `reference_stimuli =`, gives a
  median pure-noise InfoVal near 0, where the default gives the inflated value.
- **Cache isolation:** a subset call never alters `reference_norms`; a second identical call hits
  the cache; a different subset misses it.
- **Validation errors:** duplicates, out-of-range values, non-integers, and an empty vector.
- **Canonical form:** `c(3, 1, 2)` as doubles takes the full-set path with `n_trials = 3`, and a
  subset given once as doubles and once as integers hits the same cache entry.
- **Guard:**
  - the subset form fires for a one-responder subset CI scored without `reference_stimuli` or with a
    different set;
  - the unsupported-design form fires for pooled repeats, for a participant average over the full
    set, and for either scored with `reference_stimuli`, and never mentions the argument;
  - it stays silent for a one-responder full-set CI, when the sets match, and on a CI without the
    attribute;
  - the attribute survives `generateCI2IFC()` and `batchGenerateCI()`, and does not change
    `names()` of the returned list.
