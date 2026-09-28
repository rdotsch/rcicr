# Plan: let researchers match the InfoVal reference to the stimuli they scored (#349)

## Position

Matching the reference to the target CI's design is the **researcher's responsibility**. The package
makes the paper-covered case easy to ask for and changes no default. Schmitz et al. (2020, <https://doi.org/10.3758/s13428-019-01232-2>, Part I,
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
   - `NULL`, or exactly `seq_len(n_trials)` after sorting, takes today's path unchanged.
   - Otherwise the value must be distinct integers in `1..n_trials`. Duplicates are an error that
     says repeats are outside what the reference covers. The value is sorted to a canonical form,
     since the null does not depend on stimulus order.
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

`generateCI()` records the sorted distinct saved stimuli it used as **an attribute** on its return
value, `attr(ci, "reference_stimuli")`, not as a list field. The list's names are pinned
(`test-generateCI.R:29`, `test-batchGenerateCI.R:20`), and scripts may iterate its fields expecting
matrices. `generateCI2IFC()` and `batchGenerateCI()` return `generateCI()`'s value, so they carry
the attribute too.

`computeInfoVal2IFC()` prints a `message()` when the attribute exists and differs from the stimuli
the reference is built over: every saved stimulus when `reference_stimuli` is omitted, otherwise
the ones passed. The message names the number of stimuli in each and the one-line use. It never
warns or errors and changes no number, so a researcher who knows their design can ignore it or
silence it with `suppressMessages()`. A CI without the attribute (from an older version, or built
by hand) gets no message, because there is nothing to compare.

## Risk: the step most likely to fail

`generateReferenceDistribution2IFC()` re-saves **its whole frame** into the `.Rdata`, excluding a
hand-kept `internals` list, and `load()` assigns into the same frame. The name
`reference_stimuli` removes the collision with the existing local `stimuli` (the noise matrix), but
not the other two ways the argument can go wrong:

- a same-named object in a file could override it;
- it could be written back into the user's file as if it were stimulus metadata.

Mitigation:
- restore the argument from `.args` after `load()`, as every `load()` site does;
- add it to `internals`;
- test that a file carrying a `reference_stimuli` object does not override the argument, and that a
  save adds only the cache field.

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
- **Guard:** the message fires on a subset CI scored without `reference_stimuli`, and on one scored
  with a different set. It stays silent for a full-set CI, when the sets match, and on a CI without
  the attribute. The attribute survives `generateCI2IFC()` and `batchGenerateCI()`, and does not
  change `names()` of the returned list.
