# Plan — cover the reference refresh in the release gate

## Why the gate does not reach it today

This PR makes `computeInfoVal2IFC()` regenerate a cached reference distribution that carries no
`reference_norms_source`, because nothing in such a file can show it was built on the noise the
stimuli actually use. That path is the riskiest thing here: it touches every existing file rather
than only the legacy ones the fix targets, and its rule changed three times under review.

The battery cannot see it. Every configuration writes its stimulus file with the version under
test, so on the current side the file is **fresh**: the first `computeInfoVal2IFC()` finds no
cache, regenerates, and writes the marker. There is never an unmarked cache, so the refresh never
runs. Scoring twice without more would only pin that cache *reuse* is stable — worth having, but
not a test of this change.

## The change

Add one output to the InfoVal configurations in `tools/compare-harness.R`, after `out$infoval`:

1. Remove `reference_norms_source` from the `.Rdata` if it is there. On v1.0.1 the field does not
   exist and this is a no-op, so both sides run the same code.
2. Score the same CI again with the same arguments, into `out$infoval_twice`.

On the reference side both calls hit the cache, so the value repeats. On the current side the
second call takes the refresh path and rebuilds from the saved noise. The gate then asserts,
across versions, that the refreshed value equals what v1.0.1 computed — which is the property
this PR claims and the only place it is checked against the real old code rather than against an
oracle this repo wrote for itself.

This is an added output on four existing configurations, not a new configuration, so nothing is
asked of the reference version that it does not already do: `computeInfoVal2IFC()` twice, and a
`load()`/`save()` round trip between. `AGENTS.md`'s warning about the reference version being
unable to execute a configuration does not apply, and it will be confirmed by running the gate
rather than assumed.

## Expected impact

- **No new `EXPECTED` entry.** `sinusoid-64-infoval` and `defaults-*` have no InfoVal deviation
  today, so their `infoval_twice` must match too. Three of the four InfoVal configurations do
  deviate — `sinusoid-64-nscales3-infoval`, `gabor-64-sigma10-infoval` and
  `sinusoid-64-twobase-indep-infoval` — and their `infoval_twice` will deviate for the same
  reason and needs the same entry, extended by key rather than added as a new cause.
- **No change to the package.** This is gate machinery only; `R/` is untouched.
- A deviation on `sinusoid-64-infoval/infoval_twice` specifically would mean the refresh does
  not reproduce what it replaced, which is the failure this exists to catch.

## Most likely to fail

The `EXPECTED` bookkeeping rather than the property. Each existing InfoVal entry is keyed to one
output, and three of them now have a second output with the same cause; a key that stops firing
fails the run as stale, so the keys have to be extended in the same commit. `tools/compare-release-output.R`
accepts a vector of keys for one cause, which is what that is for.

## Verification

`Rscript tools/compare-release-output.R --quick` runs here — it needs the release tags fetched
and `assertthat` present for the reference version — and currently reports
`PASS (quick): this tree reproduces v1.0.1 on the configs that were run`, with 15 expected
deviations on file. The same run has to pass with the new output, and the `NOT EXERCISED` list
must not grow.
