# Plan: an unwritable archive must not turn a computable InfoVal into an error

## What was measured

R 4.3.3 in this container (`tools/setup-container-r.sh`), rcicr 1.4.0 installed from the
PR #321 head (4d844ff). Sessions run as root, which writes through read-only permission
bits, so every probe below was re-run under `setpriv --reuid=65534` on a `0444` archive.

`resolveReferenceNorms()` is reached only from `computeInfoVal2IFC()` — directly for shared
parameters, via `computeBaseInfoVal()` for independent ones. It sets

```r
readonly   <- automatic && !writableFile(rdata)
save_rdata <- is.null(response_seed) && !readonly
```

so the read-only fallback exists **only** for `automatic`, the unsolicited refresh of a stale
entry. Every other path asks `generateReferenceDistribution2IFC()` to save. As `nobody` on a
`0444` file, all three fail after the reference has already been simulated:

| case | entry | result |
|---|---|---|
| independent base, legacy unscoped `reference_norms`, no `reference_norms_by_base` | `NULL` | `ERROR: cannot open the connection` |
| shared base, no cache at all | `NULL` | `ERROR: cannot open the connection` |
| shared base, valid cache, `force_gen_ref_dist = TRUE` | present | `ERROR: cannot open the connection` |

The reported case is the first row; it is one instance of a general gap, not a per-base one.
Under the mocked `writableFile()` seam the tests use, the same call instead **writes the new
per-base cache into the archive** (md5 changes) — root defeats the permission bits, so the
mocked seam and real bits show the two halves of one defect.

The first row is also a regression for archived files: before per-base caches, an
independent-base file with an unscoped `reference_norms` scored from that cache and wrote
nothing.

## The change

In `resolveReferenceNorms()`, decouple the fallback from `automatic`:

```r
wanted_save <- is.null(response_seed)
readonly    <- wanted_save && !writableFile(rdata)
save_rdata  <- wanted_save && !readonly
```

`wanted_save` keeps a `response_seed` draw on its own "deliberately not saved" message rather
than reporting writability, which is irrelevant when no save was ever intended.

`automatic` keeps its other two jobs untouched — inheriting the cached `iter` and preserving
the random stream — so no simulated draw and no returned number moves. The change converts an
error into the value that was already computed.

The line this rests on: in `computeInfoVal2IFC()` the cache is an optimization, so an
unwritable archive must not cost the caller their InfoVal. A direct
`generateReferenceDistribution2IFC(save_rdata = TRUE)` still errors, because there the save
*is* the request.

The read-only message drops its `automatic`-specific "Rebuilt"/"rebuilt" wording for one
phrasing that fits a first build too. The separate "supersedes" message still marks a refresh
whose values moved.

## Most likely to fail

Tests that call `resolveReferenceNorms()` with a placeholder path (`"unused"`,
`"archive.Rdata"`). Those files do not exist, so `file.access()` returns `-1` and the widened
`readonly` would fire where `automatic` used to suppress it. They appear to mock
`writableFile()` already; the full suite decides, not this reading.

## Verification

1. Three regression tests — the reported independent case, the shared-base first build, and
   forced regeneration — each denying writes through the existing `deny_writes()` helper, each
   asserting a finite InfoVal, an unchanged md5, and the value a writable run returns. Run
   against the unfixed code first: they must fail.
2. `testthat::test_local()` in full.
3. `R CMD check`.
4. The release gate against `origin/main` per `CLAUDE.md`, expecting `0 expected deviations`,
   to show the change moves no number.
5. Package code changes, so the 1.4.0 tarball and both win-builder submissions are void and
   must be rebuilt from the new head.
