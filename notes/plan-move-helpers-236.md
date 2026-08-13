# Plan: move the internal helpers out of `zzz.R` (#236)

## Why

By R convention `zzz.R` holds load hooks (`.onLoad`/`.onAttach`) and package-level odds and
ends — it is named to sort last. This package's own convention is one function per file, named
after the function (`deg2rad.R`, `autoscale.R`, `generateCINoise.R`). `zzz.R` currently holds
four unrelated internal helpers and matches neither convention.

It is more scattered than `zzz.R` alone shows: `stopClusterSafely()`, the third parallel
helper and the only one called from *both* pipeline halves, sits at `R/generateCI.R:423` — a
571-line file that is not about parallelism.

## The moves

| function | from | to |
|---|---|---|
| `default_ncores()` | `zzz.R:9` | `R/parallel.R` |
| `startBackend()` | `zzz.R:34` | `R/parallel.R` |
| `progressOption()` | `zzz.R:50` | `R/parallel.R` |
| `stopClusterSafely()` | `generateCI.R:423` | `R/parallel.R` |
| `captureArgs()` | `zzz.R:69` | `R/rdata.R` |
| `rdataWriterNote()` | `zzz.R:88` | `R/rdata.R` |

`zzz.R` keeps `utils::globalVariables()` and nothing else — the one thing that genuinely
belongs there.

**Comments move verbatim with their functions.** They carry the reasoning (the `[8s/126s]`
measurement behind `registerDoSEQ()`, the connection-leak history behind
`stopClusterSafely()`, why `captureArgs()` skips required-and-absent arguments but never
defaulted ones). Rewriting them is not part of this change.

## Verified before writing this

- **No `Collate` field in `DESCRIPTION`**, so R sources `R/` alphabetically. All six are plain
  function definitions evaluated at load time, so order cannot matter. No new field needed.
- **All six are internal** — none appears in `NAMESPACE`, and none carries a roxygen block.
  So `NAMESPACE` does not change and no `man/` page appears or disappears.
- `R/` is not affected by `.Rbuildignore`; new files there ship, which is correct.

## What does not change

Nothing executable. This is cut-and-paste between files in the same package: same functions,
same bodies, same names, same (absent) exports. No numeric output, no behaviour, no argument.
The release gate should report `max|d| = 0` and the full suite should pass untouched.

## Scope boundary — what this does *not* take

`generateCI.R`'s `# Functions` section also holds `hasMask()`, `applyMask()`,
`applyScaling()`, `combine()` and `saveToImage()`. **They stay put.** They are `generateCI()`'s
own concerns and belong to #184 (splitting that file), which will want to group them by
computation/scaling/plotting rather than move them once now and again later. Only
`stopClusterSafely()` leaves, because it is the odd one out: parallelism, used by two files.

`#185` (de-duplicating mask-import logic) overlaps that same set and is a further reason to
leave the mask helpers alone here.

## Step most likely to fail

Nothing about the moves. The risk is **collision with in-flight work**:

- **#237** (#235's fix) also edited `generateCI.R`. It landed as `b2c8009` before this branch
  was written, and this branch is based on it, so the collision is already resolved rather
  than pending. Line numbers above are as of `b2c8009`.
- **#184** would move the same file wholesale. Sequence: this, then #184.

## Verification

Because "no behaviour change" is the entire claim, the check is that nothing observable moved:

1. `git diff -M --stat` should show the moved bodies as renames/moves, not rewrites — read the
   diff and confirm no line inside a moved function changed.
2. Full suite passes with **no test edited**. Any test needing a change would mean this was not
   a pure move.
3. `R CMD check` clean — in particular no new "no visible binding" NOTE, which is what a
   helper left behind or duplicated would produce.
4. `NAMESPACE` byte-identical after `roxygen2::roxygenise()`.

No `NEWS.md` entry: nothing user-visible changed, and `NEWS.md` is for what changed for users.

Docs, checked rather than assumed — only **one** of the four `zzz.R` mentions needs editing:

- `DECISIONS.md:381` — "`startBackend()` in `zzz.R`". **Needs the new path.**
- `AGENTS.md:183`, `CONTRIBUTING.md:224` — both about `globalVariables()`, which stays. No change.
- `CONTRIBUTING.md:194` — "`zzz.R` holds the `globalVariables()` declarations", given as the
  exception to one-file-per-function. Unchanged, and *more* true afterwards, since that becomes
  all it holds.
