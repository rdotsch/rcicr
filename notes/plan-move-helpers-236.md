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
The full suite should pass untouched, and the release gate should **pass unchanged** — which
is not the same as `max|d| = 0`. The gate compares the working tree against the pinned v1.0.1
(and v1.1.0) references, which already deviate deliberately: `tools/compare-release-output.R`
carries `EXPECTED` entries for the InfoVal cases against v1.0.1 and the blur-based z-maps
against v1.1.0. Those keep firing after a pure move, and must. The criterion is the gate's own:
every `EXPECTED` entry still fires, none goes stale, and nothing new deviates.

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

1. **Each helper is byte-identical before and after, compared in order.** Two approaches were
   tried and discarded first, both of which would have reported success regardless of what the
   edit did:

   - `git diff -M` detects file *renames*. Both source files stay, so a split shows as two
     modifications plus two additions and says nothing about the bodies.
   - A whole-patch line-multiset comparison fails on this very branch — `git add -A` also
     stages the plan's own deletion and the two pointer rewrites, whose lines have no
     counterpart, so it reports a difference on a *correct* move. And sorting discards order:
     swapping `makeCluster()` and `registerDoSNOW()` inside `startBackend()` — which would
     register a backend on a cluster that does not exist yet — passes it. Verified, not
     supposed.

   So compare per function, in order, scoped to the helper itself. Extract the definition
   with the comment block directly above it:

   ```sh
   extract() {  # extract <fname> < file.R
     awk -v fn="$1" '{ line[NR] = $0 }
       $0 ~ "^"fn" <- function" { start = NR }
       END { first = start; while (first > 1 && line[first-1] ~ /^#/) first--
             last = start; while (last < NR && line[last] !~ /^\}/) last++
             for (i = first; i <= last; i++) print line[i] }'
   }
   base=$(git merge-base main HEAD)
   diff <(git show $base:R/zzz.R | extract startBackend) \
        <(extract startBackend < R/parallel.R)
   ```

   Empty for all six is the invariant. Being scoped to the function, it is unaffected by the
   plan deletion and the pointer rewrites; being a plain `diff`, it catches a reorder. Both
   properties were checked against the real `startBackend()` before this was written down.

2. **Each helper exists exactly once in `R/` afterwards** — `grep -rc '^<fn> <- function' R/`
   summing to 1 per function. That is what a copy-instead-of-move, or a botched deletion,
   would break, and the per-function diff above cannot see it.
3. Full suite passes with **no test edited**. Any test needing a change would mean this was not
   a pure move.
4. `R CMD check` clean — in particular no new "no visible binding" NOTE, which is what a
   helper left behind or duplicated would produce.
5. `NAMESPACE` byte-identical after `roxygen2::roxygenise()`.

No `NEWS.md` entry: nothing user-visible changed, and `NEWS.md` is for what changed for users.

Pointers to fix, from a **repo-wide** grep (`grep -rn 'zzz\.R' . --exclude-dir=.git`) — an
earlier docs-only grep missed the code comment and reported one fewer:

- `DECISIONS.md:381` — "`startBackend()` in `zzz.R`". **Needs the new path.**
- `R/generateStimuli2IFC.R:185` — "See `startBackend()` in `zzz.R`". **Needs the new path.**
  A stale pointer in a comment is worse than in prose: it is read by whoever is editing the
  loop it annotates.
- `AGENTS.md:183`, `CONTRIBUTING.md:224` — both about `globalVariables()`, which stays. No change.
- `CONTRIBUTING.md:194` — "`zzz.R` holds the `globalVariables()` declarations", given as the
  exception to one-file-per-function. Unchanged, and *more* true afterwards, since that becomes
  all it holds.
