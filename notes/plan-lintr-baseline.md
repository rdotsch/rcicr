# Plan — stop pinning lint exclusions to line numbers (#250)

## The problem, measured

`.lintr` records every accepted lint as a line number. Line numbers do not survive an edit above
them. On #249 this fired **twice on one branch**: adding four lines to `hasMask()` shifted seven
exclusions in `R/generateCI.R` and produced six failures in `applyScaling()`, ~90 lines below the
edit; adding one more comment line later did it again.

```
0 lints  ->  edit R/generateCI.R (+4 lines at ~419)  ->  6 lints at 510-535
```

The failure is silent until CI, it accuses untouched code, and the remedy —
`tools/regenerate-lintr-baseline.R` — rebaselines from scratch while ignoring the existing
exclusions, so a genuinely new lint introduced in the same commit is absorbed without anyone
seeing it.

## Two mechanisms, both verified before choosing

Neither was taken from documentation. Measured on scratch copies outside this repo so the
repo's own `.lintr` could not confound the result.

**`Inf` in place of the line vector** — "this linter, this whole file":

```
R/generateCI.R, no exclusions:              104 lints
                infix_spaces_linter = Inf:   62 lints
                  infix_spaces remaining:     0
                  other linters still firing: 5 categories
```

Granular per linter, not a blanket file skip, and immune to line movement.

**Inline `# nolint: <linter>.`** — per line, travels with the code:

```
t.R    line 2 (bare)          -> fires
       line 3 (# nolint:)     -> suppressed
t.Rmd  same, inside a chunk   -> same result   (so the vignette is covered)
t2.R   `y = 1:length(x)  # nolint: seq_linter.`
       -> assignment_linter still fires
```

The last line matters: a targeted `# nolint` suppresses only the named linter, so it is not a
blanket "ignore this line".

## The shape of the problem

300 lints across 30 files, split by what a lint would *mean* if a new one appeared:

| class | linters | lints | file x linter pairs |
|---|---|---|---|
| cosmetic | `infix_spaces`, `commas`, `indentation`, `assignment`, `pipe_consistency`, `spaces_inside`, `spaces_left_parentheses`, `brace` | 275 | 47 |
| semantic | `object_usage` (10), `commented_code` (7), `object_name` (4), `seq` (3), `object_length` (1) | 25 | 17 |

The 25 semantic ones are spread over 15 files, at most 4 per file.

## The design

**`Inf` for the cosmetic pairs; inline `# nolint:` for the 25 semantic lines.**

Not `Inf` everywhere, because `object_usage_linter` is the one here that can find a real bug —
an undefined or unused variable — and blanket-excluding it for `tests/testthat/test-legacy-rdata.R`
means the next genuinely undefined name in that file is never reported. 25 comments buys back
exact coverage on precisely the linters where "a new one appeared" is worth knowing about. Not
inline everywhere either: 300 comments would swamp the source, against `CONTRIBUTING.md` ->
Code conventions ("comment sparingly"), and rewrite 300 lines of `git blame`.

**The regenerator gets simpler, not more complex.** A line carrying `# nolint:` never appears in
`lint_package()` output, so once the semantic lints are annotated the tool has nothing but
cosmetic pairs left and can emit `Inf` unconditionally — no per-linter policy table inside it.
The policy lives in where the comment is, which is where someone reading the code will find it.

## The step most likely to fail

**The regenerator silently coarsening a future semantic lint.** With the tool emitting `Inf` for
whatever it is handed, a new `object_usage_linter` hit introduced next year gets blanket-excluded
for that whole file on the next regeneration — quietly undoing the distinction this plan draws,
and looking exactly like the routine renumbering diffs everyone has learned to rubber-stamp.

Mitigation, and it is part of the change rather than a note in a doc: the tool carries the list
of linters that must be annotated inline and **fails** when one reaches it, naming the file, the
line and the comment to add. That is the only place the classification is encoded, and it errors
rather than warns, because a warning in a script whose output is a 150-line diff is not read.

Second risk, smaller: `R CMD check` must stay clean. The `# nolint` comments land in `R/`,
`tests/` and one `.Rmd` chunk; they are ordinary comments, but the vignette one is inside code
that executes during the build.

## Rejected alternatives

- **Fix the lints and drop the baseline.** 275 of 300 are cosmetic and mechanically fixable.
  Rejected upstream already: `MAINTENANCE.md` records that `styler` is deliberately kept out of
  pre-commit because "it would reformat nearly every file in one sweep and destroy `git blame`".
  A targeted sweep is smaller but has the same character, and the guiding constraint means `R/`
  is not reformatted for taste.
- **Keep line numbers and make the regenerator run in a pre-commit hook**, so drift is repaired
  before it reaches CI. Rejected: `MAINTENANCE.md` explains why `lintr` is not a pre-commit hook
  — `object_usage_linter` needs the package *installed*, which the language-agnostic hooks never
  need. It also repairs the symptom on every commit rather than removing the fragility.

## What this cannot break

`.lintr` and `tools/` are outside the reproducibility gate's inert allowlist, so `compare` runs
in full. It will pass: nothing here touches `R/` behaviour beyond adding comments, and no numeric
output moves. The `# nolint` comments do change lines in `R/` — comment-only, but they will show
in `git blame` for those 25 lines, which is the price named above.

## Steps

1. Annotate the 25 semantic lints inline with `# nolint: <linter>.`
2. Rewrite `tools/regenerate-lintr-baseline.R` to emit `Inf` per file x linter, and to fail on
   any linter in the must-annotate list.
3. Regenerate `.lintr`; confirm `lint_package()` is 0 and that the tool is idempotent (rerunning
   leaves the file byte-identical — it was not, until #249 fixed a trailing newline).
4. Prove the fragility is gone: insert a throwaway line near the top of `R/generateCI.R`, confirm
   `lint_package()` stays 0, remove it.
5. Prove the guard still works: introduce one new cosmetic lint in an unbaselined file and one
   new `object_usage` lint, confirm each is reported.
6. Update the `lint.yaml` row in `MAINTENANCE.md`, and record there — and only there — why the
   cosmetic guard is deliberately coarser. Not `DECISIONS.md`: that file's test is "would this
   still matter if the package were maintained somewhere else entirely", and a lint baseline
   would not. Check `wc -w MAINTENANCE.md` against its 1800-word budget afterwards.
7. Delete this file, mark the draft ready, answer the Codex review, squash.
