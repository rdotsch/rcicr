# Plan — fix the remaining lints in the 6 dense files (#250 follow-up 2)

## Why, given the last round stopped short of these

The last round left these six on a whole-file `<linter> = Inf` and justified it as "dense legacy
files nobody edits". **That justification was wrong, and checking it is what prompted this
round.** Every one has been edited in the last 90 days:

| file | commits total | in last 90d | last touched |
|---|---|---|---|
| `generateCI.R` | 41 | 24 | 2026-08-16 |
| `generateStimuli2IFC.R` | 36 | 19 | 2026-08-13 |
| `computeInfoVal2IFC.R` | 25 | 14 | 2026-08-15 |
| `generateReferenceDistribution.R` | 19 | 13 | 2026-08-16 |
| `batchGenerateCI.R` | 8 | 4 | 2026-08-07 |
| `generateNoisePattern.R` | 3 | 2 | 2026-08-16 |

Lint density tracks file size and age, not inactivity. These are the *most* active files in the
package, so switching their cosmetic linters off wholesale left the least guarded code exactly
where new code is most likely to land. `generateCI.R` was edited twice on the day it was
exempted.

## Scope

199 lints on 105 lines. 186 are whitespace; two groups need judgement.

| linter | n | treatment |
|---|---|---|
| `infix_spaces` | 144 | whitespace |
| `commas` | 18 | whitespace |
| `indentation` | 17 | whitespace; the hanging-indent ones are fixed by moving `)` to its own line, as last round |
| `assignment` | 7 | `=` -> `<-`, all in `generateNoisePattern.R` |
| `pipe_consistency` | 6 | **configure, do not rewrite** |
| `spaces_left_parentheses` | 4 | whitespace |
| `spaces_inside` | 3 | whitespace |

**The assignments: 12, not the 7 the linter reports.** Raised in review, and confirmed wider
than the finding. `assignment_linter` has blind spots — it misses the four inside the
`if (pre_0.3.0)` branches at `generateNoisePattern.R:32-36`, and one at `generateCI.R:255`
(`params = -params`). Counted from the AST rather than from the linter or a regex, by walking
the parse data for `EQ_ASSIGN`:

```
R/generateCI.R            1   (line 255)
R/generateNoisePattern.R 11   (lines 21, 24, 27, 28, 32, 33, 35, 36, 59, 62, 68)
total                    12
```

All 12 are plain statement-level, including one indexed assignment (`patchIdx[...] = idx`),
where `<-` is exactly equivalent. **Fix all 12**, not the 7 that lint: a zero-lint result would
otherwise certify a file that still violates the convention.

**And the convention's own claim is false.** `CONTRIBUTING.md:192` reads "Assignment | `<-`,
never `=` | Already universal — there are zero `=` assignments in `R/`." There are 12. Fixing
them all makes the sentence true, which is better than rewording it to admit exceptions.

**The 6 pipes** are all in `computeInfoVal2IFC.R`, in the `ref_lookup` region that
`AGENTS.md` flags as delicate. Rewriting them is the wrong fix: the package uses `%>%` 7 times
and `|>` zero times, so it is already *consistent*, which is what the linter is named for — it
merely defaults to preferring native. `%>%` comes from `dplyr` (`importFrom(dplyr,"%>%")`), which
is an Import regardless for `count`/`filter`/`summarise`, so converting would not shed a
dependency either. Set `pipe_consistency_linter(pipe = "%>%")` in `.lintr`'s `linters:` block,
which is how that block already encodes house style for `object_name_linter`. Verified: it
clears all 6, and still catches a stray `|>` —

```
a stray |> under pipe="%>%": 1 lint -- Use the %>% pipe operator instead of the |> pipe operator.
```

## The step most likely to fail

**`.lintr` ends up with an empty `exclusions:` block, and the obvious empty form breaks lintr.**
If all six files come clean there is nothing left to baseline, and
`tools/regenerate-lintr-baseline.R` builds the block by pasting file entries — with none, it
emits `exclusions: list(\n\n  )`, which does not parse: `lint_package()` aborts rather than
returning zero. A `.lintr` that errors is worse than a stale one, and CI would report it as a
lint failure with a message about the config rather than the code.

Measured, on a tree with 199 lints so a working config must report exactly that:

```
exclusions: list(\n\n  )   ->  ERROR (rlang::abort)
exclusions: list()         ->  199
key omitted entirely       ->  199
```

So the regenerator needs a no-lints branch emitting `exclusions: list()` on one line. This is
worth handling even if some file keeps an entry today, because the next cleanup would hit it.

## Verification

1. Parse-tree comparison per file, as last round — `git diff -w` is not evidence here because it
   ignores whitespace inside string literals too, and these files are full of `stop()`/`paste0()`
   message text. Expected: 4 of 6 identical, `generateNoisePattern.R` differing in exactly 11
   deparsed lines and `generateCI.R` in 1, all `=` -> `<-`.
2. Full suite including `test-regression-baseline.R`.
3. The reproducibility gate runs on its own (`R/` is not inert) and must be green with no
   `EXPECTED` entry.
4. `lint_package()` is 0, and `.lintr` parses — check by running the linter, not by reading it.
   Separately, the AST walk for `EQ_ASSIGN` over `R/` returns 0, which the linter cannot
   establish since it cannot see 5 of the 12.
5. Drift check: insert lines at the top of `generateCI.R`, confirm still 0.

## `MAINTENANCE.md`

Its baseline paragraph says six dense files keep `Inf`; that becomes false. The file is at 1799
of its 1800-word budget, so the paragraph is rewritten rather than extended.

## Steps

1. Add `pipe_consistency_linter(pipe = "%>%")` to `.lintr`'s `linters:` block.
2. Fix the 105 lines across the 6 files.
3. Teach the regenerator the empty case; regenerate; confirm `.lintr` parses and is 0.
4. Parse comparison, full suite, gate, drift check.
5. Rewrite the `MAINTENANCE.md` paragraph within budget.
6. Delete this file, mark the draft ready, answer the Codex review, squash.
