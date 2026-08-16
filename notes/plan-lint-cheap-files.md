# Plan — fix the cosmetic lints in the 14 cheap files (#250 follow-up)

## What this changes

`.lintr` currently excludes cosmetic linters per file with `Inf` for all 20 baselined files.
That is coarser than it needs to be for most of them: **10 of the 20 need 3 lines or fewer
fixed**, and the whole cheap half is 45 lines across 14 files.

Fix those 45 lines. A file with no lints leaves `exclusions:` entirely and regains **full**
strictness on every linter — better than annotating, which leaves the linter off for that line
forever and adds a comment that serves the linter rather than the reader.

The 6 dense files keep `Inf`: `generateCI.R` (32 lines), `batchGenerateCI.R` (19, deprecated),
`computeInfoVal2IFC.R` (16), `generateNoisePattern.R` (14), `generateStimuli2IFC.R` (13),
`generateReferenceDistribution.R` (11) — 105 lines whose style is frozen and which nobody should
be reformatting for taste.

Nothing is annotated with `# nolint` unless a fix turns out not to be safe. That is the
preferred outcome, not a shortcut: no comment is better than a comment.

## What is being fixed, and how risky each kind is

45 lines, 76 lint instances:

| linter | n | change | risk |
|---|---|---|---|
| `infix_spaces` | 41 | `a=b` -> `a = b` | whitespace |
| `commas` | 15 | `f(a,b)` -> `f(a, b)` | whitespace |
| `indentation` | 12 | leading spaces | whitespace |
| `spaces_inside` | 3 | `f( x )` -> `f(x)` | whitespace |
| `spaces_left_parentheses` | 2 | `if(x)` -> `if (x)` | whitespace |
| `brace` | 1 | `} else{` -> `} else {` | whitespace |
| `assignment` | 2 | `=` -> `<-` | **not whitespace** |

The two non-whitespace ones, both plain statement-level assignments to a local, where `<-` is
exactly equivalent:

```
R/generateGabor.R:19    gauss_mask = exp( -(((gauss$x^2)+(gauss$y^2)) / ...
R/generateSinusoid.R:16 sinepatch = matlab::repmat(matlab::linspace(0, cycles, img_size), ...
```

`pipe_consistency_linter` would be the dangerous one — `%>%` and `|>` are not interchangeable in
every position — and it has **no hits in the cheap files**, only in `computeInfoVal2IFC.R`, which
stays on `Inf`. Nothing in scope here requires judgement about semantics.

## The step most likely to fail

**A whitespace edit landing inside a string literal.** `git diff -w` cannot catch it: it ignores
whitespace everywhere, string contents included, so the check that looks most obviously right is
precisely the one that is blind to the only silent failure available here. A changed string would
alter an error message, or a filename, with a clean-looking diff.

The check that does work is to compare the **parsed** form, which discards whitespace outside
strings and preserves it inside:

```r
identical(parse(f, keep.source = FALSE), <same file before>)
```

Expected result: identical for all 14 files except `generateGabor.R` and `generateSinusoid.R`,
which must differ in exactly one expression each — `=` and `<-` parse to different calls. Any
other difference is a real change and must be explained before merging.

## Verification

1. Parse comparison above: 12 files identical, 2 differing in exactly the assignment operator.
2. Full test suite, including `test-regression-baseline.R`.
3. The reproducibility gate runs by itself — `R/` is outside the inert allowlist — and must be
   green with no `EXPECTED` entry added. If it deviates, this change is wrong; a formatting fix
   cannot move a number, so a deviation means an edit was not a formatting fix.
4. `lint_package()` is 0, and the 14 files are gone from `.lintr` rather than listed with fewer
   entries.
5. Re-run the drift check from #250: insert lines at the top of a fixed file, confirm 0 lints.

## No `NEWS.md` entry

Nothing changed for users: no signature, no message, no number. `NEWS.md` is what changed for
them, and a whitespace pass is not it.

## `MAINTENANCE.md` is at 1799 of its 1800 words

The paragraph describing the baseline says cosmetic linters get `Inf`; it now needs to say that
this applies to the dense legacy files, the rest being fixed outright. There is no headroom, so
something comes out before this goes in — the edit has to be net-neutral or shorter.

## Steps

1. Fix the 45 lines in the 14 cheap files.
2. Run the parse comparison; investigate anything beyond the 2 expected differences.
3. Regenerate `.lintr`; confirm the 14 files are absent and `lint_package()` is 0.
4. Full suite; check the gate in CI.
5. Update the `MAINTENANCE.md` paragraph within budget (`wc -w`).
6. Delete this file, mark the draft ready, answer the Codex review, squash.
