# Plan: enforce the documentation word budgets in CI (#287)

## What is measured today

`LC_ALL=C.UTF-8 wc -w` on `main` at 7e65245:

| file | words | budget | spare |
|---|---|---|---|
| `AGENTS.md` | 2655 | 2800 | 145 |
| `CONTRIBUTING.md` | 2627 | 2800 | 173 |
| `RELEASING.md` | 1530 | 1600 | 70 |
| `MAINTENANCE.md` | 1695 | 1800 | 105 |
| `SECURITY.md` | 325 | 600 | 275 |
| `DECISIONS.md` | 5200 | 5200 | **0** |

The budgets are stated twice: in the table in `AGENTS.md` → "Which file a thing goes in", and in a "Keep this file under N words" line at the top of each file except `AGENTS.md`, whose line is in its first paragraph. Today all of those agree.

## The four decisions the issue names

- **Where the numbers live: the `AGENTS.md` table**, parsed. It is the one place every budget is already listed, and duplicating them in the script would add a third copy to drift. Rows are matched by their shape, `` | `FILE.md` | ... | N | ``. A row whose last cell is not a number (`NEWS.md`'s "none; never trimmed") is skipped. A file named in the table but missing, or a table yielding no rows at all, fails the check, so a reformatted table cannot silently pass by matching nothing.
- **The per-file "Keep this file under N words" lines are checked against the table too.** They are the copy a reader of that file sees, and they are what drifts if only the table is edited.
- **Inclusive: a file fails only when it is over its budget.** `AGENTS.md` and `DECISIONS.md` both phrase the rule as "over budget, something comes out", and `DECISIONS.md` sits at exactly 5200 today under that reading. The "under N" wording stays as it is.
- **A failure, not a warning.** A warning in a passing check is not read; the budgets exist because a doc that only grows stops being read. The error names the file, its count, its budget and the excess, so the fix starts from a number.

Only the table's rows are checked. `README.md`, `NEWS.md`, `notes/` and `analyses/` have no budget.

## Shape

- `tools/check-word-budgets.sh`, runnable locally from the repository root, so the number a contributor sees is the one CI enforces. It counts with `LC_ALL=C.UTF-8 wc -w`, the command `AGENTS.md` already gives, rather than reimplementing word splitting in R, so the two can never disagree. The locale is set inside the script, so a contributor's own locale cannot change the count.
- A last step in `ubuntu-latest (release)`, after "CITATION.cff matches its sources". It is not a new job, for the ruleset reason `MAINTENANCE.md` records: required checks match by name, and a new job's check would report without blocking. It goes last so a failure cannot hide the check results.
- `MAINTENANCE.md` gets one short entry next to the stale-`man/` gate, and `AGENTS.md`'s instruction to run `wc -w` points to the script. Both stay within their budgets, which the script then checks.

## Verification

- Run against the current tree: passes and prints the table above.
- Mutations, each expected to fail with a message naming the cause: add 1 word to `DECISIONS.md`; change a budget in the table without changing the file's own line; change a file's line without the table; rename a file in the table to one that does not exist; break the table's format so no row matches.
- In CI, the step's output is quoted in this PR once it runs.

## The step most likely to fail

Parsing the markdown table. It is the reason the issue offered duplicating the numbers instead, and an edit to the table's layout is the likeliest way to break the check. A layout it does not recognise therefore fails the check with the instruction to update the script, never passes it.
