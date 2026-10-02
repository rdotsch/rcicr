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

- **Where the numbers live: the `AGENTS.md` table**, parsed. It is the one place every budget is already listed, and duplicating them in the script would add a third copy to drift. The table is the block of `|` lines under the `| file | its job | budget |` header, and **every row in it must parse**: a first cell naming a backticked file, and a last cell that is either a whole number or begins with `none` (`NEWS.md`'s "none; never trimmed", the only unbudgeted row). Any other row fails the check, naming the row, so a reformatted or garbled row cannot drop out unnoticed. So do a missing header, a file named in the table but absent, and a table with no budgeted row.
- **A deleted row is caught from the other side.** Every top-level `*.md` file that carries a "Keep this file under N words" line must have a budgeted row in the table. Deleting a row while the file keeps its line fails the check, and so does a new doc given a line but no row.
- **The per-file "Keep this file under N words" lines are checked against the table too.** They are the copy a reader of that file sees, and they are what drifts if only the table is edited.
- **Inclusive: a file fails only when it is over its budget.** `AGENTS.md` and `DECISIONS.md` both phrase the rule as "over budget, something comes out", and `DECISIONS.md` sits at exactly 5200 today under that reading. The "under N" wording stays as it is.
- **A failure, not a warning.** A warning in a passing check is not read; the budgets exist because a doc that only grows stops being read. The error names the file, its count, its budget and the excess, so the fix starts from a number.

Only the table's rows are checked. `README.md`, `NEWS.md`, `notes/` and `analyses/` have no budget.

## One budget changes: `DECISIONS.md` goes from 5200 to 6000

Decided with the maintainer. It is the only doc with no room (5200 of 5200), and the only one that grows by design: `AGENTS.md` asks for an entry for every decision that was not obvious. It is read by section, through its headings and links, not end to end. At its budget, each new entry has to displace older reasoning; the last one cost about 45 words of measured detail. The process docs stay as they are: they describe how to work, should stay short, and have 70 to 275 words spare. The table in `AGENTS.md` and the line at the top of `DECISIONS.md` both change to 6000 in this PR, which the check itself then confirms agree.

## Shape

- `tools/check-word-budgets.sh`, runnable locally from the repository root, so the number a contributor sees is the one CI enforces. It counts with `LC_ALL=C.UTF-8 wc -w`, the command `AGENTS.md` already gives, rather than reimplementing word splitting in R, so the two can never disagree. The locale is set inside the script, so a contributor's own locale cannot change the count.
- A last step in `ubuntu-latest (release)`, after "CITATION.cff matches its sources". It is not a new job, for the ruleset reason `MAINTENANCE.md` records: required checks match by name, and a new job's check would report without blocking. It goes last so a failure cannot hide the check results.
- `MAINTENANCE.md` gets one short entry next to the stale-`man/` gate, and `AGENTS.md`'s instruction to run `wc -w` points to the script. Both stay within their budgets, which the script then checks.

## Verification

- Run against the current tree: passes and prints the table above.
- Mutations, each expected to fail with a message naming the cause: add 1 word to `DECISIONS.md`; change a budget in the table without changing the file's own line; change a file's line without the table; rename a file in the table to one that does not exist; delete one budgeted row; make one row's budget non-numeric; drop one row's leading `|`; remove the header so no table is found.
- In CI, the step's output is quoted in this PR once it runs.

## The step most likely to fail

Parsing the markdown table. It is the reason the issue offered duplicating the numbers instead, and an edit to the table's layout is the likeliest way to break the check. A layout it does not recognise therefore fails the check with the instruction to update the script, never passes it.
