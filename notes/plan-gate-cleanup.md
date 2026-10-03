# Plan: make the release gate clean up after itself (#401)

## The defect

`tools/compare-release-output.R` registers `cleanup()` with a top-level `on.exit()`, which never runs under `Rscript` (shown on #401). Every run without `--keep` leaves its temporary directory (up to 872 MB measured) and a `git worktree` registration behind. In this session that filled a 20 GB disk and made the next gate run fail.

## The change

1. Replace `on.exit(cleanup(), add = TRUE)` with `reg.finalizer(environment(), function(e) cleanup(), onexit = TRUE)`. A finalizer with `onexit = TRUE` runs when R exits: after the last line, on `quit(status = 1L)` or `die()`'s `quit(status = 2L)`, and after an uncaught error (all three shown on #401).
2. `cleanup()` stays as it is: `--keep` still keeps the directory and says where it is; otherwise the worktree is removed and the directory unlinked.
3. A finalizer can run more than once if the environment is collected before exit. The global environment never is, so it runs once. `cleanup()` is idempotent anyway: removing a removed worktree and unlinking a missing directory are both no-ops.

## Verification, and the step most likely to fail

- **A run that passes:** `--quick --ref=<main>` leaves no `rcicr-compare-<pid>` directory, and `git worktree list` gains no entry.
- **A run that dies early:** `--ref=not-a-ref` (`die()`) leaves neither.
- **A run that fails:** a gate run that ends in `quit(status = 1L)`, with `--keep` *not* given, leaves neither. With `--keep`, the directory is still there and printed.
- **Most likely to fail:** the finalizer running before the driver has finished reading the outputs, if R collected the environment early. It is the global environment, which is never collected, and the run checks `PASS`/`FAIL` is still printed in full before the directory goes.
