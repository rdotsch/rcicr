# Plan: never overwrite a stimulus `.Rdata` file (#338)

## Problem

`generateStimuli2IFC()` names its file `<label>_seed_<seed>_time_<%b_%d_%Y_%H_%M>.Rdata`, a one-minute resolution, and `save()`s over whatever has that name. Two calls into one `stimulus_path` in the same minute, with the same `label` and `seed`, write one file: the second replaces the first. Measured on `main`, with base images named `male` then `female`: 16 PNGs, one `.Rdata`, holding only `female`. The `male` PNGs have lost the only record of their noise parameters.

The name is built at the end, after every PNG is written, so a check there would come too late to prevent anything.

## Change

1. **Reserve the name before any stimulus is generated, with a lock beside the file, never the file itself.** The lock is a directory named like the file with `.Rdata` replaced by `.rcicr-lock`, created with `dir.create()`, which is atomic and fails if it exists. Holding it, the call checks that the `.Rdata` file itself does not exist yet, generates, `save()`s the real file at the end, and removes the lock (also on error). Two calls that overlap cannot both hold the lock, however long generation takes; a check-then-write would let both through and the later `save()` would still overwrite the earlier. If the lock cannot be taken or the file exists, the call stops before writing anything, with an error that names the file and says to use a different `label` or `stimulus_path`, or to delete the file if it is a set being regenerated on purpose. A call killed outright leaves the lock, which then blocks that name. The error names the lock and says it belongs either to a call still running or to one that was interrupted, and to delete it only after confirming that no `generateStimuli2IFC()` call is still writing to that folder. rcicr cannot tell the two apart: a process ID proves nothing across machines or network drives, so the lock records no owner and the error claims none. The lock never contains "Rdata", so no `list.files()` pattern for the output, with or without `$`, can pick it up; the empty placeholder that reserving the file itself would leave could be.
2. **The time in the name becomes the start of the call** instead of the end. Nothing reads that time; it only makes names unique, and it is now fixed before anything is written.
3. **Only when `save_rdata = TRUE`.** Without an `.Rdata` file there is nothing to lose.
4. **`NEWS.md`, under "Bug fixes"**: what used to happen, what happens now, and that the time in the name is the start of the call.
5. **`?generateStimuli2IFC`**: the file name, and that an existing one stops the call.

**Rejected: adding a suffix (`_2`) on collision.** It loses nothing, but a script that picks the file with `list.files(pattern = "Rdata$")[1]` would then silently get the *older* file, since `..._13_24.Rdata` sorts before `..._13_24_2.Rdata`. With identically named base images that older file no longer matches the PNGs on disk. A stop is loud; a suffix can hand a CI the wrong parameters. Adding seconds to the name only makes a collision rarer.

## Verified before planning

- The reproduction above, on `main`.
- `dir.create()` on an existing directory returns `FALSE` (measured: `TRUE`, then `FALSE`), as `?dir.create` documents for a directory that already exists; the underlying `mkdir` is a single atomic operation on every OS CI checks.
- Every example in the repository finds the output with `list.files(stimulus_path, pattern = "\\.Rdata$")[1]`; a `.rcicr-lock` name matches neither that nor a bare `"Rdata"`.
- Every caller in the repository tolerates the stop: with `generateStimuli2IFC()` temporarily made to stop whenever its `.Rdata` name already existed, `testthat::test_local()` had no failures. `tools/compare-harness.R` empties its directory before each call.

## Tests

- two calls with the same `label` and `seed` into one directory, with the clock mocked to one minute: the second stops, and the directory holds exactly the first call's PNGs and `.Rdata`, unchanged (checksums);
- a different `label`, or `save_rdata = FALSE`, does not stop;
- a held lock stops a second call before it writes anything; its error names the lock and says to delete it only once no call is writing to the folder;
- a call that errors after taking the lock leaves neither the lock nor an `.Rdata` file;
- while a call runs, `list.files(pattern = "Rdata")` finds nothing new (checked from inside the generation step with a mocked hook);
- the name carries the start time: with the clock mocked to advance during the call, the file is named for the first reading.

The failure test must fail on the current code; checked with `git stash push -- R/`.

## The step most likely to fail

**Quick reruns.** Rerunning an identical script into the same folder within the same minute now stops, where it used to overwrite the file with identical contents. That is the cost of the stop, and the message says what to do. The suite measures how common it is in practice: nothing in it does this.

## Out of scope

- **PNG overwrites.** PNG names carry no time, so a later call with the same `label`, base names and `seed` rewrites the PNGs whatever the minute, even when the `.Rdata` files differ. That is also how regenerating a set into its own folder works, which this change should not break; it is a separate question.
- #339.
