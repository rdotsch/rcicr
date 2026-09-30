# Plan: generateStimuli2IFC() never overwrites stimulus PNGs (#350)

## What is wrong

#338 reserves the timestamped `.Rdata` name, but the PNG names (`<label>_<base>_<seed>_<trial>_ori.png`
and `_inv.png`, `R/generateStimuli2IFC.R:256-269`) have no time in them. A call in a later minute into
the same folder, with the same label, base label and seed, passes the `.Rdata` guard and overwrites the
earlier PNGs. The earlier `.Rdata` survives, and no longer describes them.

The maintainer confirmed both effects at runtime on #350. With the clock fixed to two minutes and
`nscales` changed:
- 8 of 12 PNG hashes changed, and trial 1's PNG matched the *second* file's noise;
- a shorter second run left trials 5 and 6 from the first run beside the new ones.

Both calls completed normally.

## Policy

The rule #338 set for the parameter file, extended to the PNGs: **never overwrite; stop before writing
anything**. The error says to use a different `label` or `stimulus_path`, or, to regenerate on purpose,
to delete the earlier stimulus set, its PNGs and its `.Rdata`.

No new argument. An `overwrite = TRUE` option would be a new contract that makes the bad state easy to
reach again, and deleting files is an explicit act the user already controls.

## Changes

1. **Preflight, before any generation or write.** Compute every PNG path the call will write: each
   base label and trial, `_ori` and `_inv`. If any exists, stop, releasing the `.Rdata` lock. The
   lookup is `file.exists()` on the real paths, so the file system decides what counts as the same
   name, including case and Unicode spellings. That is the concern that kept #338's lock off `label`.
   The message gives the number of existing files and the first few paths.
2. **A lock for concurrent runs.** A preflight alone still races: two same-seed calls into one folder
   in different minutes both pass it, then write together. With `save_as_png = TRUE`, the call also
   takes an atomic `dir.create()` lock keyed on the seed alone, not the minute and not the label, for
   the same reason as #338. It is held until the PNGs are written and released on error. A second
   same-seed call into the folder meanwhile stops with a message naming the lock, as #338's does.
   - **Cost:** two concurrent same-seed calls into one folder with *different* labels now also stop,
     though their files would not collide. That trade is deliberate: the lock cannot key on the label,
     and such calls are rare.
3. **Unchanged:** `save_as_png = FALSE` writes no PNGs and takes neither the preflight nor the new
   lock. The `.Rdata` reservation and its messages stay as #338 left them.

## Considered and not done

- **Time in the PNG names.** Experiment software refers to these files by name, so that would break
  every existing workflow.
- **Stopping on an earlier `.Rdata` with the same label and seed from another minute.** #338 allowed
  that on purpose, and matching it would need a pattern on `label`, which file systems do not spell
  consistently. The new error tells the user to delete the earlier `.Rdata` along with its PNGs.

## Tests, with the clock fixed through `stimulusTime()`

- The issue's reproduction: `nscales = 1` at 10:00, then `nscales = 2` at 10:01 into the same folder.
  The second call stops. Every PNG hash is unchanged, no second `.Rdata` exists, and no lock is left.
- A shorter second run (six trials, then four) stops the same way, leaving the first set intact.
- A different label, or a different `stimulus_path`, still succeeds.
- Regenerating on purpose, after deleting the earlier PNGs and `.Rdata`, succeeds.
- `save_as_png = FALSE` into a folder holding the PNGs succeeds, unchanged.
- With a planted seed lock, the call stops with the lock message, and nothing is written.
- Mutations that must fail a test: skipping the preflight; checking only the first base label; taking
  the lock after writing.

## Risk

The seed lock is the part most likely to surprise someone. A stale lock after a killed R session blocks
later same-seed calls into that folder until deleted, exactly as #338's lock does. The message says
how to check and delete it.
