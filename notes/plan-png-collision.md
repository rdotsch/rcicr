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

1. **Reserve every PNG before any generation.** Walk the PNG paths the call will write, in order:
   each base label and trial, `_ori` and `_inv`. For each one:
   - if it already exists, that is a collision, and the call stops;
   - otherwise, create it empty to reserve it.

   A file can exist either because an earlier run wrote it, or because this same call just reserved
   it under a spelling the file system treats as equal. For example, base labels `face` and `Face`
   on a case-insensitive or Unicode-normalising file system name the same file. Nothing is compared
   as strings, so what counts as the same name is always the file system's decision, which is the
   concern that kept #338's lock off `label`.

   The collision message gives the number of paths taken and the first few. It also says when the
   clash is inside the call itself, so it names equivalent base labels rather than an earlier run.
   On a collision, or on any error before the call completes, every file this call created is
   removed (placeholders and PNGs alike), and so is the `.Rdata` lock. A failed call therefore
   never leaves a partial set, or placeholders that would block the next run. The check goes through
   a small helper so that tests can model a case-insensitive file system on a case-sensitive CI
   runner.
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
4. **`NEWS.md`**, in the development version's "Bug fixes": calls that used to overwrite PNGs now stop,
   and how to regenerate on purpose.
5. **`DECISIONS.md`.** The entry "A stimulus `.Rdata` file is never overwritten" becomes "Stimulus
   files are never overwritten", and takes the policy and what was rejected:
   - time in the PNG names;
   - stopping on another minute's `.Rdata`;
   - an `overwrite` argument;
   - and the cost of keying the lock on the seed alone.

   The file is at its word budget, so the entry must fit by tightening the existing text.

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
- Base labels that alias each other within one call stop it, with nothing left behind:
  - a mocked case-insensitive file system (the existence helper), which runs everywhere;
  - `face` and `Face` for real where the file system is found to be case-insensitive, which the
    macOS runner normally is. The test checks for that by creating `a` and looking up `A`.
- An error part-way through generation leaves no placeholder and no PNG from the failed call.
- Mutations that must fail a test:
  - skipping the reservation;
  - checking only the first base label;
  - not removing reservations on failure;
  - taking the lock after writing.

## Risk

The seed lock is the part most likely to surprise someone. A stale lock after a killed R session blocks
later same-seed calls into that folder until deleted, exactly as #338's lock does. The message says
how to check and delete it.
