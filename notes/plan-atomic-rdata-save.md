# Plan: save reference distributions without risking the stimulus file (#333)

## Problem

Two `save()` calls overwrite the stimulus `.Rdata` file in place: `R/generateReferenceDistribution.R:174` (shared parameters) and `R/reference-base.R:175` (independent base images). `save()` truncates the file before writing, so an interrupted write destroys it. #333 has the reproduction: killed mid-write, a valid 733-byte file became 10,993,664 bytes that `load()` rejects.

## Approach: back up, then save in place

The file keeps being written in place, exactly as today, so its owner, group, ACLs, mode, hard links and symlinks are untouched. What changes is that a complete copy exists for the duration of the write.

Writing a temporary file and renaming it over the original was rejected: a rename replaces the file, so it takes the process's owner and group (measured: `nobody:nogroup` became `root:root`), drops ACLs, and base R can restore neither. For a shared archive that silently changes who can read the stimulus file.

## Change

1. **One internal helper, `saveRdataSafely(names, file, envir)`, in `R/rdata.R`**, used by both call sites. It runs in four phases and records which one it is in, because what an interruption must undo differs per phase. `backup` is `paste0(file, ".rcicr-backup")`.

   **Phase 0, checks (nothing written yet).**
   - List leftover `<file>.rcicr-staging-*` files and warn with their names, saying they are incomplete copies that can be deleted. They are not deleted: a file this call did not create could belong to another R session saving at the same moment.
   - **If `backup` exists, validate the original by loading it** into a fresh environment. If it loads, it is a complete file (from before the killed save, or written by it just before the kill), so the backup is stale: remove it and say so. If it does not load, stop, naming both files and the restore command below. Restoring without this check could put an older reference over a newer, complete one.
   - Generate the staging path with `tempfile(pattern = paste0(basename(file), ".rcicr-staging-"), tmpdir = dirname(file))`, which names a file without creating it. **If that name exceeds 255 bytes** (`nchar(basename(staging), type = "bytes")`; the staging suffix is 27 bytes, measured, so this means an original name over 228 bytes, where rcicr's own names are about 40), save in place without a backup and warn that the name is too long to protect. This check comes before any file is created, so the name limit cannot surface as a creation failure.

   **Phase 1, backup.** An interruption here removes this call's staging file and restores nothing: the original has not been touched.
   - Create the staging file empty. If that fails, check the directory with `file.access(dirname(file), 2)`: only on Unix-alikes, and only when it reports the directory not writable (a read-only directory holding a writable file), save in place without a backup, as today. **Any other failure stops**: a directory that reports writable but refuses a new file is out of space, inodes or quota, where saving in place is the very failure this change prevents.
   - Set it to mode 600 with `Sys.chmod(staging, "600", use_umask = FALSE)`. On Unix-alikes, check that `file.mode(staging)` is then 600 and stop if not: `Sys.chmod()` reports failure by its return value, and a staging file left at its creation mode would expose the data after a hard kill. Mode 600 lets only the user running R read it, who could already read the original; the original's mode would not do, because the new file takes the process's owner and group. **On Windows** `Sys.chmod()` only sets the read-only attribute (`?Sys.chmod`), so there the backup has the permissions its folder gives any new file, the same as the stimulus file rcicr itself writes there. The backup is kept on Windows anyway, and `NEWS.md` says so.
   - Copy the original in with `file.copy(file, staging, overwrite = TRUE, copy.mode = FALSE)`, check that its `tools::md5sum()` equals the original's, then `file.rename(staging, backup)`. **Stop if the copy or checksum fails, or if the rename returns `FALSE`** (`file.rename()` reports a refusal, such as a scanner holding the file, by its return value, not an error).

   **Phase 2, save.** `save()` in place, as today. An error or interrupt here (Esc, Ctrl-C and R errors all run `on.exit()`): if the original's checksum differs from the backup's, copy the backup back in place with `file.copy(backup, file, overwrite = TRUE, copy.mode = FALSE)`; once the checksums match, remove the backup and re-raise the error. If the restore fails, keep the backup and give the restore command. When `save()` returns, the save is **committed**.

   **Phase 3, committed.** Remove the backup. **A committed save is never rolled back**: the new file and the backup differ by design, so treating a cleanup failure like a save failure would restore the old contents over the new ones. If the removal fails, warn that the save succeeded and name the backup; the next save's phase 0 finds that the original loads and removes it.

2. **A hard kill** (the R process killed, or crashing) skips `on.exit()`. Phase 0 of the next save handles every case:
   - in phase 1, only a staging file is left, beside an untouched original: it is reported;
   - in phase 2, the backup is left beside a possibly damaged original: if the original loads it is complete and the backup is removed; if not, the save stops with the restore command;
   - in phase 3, the original is complete, so the backup is removed.

   The restore **copies the backup's contents into the existing file**, never renames it over the original, which would replace the file and change its owner, group, ACLs and hard links: `file.copy("<file>.rcicr-backup", "<file>", overwrite = TRUE, copy.mode = FALSE)`, then delete the backup. `NEWS.md` and the error message both give this command, and say to use it only when the original no longer loads.
3. **Read-only targets behave as today.** `save()` fails to open a read-only file without truncating it (measured as an unprivileged user: the file's checksum is unchanged), the checksums then match, and the backup is removed.
4. **`NEWS.md`**, under Bug fixes: an interrupted or failed save of a reference distribution no longer destroys the stimulus file. The entry states every exception and change, so it claims no more than the code does:
   - **not protected, saved as before with a warning:** a writable file in a read-only directory (Unix), and a file name over 228 bytes;
   - **a new error:** on Windows, a writable file in a folder that refuses new files, which saves today; also, on any platform, when the backup cannot be made for another reason, such as a full disk;
   - on Windows the backup has its folder's permissions rather than being private to the user;
   - power loss and operating-system crashes are not covered;
   - the hard-kill restore command. Nothing under "Reproducibility impact": no number changes.

**What this does not cover: power loss or an operating-system crash.** A verified copy only shows the backup is readable through the operating system's cache, not that it has reached the disk, and base R cannot force that (it has no `fsync`). The guarantee is therefore limited to interrupts, R errors and a killed or crashed R process, and `NEWS.md` says so rather than promising more.

## Verified before planning

- **`save()` in place keeps the file's identity**: same inode, owner `nobody:nogroup` and mode 660 before and after (measured).
- **A restore with `file.copy(overwrite = TRUE, copy.mode = FALSE)` keeps it too**: same inode, owner and mode 660, and the restored file loads with its original value (measured).
- **`file.copy(copy.mode = TRUE)` and `Sys.chmod()` both apply the umask by default**: mode 660 became 640 (measured). Hence the pre-created backup, `use_umask = FALSE`, and `copy.mode = FALSE` throughout.
- **A failed `save()` on a read-only file leaves it untouched** (measured as an unprivileged user, since sessions here run as root).
- **A writable file in a read-only directory** can be saved in place, but no file can be created beside it (both measured as an unprivileged user). Hence the fallback to today's save.
- **`file.access(dir, 2)` identifies that case**: -1 for a mode-555 directory, 0 for a writable one (measured as an unprivileged user). On Windows, where `DECISIONS.md` calls `file.access()` unreliable, a false "not writable" on a full disk would trigger the unprotected save, so there any staging failure stops. The cost is an error on Windows for a writable file in a directory that refuses new files, a case that works today.
- **Content is unchanged**: the same objects are saved to the same path. A test asserts `identical()` loaded contents against the current behaviour.

## Tests

- Both call sites: after a successful save the loaded contents match, and no `.rcicr-backup` file remains.
- **Interruption:** a mocked writer that writes part of the file and then errors. The original loads with its original contents, has the same mode, and the backup is gone.
- A leftover backup beside an original that loads is removed with a message, and the save proceeds.
- A leftover backup beside an original that does not load stops the save with an error that gives the restore command, and neither file changes.
- An interrupt during the backup copy (phase 1) removes the staging file, leaves the original unchanged and attempts no restore (mocked).
- A failing backup copy stops before the original is touched, and leaves neither a staging file nor a backup (mocked).
- A `Sys.chmod()` that leaves the staging file at another mode stops before any data is copied (mocked; Unix only).
- A file name too long for the backup suffix falls back to today's save with a warning (mocked).
- Leftover staging files for the same file are reported by name in a warning and left in place (mocked).
- A refused publishing rename stops before the original is touched, and leaves no staging file (mocked).
- A failure to delete the backup after a successful save warns, keeps the new contents, and does not restore the backup (mocked).
- A staging file that cannot be created falls back to today's in-place save when the directory is not writable, and stops without touching the original when it is (both mocked).
- The backup is mode 600 whatever the original's mode (skipped on Windows).
- File mode 600 and 660 are preserved (skipped on Windows).

## The step most likely to fail

**The restore on Windows.** Copying the backup back in place could fail there when a virus scanner or indexer holds the file open. The helper then keeps the backup and names it in the error, so the data survives, but the user has to restore it by hand. The `windows-latest (release)` job runs the interruption test.

## Out of scope

- #334 (a file without `seed`): a separate change.
- Loading the stimulus file two or three times per InfoVal call: a speed issue, not correctness.
