# Plan: save reference distributions without risking the stimulus file (#333)

## Problem

Two `save()` calls overwrite the stimulus `.Rdata` file in place: `R/generateReferenceDistribution.R:174` (shared parameters) and `R/reference-base.R:175` (independent base images). `save()` truncates the file before writing, so an interrupted write destroys it. #333 has the reproduction: killed mid-write, a valid 733-byte file became 10,993,664 bytes that `load()` rejects.

## Approach: back up, then save in place

The file keeps being written in place, exactly as today, so its owner, group, ACLs, mode, hard links and symlinks are untouched. What changes is that a complete copy exists for the duration of the write.

Writing a temporary file and renaming it over the original was rejected: a rename replaces the file, so it takes the process's owner and group (measured: `nobody:nogroup` became `root:root`), drops ACLs, and base R can restore neither. For a shared archive that silently changes who can read the stimulus file.

## Change

1. **One internal helper, `saveRdataSafely(names, file, envir)`, in `R/rdata.R`**, used by both call sites:
   - `backup <- paste0(file, ".rcicr-backup")`. **If `backup` already exists, stop** before touching anything, naming both files and the restore command below. A backup only exists at that path once it is complete (next step), so a leftover one means an earlier save was killed mid-write and the original may be damaged. Copying that original over the backup could destroy the only good copy;
   - **make the backup under a staging name, then publish it.** Create `staging <- tempfile(pattern = paste0(basename(file), ".rcicr-staging-"), tmpdir = dirname(file))` empty, set it to mode 600 with `Sys.chmod(staging, "600", use_umask = FALSE)`, **on Unix-alikes check that `file.mode(staging)` is then 600 and stop if not** (`Sys.chmod()` reports failure by its return value, and a staging file left at its creation mode would expose the data after a hard kill), and copy the original into it with `file.copy(file, staging, overwrite = TRUE, copy.mode = FALSE)`. Check that its `tools::md5sum()` equals the original's, then `file.rename(staging, backup)`, and **stop if that returns `FALSE`**: `file.rename()` reports a refusal (a scanner holding the file, a path-length limit) by its return value, not an error, and continuing would save without a backup. Mode 600 lets only the user running R read the backup, who could already read the original. The original's mode would not do: the backup takes the process's owner and group, so a group-readable original would give a group-readable backup for a different group. **On Windows** `Sys.chmod()` only sets the read-only attribute (`?Sys.chmod`), so there the backup has the permissions its folder gives any new file, the same as the stimulus file rcicr itself writes there. That is broader than the original only when the original was given tighter permissions than its folder. The backup is kept on Windows anyway: dropping it would leave Windows users without the protection this change exists for, and `NEWS.md` states the Windows behaviour;
   - **if the staging file's name would exceed 255 bytes** (`nchar(basename(staging), type = "bytes")`, the usual file-name limit; the staging name is longer than the backup name, so it decides), save in place without a backup and warn that the name is too long to protect. The original's name must then be over 228 bytes (the staging suffix is 27, measured); rcicr's own names are about 40. That limit can be recognised before anything is written, so this fallback does not mask a full disk;
   - **if the staging file cannot be created**, check the directory with `file.access(dirname(file), 2)`. Only on Unix-alikes, and only when that reports it not writable (a read-only directory holding a writable file) does the helper save in place without a backup, exactly as today: that case works today, and protecting it would need somewhere to write. **Any other failure stops before the original is touched**: a directory that reports writable but refuses the file is out of space, inodes or quota, and saving in place there is the very failure this change prevents;
   - if the copy or the checksum fails (for example, a full disk), remove the staging file and stop; the original has not been touched;
   - `save()` in place, as today, then mark the save as **committed**;
   - **after a committed save**, remove the backup. If that fails, warn that the save succeeded and name the backup to delete. **Never roll back a committed save**: the new file and the backup differ by design, so treating a cleanup failure like a save failure would restore the old contents over the new ones. The leftover backup then stops the next save until it is deleted, which the warning says;
   - **on error or interrupt before the save committed** (Esc, Ctrl-C and R errors all run `on.exit()`): if the original's checksum differs from the backup's, copy the backup back in place with `file.copy(backup, file, overwrite = TRUE, copy.mode = FALSE)`. Remove the backup once the checksums match, then re-raise the error. If the restore fails, keep the backup and give the restore command.
2. **A hard kill** (the R process killed, or crashing) skips `on.exit()`:
   - during the backup copy, it leaves only a staging file; the original is untouched and the next save proceeds. **Each save lists leftover `<file>.rcicr-staging-*` files for that file and warns with their names**, saying they are incomplete copies that can be deleted, so repeated interruptions cannot quietly fill the disk. The helper does not delete them itself: a file this call did not create could belong to another R session saving at the same moment;
   - during the save, it leaves the complete backup beside a possibly damaged original, and the leftover-backup check stops the next save.

   The restore **copies the backup's contents into the existing file**, never renames it over the original, which would replace the file and change its owner, group, ACLs and hard links: `file.copy("<file>.rcicr-backup", "<file>", overwrite = TRUE, copy.mode = FALSE)`, then delete the backup. `NEWS.md` and the error message both give this command.
3. **Read-only targets behave as today.** `save()` fails to open a read-only file without truncating it (measured as an unprivileged user: the file's checksum is unchanged), the checksums then match, and the backup is removed.
4. **`NEWS.md`**, under Bug fixes: an interrupted or failed save of a reference distribution no longer destroys the stimulus file, when the file's directory is writable. It includes the hard-kill restore command, and says that on Windows the backup has its folder's permissions rather than being private to the user. Nothing under "Reproducibility impact": no number changes.

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
- A leftover backup stops the save with an error that gives the restore command, and neither file changes.
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
