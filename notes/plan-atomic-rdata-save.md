# Plan: save reference distributions without risking the stimulus file (#333)

## Problem

Two `save()` calls overwrite the stimulus `.Rdata` file in place: `R/generateReferenceDistribution.R:174` (shared parameters) and `R/reference-base.R:175` (independent base images). `save()` truncates the file before writing, so an interrupted write destroys it. #333 has the reproduction: killed mid-write, a valid 733-byte file became 10,993,664 bytes that `load()` rejects.

## Approach: back up, then save in place

The file keeps being written in place, exactly as today, so its owner, group, ACLs, mode, hard links and symlinks are untouched. What changes is that a complete copy exists for the duration of the write.

Writing a temporary file and renaming it over the original was rejected: a rename replaces the file, so it takes the process's owner and group (measured: `nobody:nogroup` became `root:root`), drops ACLs, and base R can restore neither. For a shared archive that silently changes who can read the stimulus file.

## Change

1. **One internal helper, `saveRdataSafely(names, file, envir)`, in `R/rdata.R`**, used by both call sites:
   - `backup <- paste0(file, ".rcicr-backup")`. **If `backup` already exists, stop** before touching anything, naming both files. A leftover backup means an earlier save was killed and the original may be damaged; overwriting the backup with it could destroy the only good copy;
   - create `backup` empty, give it the original's mode with `Sys.chmod(backup, file.mode(file), use_umask = FALSE)`, then copy the original into it with `file.copy(file, backup, overwrite = TRUE, copy.mode = FALSE)`. If the copy fails (for example, a full disk), remove the backup and stop; the original has not been touched;
   - `save()` in place, as today;
   - **on success**, remove the backup;
   - **on error or interrupt** (Esc, Ctrl-C and R errors all run `on.exit()`): if the original's checksum (`tools::md5sum()`) differs from the backup's, copy the backup back in place with `file.copy(backup, file, overwrite = TRUE, copy.mode = FALSE)`. Remove the backup once the checksums match, then re-raise the error. If the restore fails, keep the backup and say where it is.
2. **A hard kill** (the process killed, a crash, power loss) skips `on.exit()`. It then leaves the complete backup beside a possibly damaged original, and the leftover-backup check stops the next save. `NEWS.md` says how to restore: rename `<file>.rcicr-backup` back to `<file>`.
3. **Read-only targets behave as today.** `save()` fails to open a read-only file without truncating it (measured as an unprivileged user: the file's checksum is unchanged), the checksums then match, and the backup is removed. No separate writability check is needed.
4. **`NEWS.md`**, under Bug fixes: an interrupted or failed save of a reference distribution no longer destroys the stimulus file, plus the hard-kill note above. Nothing under "Reproducibility impact": no number changes.

## Verified before planning

- **`save()` in place keeps the file's identity**: same inode, owner `nobody:nogroup` and mode 660 before and after (measured).
- **A restore with `file.copy(overwrite = TRUE, copy.mode = FALSE)` keeps it too**: same inode, owner and mode 660, and the restored file loads with its original value (measured).
- **`file.copy(copy.mode = TRUE)` and `Sys.chmod()` both apply the umask by default**: mode 660 became 640 (measured). Hence the pre-created backup, `use_umask = FALSE`, and `copy.mode = FALSE` throughout.
- **A failed `save()` on a read-only file leaves it untouched** (measured as an unprivileged user, since sessions here run as root).
- **Content is unchanged**: the same objects are saved to the same path. A test asserts `identical()` loaded contents against the current behaviour.

## Tests

- Both call sites: after a successful save the loaded contents match, and no `.rcicr-backup` file remains.
- **Interruption:** a mocked writer that writes part of the file and then errors. The original loads with its original contents, has the same mode, and the backup is gone.
- A leftover backup stops the save with an error, and neither file changes.
- A failing backup copy stops before the original is touched (mocked).
- File mode 600 and 660 are preserved (skipped on Windows).

## The step most likely to fail

**The restore on Windows.** Copying the backup back in place could fail there when a virus scanner or indexer holds the file open. The helper then keeps the backup and names it in the error, so the data survives, but the user has to restore it by hand. The `windows-latest (release)` job runs the interruption test.

## Out of scope

- #334 (a file without `seed`): a separate change.
- Loading the stimulus file two or three times per InfoVal call: a speed issue, not correctness.
