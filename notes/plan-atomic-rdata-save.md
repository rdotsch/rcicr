# Plan: save reference distributions without risking the stimulus file (#333)

## Problem

Two `save()` calls overwrite the stimulus `.Rdata` file in place: `R/generateReferenceDistribution.R:174` (shared parameters) and `R/reference-base.R:175` (independent base images). `save()` truncates the file before writing, so an interrupted write destroys it. #333 has the reproduction: killed mid-write, a valid 733-byte file became 10,993,664 bytes that `load()` rejects.

## Change

1. **One internal helper, `saveRdataAtomically(names, file, envir)`, in `R/rdata.R`**, used by both call sites:
   - resolve the target with `normalizePath(file, mustWork = TRUE)`;
   - `save()` to `tempfile(pattern = ".rcicr-", tmpdir = dirname(target), fileext = ".Rdata")`;
   - copy the original's mode onto the temporary file with `Sys.chmod(tmp, file.mode(target))`;
   - `file.rename(tmp, target)`; if that returns `FALSE`, stop with an error saying the original is unchanged;
   - `on.exit(unlink(tmp))`, so a failed save leaves no temporary file behind.
2. **Readonly detection also checks the directory.** A rename needs a writable directory, not just a writable file. `resolveReferenceNorms()` treats a file in an unwritable directory as read-only: it scores from memory and says so, as it already does for an unwritable file. Previously such a file was saved to; this is the one behaviour change.
3. **`NEWS.md`**, under Bug fixes: an interrupted or failed save of a reference distribution no longer destroys the stimulus file. Nothing under "Reproducibility impact": no number changes.

## Verified before planning

- **Rename across filesystems fails** (`?file.rename`), so the temporary file goes in the target's own directory, not `tempdir()`.
- **A rename through a symlink replaces the link** with a regular file and leaves the real target unchanged (measured here). Resolving with `normalizePath()` first keeps the link and updates the target (measured).
- **A rename resets the permissions**: a 640 file came back 644 (measured). Hence the `Sys.chmod()`.
- **Content is unchanged**: the same objects go to `save()`, only to a different path first. A test asserts that `load()` of the result is `identical()` to a direct `save()`.

## Tests

- The saved file's contents equal those of a direct `save()`, for both call sites.
- **Interruption:** mock the internal writer to write part of a file and then error. The original still loads with its original contents, and no `.rcicr-*` file is left in the directory.
- Symlink preserved, mode preserved (both skipped on Windows, where creating symlinks needs Developer Mode).
- A file in an unwritable directory is scored from memory (mocked, like the existing `writableFile()` tests, because sessions here run as root).

## The step most likely to fail

**`file.rename()` over an existing file on Windows.** `?file.rename` says it overwrites "where file permissions allow", but that has only been verified on Linux here. A virus scanner or indexer holding the target open can also make the rename fail. The `windows-latest (release)` job runs the new tests. If the rename fails, the error leaves the original intact; that is the property the fix is for.

## Out of scope

- #334 (a file without `seed`): a separate change.
- Loading the stimulus file two or three times per InfoVal call: a speed issue, not correctness.
