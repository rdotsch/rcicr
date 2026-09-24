# generateStimuli2IFC() never overwrites a stimulus .Rdata file (#338). Before
# anything is generated it takes an atomic dir.create() lock beside the file,
# keyed on seed and minute, never on label; DECISIONS.md, "A stimulus .Rdata
# file is never overwritten", says why.

# Named so tests can fix the clock.
stimulusTime <- function() Sys.time()

stimulusRdataPath <- function(stimulus_path, label, seed, time) {
  paste(stimulus_path, paste(label, "seed", seed, "time",
                             format(time, format = "%b_%d_%Y_%H_%M.Rdata"), sep = "_"), sep = "/")
}

stimulusLockPath <- function(rdata_file, seed, time) {
  key <- paste0("seed_", paste(seed, collapse = ","), "_", format(time, format = "%b_%d_%Y_%H_%M"))
  # R 4.1 has no string hash.
  holder <- tempfile()
  on.exit(unlink(holder), add = TRUE)
  writeLines(key, holder, useBytes = TRUE)
  # Not dirname(), which fails on a label the native encoding cannot represent.
  target_dir <- sub("[/\\\\][^/\\\\]*$", "", rdata_file)
  file.path(target_dir, paste0(".rcicr-lock-", unname(tools::md5sum(holder))))
}

# Returns the lock, which the caller removes once the file is written or the call fails.
acquireStimulusLock <- function(rdata_file, seed, time) {
  lock <- stimulusLockPath(rdata_file, seed, time)
  if (!dir.create(lock, showWarnings = FALSE)) {
    if (!dir.exists(lock)) {
      stop("Could not create ", lock, " to reserve ", rdata_file,
           "; check that the folder exists and is writable.", call. = FALSE)
    }
    stop(rdata_file, " is reserved by ", lock, ". It belongs to a generateStimuli2IFC() call ",
         "into this folder with the same seed, started in the same minute, which is either ",
         "still running or was interrupted. Delete it only after confirming that no ",
         "generateStimuli2IFC() call is still writing to this folder.", call. = FALSE)
  }
  if (file.exists(rdata_file)) {
    unlink(lock, recursive = TRUE)
    stop(rdata_file, " already exists, and a stimulus file is never overwritten. Use a ",
         "different label or stimulus_path, or delete the file if you are regenerating this ",
         "stimulus set on purpose.", call. = FALSE)
  }
  lock
}
