# generateStimuli2IFC() never overwrites a stimulus .Rdata file (#338). Before
# anything is generated it takes an atomic dir.create() lock beside the file,
# keyed on seed and minute, never on label; DECISIONS.md, "Stimulus files are
# never overwritten", says why.

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
  file.path(targetDir(rdata_file), paste0(".rcicr-lock-", unname(tools::md5sum(holder))))
}

# dirname() with the host's separators and roots, without its translation to the
# native encoding, which fails on a label that encoding cannot represent.
targetDir <- function(path) {
  windows <- .Platform$OS.type == "windows"
  seps <- if (windows) "/\\\\" else "/"
  if (!grepl(sprintf("[%s]", seps), path)) return(".")
  dir <- sub(sprintf("[%s][^%s]*$", seps, seps), "", path)
  if (!nzchar(dir) || (windows && grepl("^[A-Za-z]:$", dir))) paste0(dir, "/") else dir
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
  if (pathTaken(rdata_file)) {
    unlink(lock, recursive = TRUE)
    stop(rdata_file, " already exists, and a stimulus file is never overwritten. Use a ",
         "different label or stimulus_path, or delete the file if you are regenerating this ",
         "stimulus set on purpose.", call. = FALSE)
  }
  lock
}

# The PNGs are never overwritten either (#350). Their names carry no time, so a
# later-minute call passes the .Rdata reservation above.
stimulusPngPath <- function(stimulus_path, label, base_face, seed, trial, side) {
  paste(stimulus_path, paste(label, base_face, seed, sprintf("%05d_%s.png", trial, side), sep = "_"),
        sep = "/")
}

stimulusPngPaths <- function(stimulus_path, label, base_labels, seed, n_trials) {
  unlist(lapply(base_labels, function(base_face) {
    lapply(seq_len(n_trials), function(trial) {
      vapply(c("ori", "inv"), function(side) {
        stimulusPngPath(stimulus_path, label, base_face, seed, trial, side)
      }, character(1), USE.NAMES = FALSE)
    })
  }), use.names = FALSE)
}

# Named so tests can model a case-insensitive file system. A dangling symlink
# is taken too: file.exists() follows it and says no, and file.create() would
# create its target, which may lie outside the folder. Sys.readlink() gives NA
# for a path that does not exist, and "" for one that is not a link.
pathTaken <- function(path) {
  link <- Sys.readlink(path)
  file.exists(path) || (!is.na(link) && nzchar(link))
}

# Keyed on the seed alone. Same-seed calls into one folder write the same PNG
# names whatever the minute, and the label cannot be part of the key, for the
# reason #338's lock gives.
acquirePngLock <- function(stimulus_path, seed) {
  holder <- tempfile()
  on.exit(unlink(holder), add = TRUE)
  writeLines(paste0("png_seed_", paste(seed, collapse = ",")), holder, useBytes = TRUE)
  lock <- file.path(stimulus_path, paste0(".rcicr-png-lock-", unname(tools::md5sum(holder))))
  if (!dir.create(lock, showWarnings = FALSE)) {
    if (!dir.exists(lock)) {
      stop("Could not create ", lock, " to reserve the stimulus PNGs; check that ", stimulus_path,
           " is writable.", call. = FALSE)
    }
    stop("The stimulus PNGs in ", stimulus_path, " are reserved by ", lock, ". It belongs to a ",
         "generateStimuli2IFC() call into this folder with the same seed, which is either still ",
         "running or was interrupted. Delete it only after confirming that no generateStimuli2IFC() ",
         "call is still writing to this folder.", call. = FALSE)
  }
  lock
}

# Each path is looked up, then created empty, in order, so the file system
# decides what is the same name: one taken by an earlier run, or by this call
# a moment ago under a spelling the file system treats as equal. Returns the
# reserved paths; the caller removes them if it does not finish.
reserveStimulusPngs <- function(paths) {
  shown <- function(x) paste0(paste(utils::head(x, 3), collapse = ", "), if (length(x) > 3) ", ...")
  existing <- paths[vapply(paths, pathTaken, logical(1))]
  if (length(existing) > 0) {
    # Only a call that never ran its cleanup (R was killed) leaves these.
    leftover <- isTRUE(all(file.size(existing) == 0))
    stop(length(existing), " of the ", length(paths), " PNG files this call would write already ",
         "exist (", shown(existing), "), and stimulus files are never overwritten. ",
         if (leftover) paste0("They are all empty: placeholders left by a generateStimuli2IFC() ",
                              "call that was killed before it could clean up, safe to delete. "),
         "Use a different label or stimulus_path, or, to regenerate this stimulus set on purpose, ",
         "delete its PNGs and its .Rdata file.", call. = FALSE)
  }
  reserved <- character()
  # Rolled back if this loop is interrupted, since the caller only learns what
  # was reserved once it returns.
  done <- FALSE
  on.exit(if (!done) unlink(reserved), add = TRUE)
  for (path in paths) {
    if (pathTaken(path) || !file.create(path, showWarnings = FALSE)) {
      aliased <- pathTaken(path)
      unlink(reserved)
      if (aliased) {
        stop(path, " is the same file, on this file system, as another PNG this call writes: two ",
             "base image labels differ only in a way the file system ignores, such as case. Give ",
             "the base images labels that differ in more than that.", call. = FALSE)
      }
      stop("Could not create ", path, "; check that the folder is writable.", call. = FALSE)
    }
    reserved <- c(reserved, path)
  }
  done <- TRUE
  reserved
}

# The exit handler of generateStimuli2IFC(). Workers are killed rather than
# stopped: a stopped worker still finishes the trial it is writing.
releaseStimulusCall <- function(finished, cl, worker_pids, reserved_pngs, owned_rdata, locks) {
  if (!finished) {
    # pskill() takes the whole vector; try() keeps a failure from skipping the rest.
    if (length(worker_pids) > 0) try(tools::pskill(worker_pids, tools::SIGKILL), silent = TRUE)
    stopClusterSafely(cl)
    unlink(c(reserved_pngs, owned_rdata))
  }
  unlink(locks, recursive = TRUE)
  invisible(NULL)
}

# Named so tests can make the save fail part-way.
saveStimulusFile <- function(names, file, envir) {
  save(list = names, file = file, envir = envir)
}
