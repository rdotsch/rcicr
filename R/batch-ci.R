# The loop behind batchGenerateCI() and batchGenerateCI2IFC(), which differ
# only in argument order (#87). Called by name from both: the two put `label`
# in different positions, so a positional call would bind it wrongly in one.
batchCIs <- function(data, by, stimuli, responses, baseimage, rdata, save_as_png, targetpath,
                     label, antiCI, scaling, constant, participants) {

  # targetpath is required, not defaulted: a default path writes to the user's
  # filespace uninvited, which CRAN policy does not allow.
  if (save_as_png && missing(targetpath)) {
    stop(paste0('save_as_png is TRUE but no targetpath was given. Supply ',
      'targetpath = <a directory> to say where the PNGs should go, ',
      'or set save_as_png = FALSE to compute the classification ',
      'images without writing them. Use tempdir() if you only want ',
      'to try the function out.'
    ))
  }

  if (scaling == 'autoscale') {
    doAutoscale <- TRUE
    scaling <- 'none'
  } else {
    doAutoscale <- FALSE
  }

  # Rows without a grouping value cannot name an output CI and should not be
  # turned into a spurious NA group. Columns are read with [[ ]], never
  # [, col]: on a tibble the latter stays a one-column tibble, which the loop
  # and the == below treat as one unit, mixing every group into a single CI
  # (#336).
  data <- data[!is.na(data[[by]]), , drop = FALSE]
  by.levels <- unique(data[[by]]) # nolint: object_name_linter.
  if (!is.null(participants)) requireParticipantIds(data, by, participants)

  # dplyr::progress_estimated() is deprecated; use the base R progress bar
  pb <- txtProgressBar(min = 0, max = length(by.levels), style = 3)
  cis <- list()
  pb_i <- 0

  for (unit in by.levels) {
    pb_i <- pb_i + 1
    setTxtProgressBar(pb, pb_i)

    unitdata <- data[data[[by]] == unit, , drop = FALSE]

    if (label == '') {
      filename <- paste0(baseimage, '_', by, '_', unitdata[[by]][1])
    } else {
      filename <- paste0(baseimage, '_', label, '_', by, '_', unitdata[[by]][1])
    }

    cis[[filename]] <- generateCI(
      stimuli = unitdata[[stimuli]],
      responses = unitdata[[responses]],
      baseimage = baseimage,
      rdata = rdata,
      save_as_png = save_as_png,
      filename = filename,
      targetpath = targetpath,
      antiCI = antiCI,
      scaling = scaling,
      scaling_constant = constant,
      participants = if (is.null(participants)) NA else unitdata[[participants]]
    )
  }

  if (doAutoscale) {
    cis <- autoscale(cis, save_as_pngs = save_as_png, targetpath = targetpath)
  }

  close(pb)
  return(cis)
}

# Checked for the whole table before any CI: generateCI() reads an all-NA
# participants vector as "no participants" and pools the trials, so a group
# whose IDs were all missing would be pooled while the others were averaged
# per participant.
requireParticipantIds <- function(data, by, participants) {
  if (!is.character(participants) || length(participants) != 1L || is.na(participants) ||
        !participants %in% names(data)) {
    stop('participants must name one column of data; ',
         encodeString(as.character(participants)[1], quote = "'"), ' is not one.', call. = FALSE)
  }
  missing_ids <- is.na(data[[participants]])
  if (any(missing_ids)) {
    groups <- unique(data[[by]][missing_ids])
    shown <- paste(utils::head(groups, 5), collapse = ', ')
    if (length(groups) > 5) shown <- paste0(shown, ', ...')
    stop(sum(missing_ids), ' rows have no participant ID in column ', participants, ' (', by,
         ' ', shown, '). Give every trial an ID, or remove the trials without one.',
         call. = FALSE)
  }
  invisible(NULL)
}
