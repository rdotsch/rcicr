# Reading a base image is one operation shared by both pipelines: the
# pixel-noise stimuli need one base face per label, the PCA generator in
# latentGeneratorPCA() needs a whole set, and both reject the same malformed
# inputs for the same reasons. The advice closing each message names the
# caller, so a reader is told which of the two calls to change.

# Which reader a base image needs: 'png', 'jpeg', or NA.
#
# Anchored to the extension. The old grepl('png|PNG', filename) matched anywhere
# in the path, so a JPEG under a directory named "png" went to png::readPNG().
baseImageFormat <- function(filename) {
  if (grepl('\\.png$', filename, ignore.case = TRUE)) {
    return('png')
  }
  if (grepl('\\.jpe?g$', filename, ignore.case = TRUE)) {
    return('jpeg')
  }
  NA_character_
}

readBaseImage <- function(filename, label, img_size, maximize_contrast, caller) {
  img_format <- baseImageFormat(filename)

  img <- tryCatch(
    if (img_format == 'png') png::readPNG(filename) else jpeg::readJPEG(filename),
    error = function(e) {
      stop(paste0('Base image "', label, '" (', filename, ') could not ',
          'be read as ', img_format, ': ', conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  # Check if base face is square. If not, throw an error
  if (dim(img)[1] != dim(img)[2]) {
    stop(paste0('Base image "', label, '" (', filename, ') is not ',
      'square! It\'s ', dim(img)[1], ' by ', dim(img)[2],
      ' pixels. Please use a square base face.'
    ))
  }

  # Change base face to greyscale if necessary.
  #
  # This averages alpha in along with the colour channels, which is wrong: #295
  # fixes it in PR #296, against main, where the bug has been since v1.0.1.
  # Left alone here so this branch stays a pure extraction and the two changes
  # do not arrive in one diff.
  if (length(dim(img)) == 3) {
    img <- apply(img, c(1, 2), mean)
  }

  # Check that the base face matches the requested stimulus size. Automatic
  # resizing used to happen here via biOps, but that dependency was dropped
  # and never replaced, so a mismatch would otherwise surface much later as
  # an opaque "non-conformable arrays" error from inside a parallel worker
  # (when the noise is added to the base image).
  # img_size = NULL skips the check, for a caller that is reading the first
  # image in order to learn what size the set is.
  if (!is.null(img_size) && nrow(img) != img_size) {
    stop(paste0('Base image "', label, '" (', filename, ') is ',
      nrow(img), ' by ', ncol(img), ' pixels, but img_size is ',
      img_size, '. rcicr does not resize base images: please ',
      'either resize the image to ', img_size, ' by ', img_size,
      ' pixels, or call ', caller, '() with img_size = ',
      nrow(img), '.'
    ))
  }

  # If necessary, rescale to maximize contrast
  if (maximize_contrast) {
    # (img - min) / (max - min) is 0/0 on a uniform image, and the all-NaN
    # base face went into the .Rdata unremarked (#176). Only an error here:
    # with the rescale off a flat base image is usable and produces valid
    # stimuli, so rejecting it outright would break a legitimate call.
    if (max(img) == min(img)) {
      stop(paste0('Base image "', label, '" (', filename, ') has no ',
        'contrast: every pixel is ', min(img), '. Contrast cannot ',
        'be maximized on a uniform image, and doing so would make ',
        'the base image entirely NaN. Use a base image with some ',
        'variation, or call ', caller, '() with ',
        'maximize_baseimage_contrast = FALSE.'
      ))
    }

    img <- (img - min(img)) / (max(img) - min(img))
  }

  return(img)
}

# Indexed by position rather than by name. Looking each image up by its label
# returns the first entry carrying that label however many share it, and then
# writes every one of them into the same slot, so a duplicated label silently
# drops images from the set. generateStimuli2IFC() is protected from that by
# validateBaseFaceFiles(), which requires unique names; nothing protected the
# latent pipeline, which builds its labels from file names.
readBaseImages <- function(base_face_files, img_size, maximize_contrast, caller) {
  labels <- names(base_face_files)
  base_faces <- vector('list', length(base_face_files))
  names(base_faces) <- labels

  for (i in seq_along(base_face_files)) {
    base_faces[[i]] <- readBaseImage(
      base_face_files[[i]], labels[i], img_size, maximize_contrast, caller
    )
  }

  return(base_faces)
}
