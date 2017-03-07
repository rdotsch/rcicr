# CRAN Note avoidance


if (getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      # Suppress checking notes for variables loaded at runtime from .RData files
      "p", "s", "base_faces", "stimuli_params", "img_size",

      # Suppress checking notes for variables in foreach loop (parallel runs)
      "obs"
      )
  )