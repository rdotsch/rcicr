# CRAN Note avoidance

if (getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      # Suppress checking notes for variables loaded at runtime from .RData files
      "p", "s", "base_faces", "stimuli_params", "img_size", "base_face_files", "n_trials", "seed", "noise_type", "reference_norms", "reference_norms_seed",

      # Latent-space module: the fields generateStimuliLatent2IFC() saves
      "latent_params", "base_latent", "latent_sigma", "generator_spec",

      # Suppress checking notes for variables in foreach loop (parallel runs)
      "obs"
    )
  )
