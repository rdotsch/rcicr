
# Suppress checking notes for variables loaded at runtime from .RData files
if(getRversion() >= "2.15.1")  utils::globalVariables(c("p", "s", "base_faces", "stimuli_params", "img_size"))
