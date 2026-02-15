
# install.packages('devtools')
# devtools::install_github('rdotsch/rcicr')

# Load reverse correlation toolbox
library(rcicr)

# Set base face
base = list('W_M_O'='W_M_O.jpg')

# Generate and save stimuli
generateStimuli2IFC(base_face_files = base, n_trials=20, stimulus_path = "./stimuli", label='preconf', nscales=5, noise_type='gabor', sigma=25)