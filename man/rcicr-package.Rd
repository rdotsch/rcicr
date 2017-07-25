\name{rcicr-package}
\alias{rcicr-package}
\alias{rcicr}
\docType{package}
\title{
Reverse-Correlation Image-Classification Toolbox
}
\description{
Toolbox with functions to generate stimuli and analyze data of reverse correlation image classification experiments. Reverse correlation is a psychophysical technique originally derived from signal detection theory. This package focuses on visualizing internal representations of participants using visual stimuli in a perceptual taks.
}
\details{
\tabular{ll}{
Package: \tab rcicr\cr
Type: \tab Package\cr
Version: \tab 0.4.0\cr
Date: \tab 2017-07-25\cr
License: \tab GPL-2\cr
}
\bold{Generating stimuli}

Load the package with \code{library(rcicr)}. Then generate stimuli with: 

\code{generateStimuli2IFC(base_face_files, n_trials = 770)}

This will generate stimuli for 770 trials of a 2 images forced choice reverse correlation image classification task with sinusoid noise. By default the stimuli will have a resolution of 512 x 512 pixels. The stimuli will be saved as jpegs to a folder called stimuli in your current working directory, and an .Rdata file will be saved that contains the stimulus parameters necessary for analysis.

The \code{base_face_files} argument is a list of jpegs that should be used as base images for the stimuli. The base_face_files variable might look like this:

\code{base_face_files <- list('male'='male.jpg', 'female'='female.jpg')}

For each jpeg a set of stimuli will be created using the same noise patterns as for the other sets. The jpeg should have the resolution that you want the stimuli to have. By default this should be 512 x 512 pixels. If you want a different size, resize your base image to either 128 x 128 or 256 x 256 for smaller stimuli, or 1024 x 1024 for bigger stimuli. In that case, also set the \code{img_size} parameter accordingly.

You are now ready to collect data with the stimuli you just created. The stimuli are named according to their sequence number when generating and whether the original noise is superimposed or the negative/inverted noise. Stimuli with the same sequence number should be presented side by side in the same trial. Record which stimulus a participant selected at any given trial (the original, or the inverted). At the very least be sure that in your data file the connection can be made between the response key of the participant and which stimulus was selected on each trial. Use any presentation software you like (I recommend python-based open source alternatives, like PsychoPy, Expyriment, or OpenSesame).

\bold{Data analysis}

Analyzing reverse correlation data is all about computing classification images. Use the following function for your data collected using the 2 images forced choice stimuli:

\code{ci <- generateCI2IFC(stimuli, responses, baseimage, rdata)}

The \code{stimuli} paramater should be a vector containing the sequence numbers of the stimuli that were presented in the trials of the task. The \code{responses} parameter contains, in the order of the stimuli vector, the response of a participant to those stimuli (coded 1 if the original stimulus was selected and -1 if the inverted stimulus was selected). The \code{baseimage} parameter is a string specifying which base image was used (not the file name, but the name in the list of \code{base_face_files}. So for the stimuli generated above, either \code{'male'} or \code{'female'}, depending on which set of stimuli was presented to the participant whose data you're analyzing). Finally, rdata is a string pointing to the .RData file that was created automatically when you generated the stimuli. It contains the parameters for each stimulus, necessary to create the classification image.

By default JPEGs of the classification images will be saved automatically. The returned values can be used later to optimally rescale the noise relative to the base image. For instance, if you have a list of cis from various participants (i.e., a list of the values returned by several calls to \code{generateCI2IFC}, one for each participant), you can use the autoscale function to generate classification images that are scaled identically and therefore straightforward to compare:

\code{scaled_cis <- autoscale(cis, saveasjpegs = TRUE)}

\bold{Computing CIs for many participants or conditions}

Data analysis as described above can be automatized for a batch of participants or conditions using \code{batchGenerateCI2IFC}. Please see instructions for that function. 

\bold{Note}

Currently, the package is still in alpha stage. Much may still change. It only supports 2 Image Forced Choice tasks, although the underlying functions can be used for other versions of the reverse correlation task. It also only supports sinusoid noise. In the future, it will support Gaussian white noise, as well as additional variants of the task.

If you use this package for your experiments, please cite the package in your publications. Use \code{citation('rcicr')} to print the appropriate citation for the current version of the package.
}
\author{
Ron Dotsch <rdotsch@gmail.com> (http://ron.dotsch.org/)
Maintainer: Ron Dotsch <rdotsch@gmail.com>
}
\references{
Dotsch, R., & Todorov, A. (2012). Reverse correlating social face perception. Social Psychological and Personality Science, 3 (5), 562-571.

Dotsch, R., Wigboldus, D. H. J., Langner, O., & Van Knippenberg, A. (2008). Ethnic out-group faces are biased in the prejudiced mind. Psychological Science, 19, 978-980.
}
\keyword{ package }

\examples{
#simple examples will be added soon.
}
