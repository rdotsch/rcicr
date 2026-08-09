# Package index

## Generating stimuli

Build a stimulus set for a 2IFC task. The .Rdata file written here is
the only link to the analysis half, so keep it.

- [`generateStimuli2IFC()`](https://rdotsch.github.io/rcicr/reference/generateStimuli2IFC.md)
  : Generates 2IFC stimuli
- [`simulateNoiseIntensities()`](https://rdotsch.github.io/rcicr/reference/simulateNoiseIntensities.md)
  : Simulate pixel intensity range for noise

## Classification images

Turn collected responses into classification images, and render them.

- [`generateCI()`](https://rdotsch.github.io/rcicr/reference/generateCI.md)
  : Generates classification image
- [`generateCI2IFC()`](https://rdotsch.github.io/rcicr/reference/generateCI2IFC.md)
  : Generates 2IFC classification image
- [`batchGenerateCI()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI.md)
  : Generates multiple classification images by participant or condition
- [`batchGenerateCI2IFC()`](https://rdotsch.github.io/rcicr/reference/batchGenerateCI2IFC.md)
  : Generates multiple 2IFC classification images by participant or
  condition
- [`autoscale()`](https://rdotsch.github.io/rcicr/reference/autoscale.md)
  : Determines optimal scaling constant for a list of CIs
- [`plotZmap()`](https://rdotsch.github.io/rcicr/reference/plotZmap.md)
  : Plots a Z-map
- [`computeCumulativeCICorrelation()`](https://rdotsch.github.io/rcicr/reference/computeCumulativeCICorrelation.md)
  : Computes cumulative trial CIs correlations with final/target CI

## Informational value

How much signal a classification image carries, against a simulated
null.

- [`computeInfoVal2IFC()`](https://rdotsch.github.io/rcicr/reference/computeInfoVal2IFC.md)
  : Computes Informational Value
- [`generateReferenceDistribution2IFC()`](https://rdotsch.github.io/rcicr/reference/generateReferenceDistribution2IFC.md)
  : Generates reference distribution

## Noise basis

The sinusoid and Gabor machinery the stimuli are built from. Exported so
the noise can be inspected and reused; most analyses never call these
directly.

- [`generateNoisePattern()`](https://rdotsch.github.io/rcicr/reference/generateNoisePattern.md)
  : Generate sinusoid noise pattern
- [`generateNoiseImage()`](https://rdotsch.github.io/rcicr/reference/generateNoiseImage.md)
  : Generate single noise image based on parameter vector
- [`generateCINoise()`](https://rdotsch.github.io/rcicr/reference/generateCINoise.md)
  : Generate classification image noise pattern based on set of stimuli
  (matrix: trials, parameters), responses (vector), and sinusoid
- [`generateSinusoid()`](https://rdotsch.github.io/rcicr/reference/generateSinusoid.md)
  : Generate single sinusoid patch
- [`generateGabor()`](https://rdotsch.github.io/rcicr/reference/generateGabor.md)
  : Generate single gabor patch
- [`deg2rad()`](https://rdotsch.github.io/rcicr/reference/deg2rad.md) :
  Convert angle in degrees to radians

## Package

- [`rcicr`](https://rdotsch.github.io/rcicr/reference/rcicr-package.md)
  [`rcicr-package`](https://rdotsch.github.io/rcicr/reference/rcicr-package.md)
  : rcicr: Reverse-Correlation Image-Classification Toolbox
