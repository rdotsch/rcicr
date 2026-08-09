# Generate classification image noise pattern based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid

Generate classification image noise pattern based on set of stimuli
(matrix: trials, parameters), responses (vector), and sinusoid

## Usage

``` r
generateCINoise(stimuli, responses, p)
```

## Arguments

- stimuli:

  Matrix with one row per trial, each row containing the 4092 parameters
  for the original stimulus.

- responses:

  Vector containing the response to each trial (1 if participant
  selected original, -1 if participant selected inverted; this can be
  changed into a scale).

- p:

  3D patch matrix (generated using
  [`generateNoisePattern()`](https://rdotsch.github.io/rcicr/reference/generateNoisePattern.md)).

## Value

The classification image as pixel matrix.

## Examples

``` r
p <- generateNoisePattern(img_size = 32, nscales = 1)
nparams <- max(p$patchIdx)

# two trials, one chosen original (1) and one inverted (-1)
stimuli <- matrix(runif(2 * nparams, -1, 1), nrow = 2)
responses <- c(1, -1)

ci <- generateCINoise(stimuli, responses, p)
```
