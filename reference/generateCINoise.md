# Generate classification image noise pattern based on set of stimuli (matrix: trials, parameters), responses (vector), and sinusoid

Generate classification image noise pattern based on set of stimuli
(matrix: trials, parameters), responses (vector), and sinusoid

## Usage

``` r
generateCINoise(stimuli, responses, p)
```

## Arguments

- stimuli:

  Matrix with one parameter row per response. A parameter vector is also
  accepted for a single trial with exactly one response.

- responses:

  Vector with the response to each trial: 1 where the participant chose
  the original, -1 where they chose the inverted image. Other values,
  such as ratings on a scale, weight the trials accordingly.

- p:

  Noise basis, as returned by
  [`generateNoisePattern`](https://rdotsch.github.io/rcicr/reference/generateNoisePattern.md).

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
