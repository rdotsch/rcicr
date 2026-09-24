# Generate single noise image based on parameter vector

Generate single noise image based on parameter vector

## Usage

``` r
generateNoiseImage(params, p)
```

## Arguments

- params:

  Vector of contrast weights, one per patch in the noise basis.

- p:

  Noise basis, as returned by
  [`generateNoisePattern`](https://rdotsch.github.io/rcicr/reference/generateNoisePattern.md).

## Value

The noise pattern as pixel matrix.

## Examples

``` r
p <- generateNoisePattern(img_size = 32, nscales = 2)

# one contrast weight per patch index
params <- rnorm(max(p$patchIdx))

noise <- generateNoiseImage(params, p)
```
