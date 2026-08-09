# Generate single noise image based on parameter vector

Generate single noise image based on parameter vector

## Usage

``` r
generateNoiseImage(params, p)
```

## Arguments

- params:

  Vector with each value specifying the contrast of each patch in noise.

- p:

  3D patch matrix (generated using
  [`generateNoisePattern()`](https://rdotsch.github.io/rcicr/reference/generateNoisePattern.md)).

## Value

The noise pattern as pixel matrix.

## Examples

``` r
p <- generateNoisePattern(img_size = 32, nscales = 2)

# one contrast weight per patch index
params <- rnorm(max(p$patchIdx))

noise <- generateNoiseImage(params, p)
```
