# Plan: drop the matlab and scales imports (#208)

## Where they are used

On `main` at bdafa78, `grep -n "matlab::\|scales::" R/*.R` finds 7 calls in 4 files:

| call | file |
|---|---|
| `matlab::zeros(c(img_size, img_size, nrPatches))` ×2 | `generateNoisePattern.R` |
| `matlab::repmat(p, scale)` | `generateNoisePattern.R` |
| `matlab::zeros(nrep, 2)` | `simulateNoiseIntensities.R` |
| `matlab::repmat(matlab::linspace(0, cycles, img_size), img_size, 1)` | `generateSinusoid.R` |
| `matlab::meshgrid(x0, x0)` | `generateGabor.R` |
| `scales::rescale(1:img_size, to = c(-.5, .5))` | `generateGabor.R` |

The three-argument `meshgrid` #208 called the risky site is already gone (#386, checked identical).

## The replacements, from the packages' own code

Read from the installed `matlab` 1.0.4 and `scales` sources:

- `zeros(d)` is `array(0, d)`.
- `repmat(A, n)` on a matrix is `kronecker(array(1, c(n, n)), A)`, or `A` itself when `n == 1`.
- `repmat(v, n, 1)` on a vector fills rows, so it becomes `matrix(v, n, n, byrow = TRUE)`, or `v` itself when `n == 1`.
- `linspace(a, b, n)` is `seq(a, b, length.out = n)`, but returns `b` when `n < 2`.
- `meshgrid(x, x)` is `matrix(x, n, n, byrow = TRUE)` and `matrix(x, n, n)`.
- `rescale(1:n, to = c(-.5, .5))` is `(seq_len(n) - 1) / (n - 1) - 0.5`, but returns 0 when `n == 1`.

The `n == 1` and `n < 2` branches are kept on purpose. A patch of size 1 occurs when `img_size == 2^(nscales - 1)`; `seq()` alone would give 0 there instead of `b`, and the rescale would give `NaN` instead of 0.

## Verified before writing the plan

The candidate versions were compared with `identical()` against `generateSinusoid()` and `generateGabor()` on `main`: every `img_size` from 1 to 64 plus 128, 256 and 512, all six orientations, both phases, two `cycles` and two `sigma` values. Also compared: `repmat` at patch sizes 1, 4, 16 and 64 tiled 1, 2, 4, 8 and 16 times, and `zeros` in both call forms. **3,240 cases, 0 differ.**

## The change

1. The six replacements above. `matlab` and `scales` come out of `DESCRIPTION` Imports. `tools/setup-container-r.sh` keeps installing them, because the gate's reference versions import them.
2. `test-generateGabor.R` built its expected value with `scales::rescale()` and `matlab::meshgrid()`, which a test may no longer use. It gets an oracle from the definition instead: a Gaussian over coordinates from -0.5 to 0.5, `seq(-0.5, 0.5, length.out = n)`, compared with `expect_equal()`, so it does not mirror the implementation. No other test changes.
3. New tests pin the size-1 patch of `generateSinusoid()` and `generateGabor()` to its value on `main`, since those are the branches most likely to be simplified away later.
4. NEWS, under "Performance and dependencies": `matlab` and `scales` are no longer dependencies, and every noise pattern is unchanged.

## Verification, and the step most likely to fail

- **Full** release gate against `main`: `0 deviations`. It compares `patchIdx` exactly, and the patches and every stimulus PNG hash.
- The golden master (`test-regression-baseline.R`) unchanged and passing. Full suite, `R CMD check` in CI.
- **Most likely to fail:** the `generateNoisePattern()` tiling at `nscales = 1`, where `scale == 1` and the `n == 1` branch returns `p` as is. The gate's `sinusoid-128-nscales3` and default configurations do not reach `scale == 1` alone, so a test compares `generateNoisePattern(nscales = 1)`'s `patches` with its value on `main`.
