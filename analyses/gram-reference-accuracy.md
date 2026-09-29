What computing the InfoVal reference through the Gram matrix changes
================

- [The two routes](#the-two-routes)
- [Stimulus sets and classification
  images](#stimulus-sets-and-classification-images)
- [Measurements](#measurements)
- [Against the Monte Carlo error InfoVal already
  carries](#against-the-monte-carlo-error-infoval-already-carries)
- [Time and memory](#time-and-memory)
- [What this shows and what it does
  not](#what-this-shows-and-what-it-does-not)

The InfoVal reference distribution is the norms of classification images
built from random responses. rcicr has always computed them by rendering
every saved stimulus into a pixels-by-trials matrix and multiplying it
by each simulated response vector ([issue
\#354](https://github.com/rdotsch/rcicr/issues/354)). The same norms
follow from the stimulus Gram matrix, a trials-by-trials matrix of dot
products between the noise images, without rendering any image. The two
routes sum in a different order, so their results need not be
bit-identical.

This document measures by how much they differ, and how often the
difference changes a decision at the conventional cut-off of 1.96. It
implements both routes itself, so it measures the same thing whichever
route the package uses.

``` r
library(rcicr)

cutoff <- 1.96
response_seed <- 3
n_cis <- 400
```

## The two routes

A reference CI is `S r / n`: `S` holds the `n` saved noise images as
columns, and `r` the simulated +1/-1 responses. The **rendered route**
computes that product and takes its Frobenius norm, as rcicr always has.
The **Gram route** uses `norm(S r / n) = sqrt(t(r) G r) / n`, with
`G = t(S) S`.

`G` needs no rendered image. A noise image is linear in its parameters,
`s = P x`, where `P` is the basis: one weighted patch per pixel per
layer, divided by the number of layers, as `generateNoiseImage()`
averages them. So `G = X t(P) P t(X)`, with `X` the saved parameter
matrix, and `t(P) P` is computed once from the sparse basis.

Both routes consume the random stream exactly as the package does: one
`runif()` per trial, one iteration after another.

``` r
responses <- function(n, k) matrix(((runif(n * k) > 0.5) * 2) - 1, nrow = n, ncol = k)

rendered_route <- function(params, p, iter, seed) {
  noise <- vapply(seq_len(nrow(params)), function(i) as.vector(generateNoiseImage(params[i, ], p)),
                  numeric(length(p$patches[, , 1])))
  set.seed(seed)
  vapply(seq_len(iter), function(i) {
    r <- responses(nrow(params), 1)
    norm((noise %*% r) / ncol(noise), "f")
  }, numeric(1))
}

basis_gram <- function(p) {
  d <- dim(p$patches)
  npix <- d[1] * d[2]
  basis <- Matrix::sparseMatrix(i = rep(seq_len(npix), d[3]), j = as.vector(p$patchIdx),
                                x = as.vector(p$patches) / d[3], dims = c(npix, max(p$patchIdx)))
  as.matrix(Matrix::crossprod(basis))
}

gram_route <- function(params, p, iter, seed, block = 1000L) {
  G <- params %*% basis_gram(p) %*% t(params)
  n <- nrow(params)
  set.seed(seed)
  out <- numeric(iter)
  done <- 0L
  while (done < iter) {
    k <- min(block, iter - done)
    R <- responses(n, k)
    out[done + seq_len(k)] <- sqrt(pmax(colSums(R * (G %*% R)), 0)) / n
    done <- done + k
  }
  out
}
```

## Stimulus sets and classification images

Each configuration generates its own stimulus set, and every
configuration draws 10,000 references, the package default. The largest
difference is a maximum over draws, and the median and MAD are sampled,
so a smaller draw count would measure less, not the same thing more
cheaply.

The CIs scored against each pair of references come from `generateCI()`.
Each responder follows a fixed template on a random share of trials,
between 0 and 60%, and responds at random on the rest. Their InfoVals
therefore spread across the cut-off.

``` r
quiet <- function(expr) {
  out <- NULL
  suppressWarnings(utils::capture.output(out <- expr))
  out
}

stimulus_set <- function(size, n, noise_type, nscales) {
  path <- tempfile("stim")
  dir.create(path)
  base <- file.path(path, "base.png")
  set.seed(1)
  png::writePNG(matrix(runif(size^2), size, size), base)
  quiet(generateStimuli2IFC(list(base = base), n_trials = n, img_size = size, nscales = nscales,
                            noise_type = noise_type, sigma = if (noise_type == "gabor") 10 else 25,
                            stimulus_path = path, seed = 7, ncores = 1, save_as_png = FALSE))
  list.files(path, pattern = "\\.Rdata$", full.names = TRUE)[1]
}

ci_norms <- function(rdata, n) {
  set.seed(11)
  template <- sample(c(-1, 1), n, replace = TRUE)
  vapply(seq_len(n_cis), function(i) {
    follows <- runif(n) < runif(1, 0, 0.6)
    r <- ifelse(follows, template, sample(c(-1, 1), n, replace = TRUE))
    norm(matrix(quiet(generateCI(seq_len(n), r, "base", rdata, save_as_png = FALSE,
                                 n_cores = 1))$ci), "f")
  }, numeric(1))
}
```

## Measurements

InfoVal is linear in the CI norm: `z = (norm - median) / MAD`. For one
pair of references, the norms whose call at 1.96 differs between them
therefore form a single interval, whose width in InfoVal units is
`dz_at_cutoff`: the difference between the two routes’ InfoVals for a CI
at exactly 1.96. A CI’s call can change only if its InfoVal lies within
that distance of the cut-off. The table also counts the calls that
actually changed among the 400 CIs.

``` r
z_of <- function(norms, reference) (norms - median(reference)) / mad(reference)

measure <- function(size, n, iter, noise_type = "sinusoid", nscales = 5) {
  rdata <- stimulus_set(size, n, noise_type, nscales)
  saved <- new.env()
  load(rdata, envir = saved)
  params <- saved$stimuli_params$base

  rendered_time <- system.time(rendered <- rendered_route(params, saved$p, iter, response_seed))
  gram_time <- system.time(gram <- gram_route(params, saved$p, iter, response_seed))
  package <- quiet(generateReferenceDistribution2IFC(rdata, iter = iter, ncores = 1,
                                                     response_seed = response_seed,
                                                     save_rdata = FALSE))

  norms <- ci_norms(rdata, n)
  z_rendered <- z_of(norms, rendered)
  z_gram <- z_of(norms, gram)
  at_cutoff <- median(rendered) + cutoff * mad(rendered)

  data.frame(
    config = sprintf("%dpx, %d trials, %s, nscales %d", size, n, noise_type, nscales),
    iter = iter,
    identical = identical(rendered, gram),
    rel_norm = max(abs(gram - rendered) / rendered),
    d_median = abs(median(gram) - median(rendered)),
    d_mad = abs(mad(gram) - mad(rendered)),
    max_dz = max(abs(z_gram - z_rendered)),
    dz_at_cutoff = abs(z_of(at_cutoff, gram) - cutoff),
    near_cutoff = sum(abs(z_rendered - cutoff) < 0.5),
    flips = sum((z_rendered > cutoff) != (z_gram > cutoff)),
    speedup = rendered_time[["elapsed"]] / gram_time[["elapsed"]],
    package_vs_routes = max(abs(package - rendered) / rendered, abs(package - gram) / gram),
    rendered_mb = size^2 * n * 8 / 2^20,
    gram_mb = n^2 * 8 / 2^20
  )
}

results <- rbind(
  measure(64, 100, 10000),
  measure(64, 100, 10000, nscales = 3),
  measure(64, 100, 10000, noise_type = "gabor"),
  measure(128, 300, 10000),
  measure(256, 300, 10000),
  measure(256, 770, 10000),
  measure(512, 300, 10000)
)
```

The package’s own reference must agree with both routes, to a relative
1e-12. Otherwise this document is not measuring the package.

``` r
max(results$package_vs_routes)
[1] 4.784245e-14
stopifnot(max(results$package_vs_routes) < 1e-12)
```

``` r
shown <- results[, c("config", "iter", "identical", "rel_norm", "d_median", "d_mad", "max_dz",
                     "dz_at_cutoff", "near_cutoff", "flips")]
shown[, 4:8] <- lapply(shown[, 4:8], function(x) sprintf("%.1e", x))
knitr::kable(shown, row.names = FALSE)
```

| config | iter | identical | rel_norm | d_median | d_mad | max_dz | dz_at_cutoff | near_cutoff | flips |
|:---|---:|:---|:---|:---|:---|:---|:---|---:|---:|
| 64px, 100 trials, sinusoid, nscales 5 | 10000 | FALSE | 3.8e-15 | 3.9e-16 | 1.7e-16 | 5.2e-14 | 4.4e-15 | 13 | 0 |
| 64px, 100 trials, sinusoid, nscales 3 | 10000 | FALSE | 5.4e-15 | 4.4e-16 | 1.3e-15 | 1.3e-13 | 8.7e-14 | 21 | 0 |
| 64px, 100 trials, gabor, nscales 5 | 10000 | FALSE | 3.6e-15 | 5.6e-17 | 2.9e-16 | 1.4e-13 | 8.8e-14 | 23 | 0 |
| 128px, 300 trials, sinusoid, nscales 5 | 10000 | FALSE | 1.1e-14 | 2.2e-16 | 9.1e-16 | 2.1e-13 | 1.1e-13 | 49 | 0 |
| 256px, 300 trials, sinusoid, nscales 5 | 10000 | FALSE | 1.9e-14 | 3.7e-15 | 3.3e-16 | 1.2e-13 | 8.2e-14 | 47 | 0 |
| 256px, 770 trials, sinusoid, nscales 5 | 10000 | FALSE | 1.9e-14 | 1.4e-15 | 1.3e-15 | 2.0e-13 | 5.2e-14 | 28 | 0 |
| 512px, 300 trials, sinusoid, nscales 5 | 10000 | FALSE | 4.8e-14 | 2.5e-14 | 8.1e-15 | 7.9e-13 | 5.6e-13 | 47 | 0 |

`rel_norm` is the largest relative difference in any single norm.
`d_median` and `d_mad` are the absolute differences in the two
statistics InfoVal is built from. `max_dz` is the largest InfoVal
difference among the CIs, and `near_cutoff` how many of them lie within
0.5 of 1.96.

## Against the Monte Carlo error InfoVal already carries

A reference of 10,000 draws is itself a sample, so a second draw of the
same size moves InfoVal. The yardstick is the spread of InfoVal at 1.96
across 40 independent 10,000-draw references for one stimulus set, with
a 200,000-draw reference standing in for the population. It is computed
with the Gram route, whose values the table above shows agree with the
rendered route to rounding.

``` r
rdata <- stimulus_set(512, 300, "sinusoid", 5)
saved <- new.env()
load(rdata, envir = saved)
params <- saved$stimuli_params$base
population <- gram_route(params, saved$p, 200000, 1)
at_cutoff <- median(population) + cutoff * mad(population)
spread <- sd(vapply(1:40, function(i) z_of(at_cutoff, gram_route(params, saved$p, 10000, 100 + i)),
                    numeric(1)))
spread
[1] 0.03135939

largest <- max(results$dz_at_cutoff)
spread / largest
[1] 55558662879
```

## Time and memory

The rendered route holds the pixels-by-trials noise matrix; the Gram
route holds a trials-by-trials matrix, plus the basis cross-product
once. `speedup` is the rendered route’s time over the Gram route’s, both
as implemented above and timed in this run. The rendered route’s
implementation here renders serially, as the package does with
`ncores = 1`.

``` r
cost <- results[, c("config", "iter", "speedup", "rendered_mb", "gram_mb")]
cost$speedup <- round(cost$speedup, 1)
cost[, 4:5] <- lapply(cost[, 4:5], function(x) round(x, 1))
knitr::kable(cost, row.names = FALSE)
```

| config                                 |  iter | speedup | rendered_mb | gram_mb |
|:---------------------------------------|------:|--------:|------------:|--------:|
| 64px, 100 trials, sinusoid, nscales 5  | 10000 |     3.4 |         3.1 |     0.1 |
| 64px, 100 trials, sinusoid, nscales 3  | 10000 |    32.7 |         3.1 |     0.1 |
| 64px, 100 trials, gabor, nscales 5     | 10000 |     3.5 |         3.1 |     0.1 |
| 128px, 300 trials, sinusoid, nscales 5 | 10000 |    13.0 |        37.5 |     0.7 |
| 256px, 300 trials, sinusoid, nscales 5 | 10000 |    43.5 |       150.0 |     0.7 |
| 256px, 770 trials, sinusoid, nscales 5 | 10000 |    41.0 |       385.0 |     4.5 |
| 512px, 300 trials, sinusoid, nscales 5 | 10000 |   106.6 |       600.0 |     0.7 |

At the package defaults, 512 pixels and 770 trials, the rendered noise
matrix alone is 1.5 GB. The Gram matrix is 4.5 MB, and the basis
cross-product for five scales is 128 MB.

## What this shows and what it does not

None of the tested configurations gave bit-identical references. The
largest relative difference in any single norm is 4.8e-14. The largest
InfoVal difference among 2800 CIs is 7.9e-13, at 512px, 300 trials,
sinusoid, nscales 5. At 1.96 the largest difference is 5.6e-13: under
these references, a CI’s call can change only if its InfoVal lies within
that distance of the cut-off. None of the 228 CIs within 0.5 of 1.96
changed its call (0 of 2800 overall).

That distance is 5.6e+10 times smaller than the Monte Carlo spread of
InfoVal at 1.96 for a 10,000-draw reference, 0.031. So the band of
InfoVals whose call the route can change is that many times narrower
than the band a fresh draw of the reference already moves.

This document does not measure:

- the 512-pixel, 770-trial default with the rendered route, which needs
  more memory than the machine that knitted it had; the table stops at
  512 pixels and 300 trials;
- references already stored in `.Rdata` files, which are reused as
  stored and do not change;
- optimised BLAS libraries, which reorder sums in either route and can
  differ from these values by similar amounts.
