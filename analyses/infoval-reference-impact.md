What the wrong-base InfoVal reference does to a conclusion
================

- [The reference distribution ignores which base you
  score](#the-reference-distribution-ignores-which-base-you-score)
- [Base 2’s own reference](#base-2s-own-reference)
- [Scatter and the estimated directional
  component](#scatter-and-the-estimated-directional-component)
- [How large the scatter is](#how-large-the-scatter-is)
- [Trial count and resolution](#trial-count-and-resolution)
- [What it takes to reach a different
  conclusion](#what-it-takes-to-reach-a-different-conclusion)
- [What is not affected](#what-is-not-affected)
- [Limits](#limits)

When stimuli are generated for several base images with
`use_same_parameters = FALSE`, each base gets its own parameter matrix,
but `generateReferenceDistribution2IFC()` reconstructs the stimuli
without forwarding that setting. The null it builds is therefore always
the first base image’s ([issue
\#299](https://github.com/rdotsch/rcicr/issues/299)), and every base
after the first is scored against another base’s noise.

This document characterizes the discrepancy using synthetic stimulus
sets. A reference from the same generator is not necessarily calibrated
for another base’s saved noise. The measured shifts vary across pairs,
with a small estimated directional component in the configurations
examined. The median component decreases across the sampled trial
counts; the scale component’s trend is unresolved. Neither the
prevalence of affected studies nor their changes in significance calls
is measured here.

Knitting runs the full sweep of 100 base image pairs, including ten at
512 pixels. Every pair uses three 10,000-draw references; the fifty 64px
pairs at 300 or 770 trials also use four 100,000-draw references.
Runtime depends on the available hardware.

``` r
library(rcicr)

img_size <- 64
nscales <- 5
iter <- 10000
response_seed <- 97
cutoff <- 1.96
```

``` r
scratch <- tempfile("infoval-impact")
dir.create(scratch)

base_image <- function(name, seed, size) {
  f <- file.path(scratch, sprintf("%s-%d.png", name, size))
  set.seed(seed)
  png::writePNG(matrix(runif(size^2), size, size), f)
  f
}

# capture.output() here and below only keeps progress bars out of the rendered
# document; the package's functions are called exactly as a user would call them.
# ?generateReferenceDistribution2IFC guarantees the stimuli, and so everything
# measured here, do not depend on ncores. mclapply() below rejects more than one
# core on Windows, which knits there serially rather than not at all.
cores <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)

generate <- function(bases, n_trials, seed, size, same_parameters) {
  path <- tempfile("stim", tmpdir = scratch)
  dir.create(path)
  invisible(capture.output(
    generateStimuli2IFC(bases, n_trials = n_trials, img_size = size, stimulus_path = path,
                        seed = seed, use_same_parameters = same_parameters, nscales = nscales,
                        ncores = cores, save_as_png = FALSE)
  ))
  list.files(path, pattern = "\\.Rdata$", full.names = TRUE)[1]
}
```

The base images are synthetic noise rather than faces. Nothing here
depends on their content: the classification image and the reference
norms are built from the noise alone, and the base image is added only
for display.

## The reference distribution ignores which base you score

Two base images, independent parameters, against the same stimulus set
containing only the first of them. If the reference carried any
information about which base is being scored, these two runs would
differ.

``` r
a <- base_image("a", 11, img_size)
b <- base_image("b", 12, img_size)

two_bases <- generate(list(a = a, b = b), n_trials = 24, seed = 31, size = img_size,
                      same_parameters = FALSE)
first_only <- generate(list(a = a), n_trials = 24, seed = 31, size = img_size,
                       same_parameters = TRUE)

# iter far below the recommended 10,000 here: this compares two runs against each
# other, not against a cut-off, and the warning is the function's own.
demo_reference <- function(rdata) {
  invisible(capture.output(norms <- suppressWarnings(generateReferenceDistribution2IFC(
    rdata, iter = 200, ncores = 1, response_seed = response_seed, save_rdata = FALSE))))
  norms
}

reference_two <- demo_reference(two_bases)
reference_one <- demo_reference(first_only)

saved <- new.env()
load(two_bases, envir = saved)

max(abs(saved$stimuli_params[["a"]] - saved$stimuli_params[["b"]]))
[1] 1.988195
max(abs(reference_two - reference_one))
[1] 0
```

The two bases hold genuinely different parameters, and the reference
distribution is bit-identical to the one for a file that does not
contain the second base at all.

## Base 2’s own reference

The correction requires the null the second base should have had: its
noise built from the parameters actually saved for it.
`generateNoiseImage()` does the reconstruction, and the loop below is
the arithmetic of `generateReferenceDistribution2IFC()`, batched so that
the response draws are consumed in the same order.

``` r
noise_matrix <- function(params, p, size) {
  matrix(unlist(parallel::mclapply(seq_len(nrow(params)), mc.cores = cores,
                function(trial) as.vector(generateNoiseImage(params[trial, ], p)))),
         nrow = size^2)
}

# A block holds one column per iteration, so its width comes from the image size,
# keeping the intermediate near 200 MB at any resolution.
reference_norms <- function(noise, iter, seed, block = max(50L, as.integer(2.5e7 / nrow(noise)))) {
  set.seed(seed)
  n_trials <- ncol(noise)
  out <- numeric(iter)
  done <- 0L
  while (done < iter) {
    k <- min(block, iter - done)
    responses <- matrix(((runif(n_trials * k) > 0.5) * 2) - 1, nrow = n_trials, ncol = k)
    out[(done + 1):(done + k)] <- sqrt(colSums(((noise %*% responses) / n_trials)^2))
    done <- done + k
  }
  out
}
```

The first pair of each configuration below 512px checks this
reimplementation against `generateReferenceDistribution2IFC()` on the
first base. At 512px that comparison is skipped because of the repeated
matrix conversion in \#306. Agreement in the checked configurations
supports the reimplementation and the first-base reconstruction; it is
not a direct validation at 512px.

``` r
measure_pair <- function(n_trials, seed, size, validate) {
  rdata <- generate(list(a = base_image("a", 11, size), b = base_image("b", 12, size)),
                    n_trials = n_trials, seed = seed, size = size, same_parameters = FALSE)
  saved <- new.env()
  load(rdata, envir = saved)

  noise_a <- noise_matrix(saved$stimuli_params[["a"]], saved$p, size)
  noise_b <- noise_matrix(saved$stimuli_params[["b"]], saved$p, size)

  scored_against <- reference_norms(noise_a, iter, response_seed)   # what the bug uses
  own <- reference_norms(noise_b, iter, response_seed)              # what base 2 should have
  # Same base, an independent draw of the simulated responses: the Monte Carlo
  # wobble InfoVal already carries at this iter, as a yardstick for the shift.
  monte_carlo <- reference_norms(noise_a, iter, response_seed + 1L)
  # Finite-reference proxy for the population, restricted to the configurations
  # used by the comparison and trend diagnostics below.
  decompose <- size == 64 && n_trials %in% c(300, 770)
  precise <- if (decompose) reference_norms(noise_a, 10L * iter, response_seed + 2L) else NULL
  # Match response draws across bases to retain the main comparison's coupling.
  # Finite-reference error remains at this larger iteration count.
  precise_b <- if (decompose) reference_norms(noise_b, 10L * iter, response_seed + 2L) else NULL
  # Replicate the same between-base estimator with a second shared response draw.
  # Differencing removes its common target and bias, not just its population target.
  precise_alt <- if (decompose) reference_norms(noise_a, 10L * iter, response_seed + 5L) else NULL
  precise_b_alt <- if (decompose) reference_norms(noise_b, 10L * iter, response_seed + 5L) else NULL

  agreement <- NA_real_
  if (validate) {
    invisible(capture.output(package_reference <- generateReferenceDistribution2IFC(
      rdata, iter = iter, ncores = 1, response_seed = response_seed, save_rdata = FALSE)))
    agreement <- max(abs(scored_against - package_reference))
    stopifnot(agreement < 1e-12)
  }

  data.frame(
    n_trials = n_trials, img_size = size, seed = seed, agreement = agreement,
    med_own = median(own), mad_own = mad(own),
    med_scored = median(scored_against), mad_scored = mad(scored_against),
    d0 = (median(own) - median(scored_against)) / mad(scored_against),
    rho = mad(own) / mad(scored_against) - 1,
    mc_shift = (median(monte_carlo) - median(scored_against)) / mad(scored_against),
    mc_rho = mad(monte_carlo) / mad(scored_against) - 1,
    # Oriented like d0 and rho above: the reference that defines the InfoVal on
    # top, the one actually used underneath, so the two are comparable at a z
    # meaning the same thing.
    single_shift = if (is.null(precise)) NA_real_ else
      (median(precise) - median(scored_against)) / mad(scored_against),
    single_rho = if (is.null(precise)) NA_real_ else mad(precise) / mad(scored_against) - 1,
    between_d0 = if (is.null(precise_b)) NA_real_ else
      (median(precise_b) - median(precise)) / mad(precise),
    between_rho = if (is.null(precise_b)) NA_real_ else mad(precise_b) / mad(precise) - 1,
    between_d0_alt = if (is.null(precise_alt)) NA_real_ else
      (median(precise_b_alt) - median(precise_alt)) / mad(precise_alt),
    between_rho_alt = if (is.null(precise_alt)) NA_real_ else
      mad(precise_b_alt) / mad(precise_alt) - 1
  )
}

# Validated on the first pair of each configuration, except at 512px, where
# generateReferenceDistribution2IFC() rebuilds a matrix from its stimulus frame on
# every one of the 10,000 iterations (#306). That frame carries 64 times the pixels
# it does at 64px, where the check is already the slowest thing in this document,
# so the cost there is prohibitive. The reimplementation does not vary with
# resolution.
sweep <- function(n_trials, seeds, size = img_size) {
  do.call(rbind, Map(function(seed, validate) measure_pair(n_trials, seed, size, validate),
                     seeds, seq_along(seeds) == 1 & size < 512))
}
```

`d0` is the shift in the reference median, expressed in the units
InfoVal is reported in, and `rho` is the error in its scale. An InfoVal
that should have been `z` is reported as `z + d0 + z * rho`, so both
terms reach a threshold call — `rho` multiplied by the threshold — and
the scale term grows in weight as the InfoVal does.

``` r
results <- rbind(
  sweep(770, 1:30),
  sweep(300, 101:120),
  sweep(100, 201:220),
  sweep(300, 301:312, size = 128),
  sweep(770, 501:508, size = 128),
  sweep(300, 401:410, size = 512)
)
```

``` r
max(results$agreement, na.rm = TRUE)
[1] 2.609024e-15
```

``` r
# An InfoVal that should read z is reported as z + shift(). The Monte Carlo
# yardstick has to carry the same scale term, or it is not the same quantity.
shift <- function(rows, z) rows$d0 + z * rows$rho
monte_carlo <- function(rows, z) rows$mc_shift + z * rows$mc_rho

centre <- function(x) {
  se <- sd(x) / sqrt(length(x))
  c(mean = mean(x), se = se, ses_from_zero = abs(mean(x)) / se)
}

at_64px <- subset(results, img_size == 64)
```

Standard errors across base pairs describe variation conditional on the
selected response streams. The same response seeds are reused across
pairs, and response seeds are not independently resampled for each row.
These standard errors do not include all uncertainty over response
simulations.

## Scatter and the estimated directional component

``` r
# 770 trials at 64 pixels. Not the package's default configuration, whose
# 512 pixels are measured separately below and at fewer trials.
at_770 <- subset(results, n_trials == 770 & img_size == 64)

knitr::kable(round(rbind(
  "reference median only (d0)" = centre(at_770$d0),
  "scale only (rho)" = centre(at_770$rho),
  "full shift at the cut-off" = centre(shift(at_770, cutoff)),
  "full shift at z = 10" = centre(shift(at_770, 10))
), 4))
```

|                            |    mean |     se | ses_from_zero |
|:---------------------------|--------:|-------:|--------------:|
| reference median only (d0) | -0.0066 | 0.0065 |        1.0074 |
| scale only (rho)           | -0.0020 | 0.0038 |        0.5409 |
| full shift at the cut-off  | -0.0106 | 0.0114 |        0.9282 |
| full shift at z = 10       | -0.0269 | 0.0403 |        0.6684 |

Base 2’s parameters come from the same generator as base 1’s, but its
reference should be conditional on its own saved noise. Across 30 base
image pairs at 770 trials and 64 pixels the full shift at the cut-off
averages -0.011, which is 0.9 standard errors from zero against a
scatter 6 times the size of that mean.

The scale ratio need not be centred at zero. `rho` is the ratio
`mad(own) / mad(scored_against)` less one, and for exchangeable positive
scales `B/A + A/B >= 2` forces `E[B/A] >= 1`, strictly when they differ
with positive probability (assuming finite expectations). This
inequality does not determine the sign of the median term: `d0` also
depends on how the two medians relate to the denominator MAD. Shared
response draws can couple the two reference estimates, so cross-base
dependence cannot simply be dropped.

Averaging each pair’s two labellings gives a symmetric estimator of the
directional component under exchangeability. It uses both possible
orderings of the same pair, with the standard error computed across
pairs rather than treating the two orderings as independent
observations.

``` r
directional <- function(rows, z) {
  d0 <- ((rows$med_own - rows$med_scored) / rows$mad_scored +
         (rows$med_scored - rows$med_own) / rows$mad_own) / 2
  ratio <- 1 + rows$rho
  d0 + z * ((ratio + 1 / ratio) / 2 - 1)
}

knitr::kable(signif(do.call(rbind, lapply(split(at_64px, at_64px$n_trials), function(rows)
  c("median term" = centre(directional(rows, 0))[c("mean", "se")],
    "at the cut-off" = centre(directional(rows, cutoff))[c("mean", "se")],
    "at z = 10" = centre(directional(rows, 10))[c("mean", "se")]))), 3))
```

|     | median term.mean | median term.se | at the cut-off.mean | at the cut-off.se | at z = 10.mean | at z = 10.se |
|:----|-----------------:|---------------:|--------------------:|------------------:|---------------:|-------------:|
| 100 |         1.37e-03 |       5.11e-04 |            0.002380 |          0.000859 |        0.00650 |     0.002370 |
| 300 |         4.03e-05 |       1.04e-04 |            0.000452 |          0.000131 |        0.00214 |     0.000475 |
| 770 |         1.29e-04 |       4.55e-05 |            0.000531 |          0.000142 |        0.00218 |     0.000705 |

At 770 trials and 64 pixels the estimated directional part of the shift
is 5.3e-04 at the cut-off and 0.002 at an InfoVal of 10 — around 1% of
the mean absolute shift at that value. These are estimates for the
sampled configurations and the specified InfoVals, not evidence that
directional effects are absent or irrelevant in every study. The tables
show their uncertainty alongside their magnitude.

## How large the scatter is

``` r
summary_row <- function(rows) {
  dz <- shift(rows, cutoff)
  data.frame(
    n_trials = rows$n_trials[1], img_size = rows$img_size[1], pairs = nrow(rows),
    mean_abs_d0 = round(mean(abs(rows$d0)), 4),
    mean_abs_dz = round(mean(abs(dz)), 4),
    q90_abs_dz = round(unname(quantile(abs(dz), 0.9)), 4),
    max_abs_dz = round(max(abs(dz)), 4),
    monte_carlo = round(mean(abs(monte_carlo(rows, cutoff))), 4)
  )
}

by_config <- split(results, list(results$n_trials, results$img_size), drop = TRUE)
knitr::kable(do.call(rbind, lapply(by_config, summary_row)), row.names = FALSE)
```

| n_trials | img_size | pairs | mean_abs_d0 | mean_abs_dz | q90_abs_dz | max_abs_dz | monte_carlo |
|---------:|---------:|------:|------------:|------------:|-----------:|-----------:|------------:|
|      100 |       64 |    20 |      0.0788 |      0.1206 |     0.3166 |     0.3473 |      0.0340 |
|      300 |       64 |    20 |      0.0449 |      0.0577 |     0.1033 |     0.1197 |      0.0265 |
|      770 |       64 |    30 |      0.0299 |      0.0531 |     0.1003 |     0.1216 |      0.0332 |
|      300 |      128 |    12 |      0.0512 |      0.0729 |     0.1263 |     0.1318 |      0.0211 |
|      770 |      128 |     8 |      0.0268 |      0.0553 |     0.0936 |     0.1007 |      0.0285 |
|      300 |      512 |    10 |      0.0530 |      0.0812 |     0.1523 |     0.2825 |      0.0261 |

`monte_carlo` is the same quantity for a reference built from the same
base image with a different draw of the simulated responses. It is
rerun-to-rerun variability, the difference between two references of
`iter` = 10,000 draws each, so it carries the sampling error of both.
Error against an exact population reference is a different quantity. The
next comparison approximates it with a 100,000-draw reference, which
still has estimation error; the table’s population label denotes this
finite-reference proxy.

``` r
single_error <- function(rows, z) rows$single_shift + z * rows$single_rho

knitr::kable(signif(rbind(
  `the bug` = c(mean = mean(abs(shift(at_770, cutoff)))),
  `rerun to rerun` = c(mean = mean(abs(monte_carlo(at_770, cutoff)))),
  `one reference against the population` = c(mean = mean(abs(single_error(at_770, cutoff))))
), 3))
```

|                                      |   mean |
|:-------------------------------------|-------:|
| the bug                              | 0.0531 |
| rerun to rerun                       | 0.0332 |
| one reference against the population | 0.0283 |

So at 770 trials and 64 pixels the bug moves an InfoVal about 1.9 times
as far as the discrepancy against that finite-reference proxy, and about
1.6 times the rerun-to-rerun discrepancy. Neither comparison measures
the exact population error.

The two terms contribute about equally at the cut-off: the median term
averages 0.0299 and the scale term 0.0292 once multiplied by 1.96. The
scale term’s closeness to its own Monte Carlo control (0.0149 against
0.0134) does not establish how much comes from estimation error and how
much from a difference in the underlying reference spreads. The
diagnostics below explore this without assigning component shares.

``` r
at_z <- function(z, f) vapply(z, function(zz) mean(abs(f(at_770, zz))), numeric(1))

z <- c(2, 5, 10)
knitr::kable(data.frame(
  z = z,
  mean_abs_dz = round(at_z(z, shift), 4),
  as_share_of_z = round(at_z(z, shift) / z, 4),
  monte_carlo = round(at_z(z, monte_carlo), 4),
  ratio = round(at_z(z, shift) / at_z(z, monte_carlo), 2)
))
```

|   z | mean_abs_dz | as_share_of_z | monte_carlo | ratio |
|----:|------------:|--------------:|------------:|------:|
|   2 |      0.0537 |        0.0268 |      0.0337 |  1.59 |
|   5 |      0.0963 |        0.0193 |      0.0736 |  1.31 |
|  10 |      0.1694 |        0.0169 |      0.1407 |  1.20 |

In these sampled rows the ratio narrows at the displayed larger
InfoVals. Both discrepancies include a scale component multiplied by z,
but those components are different random quantities. Similarity of
their mean magnitudes is not an equivalence test.

## Trial count and resolution

``` r
by_trials <- do.call(rbind, lapply(split(at_64px, at_64px$n_trials), summary_row))
reference <- by_trials[by_trials$n_trials == 100, "mean_abs_d0"]
knitr::kable(data.frame(
  n_trials = by_trials$n_trials,
  observed = by_trials$mean_abs_d0,
  predicted_from_100 = round(reference * sqrt(100 / by_trials$n_trials), 4)
), row.names = FALSE)
```

| n_trials | observed | predicted_from_100 |
|---------:|---------:|-------------------:|
|      100 |   0.0788 |             0.0788 |
|      300 |   0.0449 |             0.0455 |
|      770 |   0.0299 |             0.0284 |

The sampled mean absolute median component is close to the displayed
inverse-square-root curve. This is a descriptive comparison at 64px, not
an established law for every trial count or resolution.

At 300 trials, the sampled mean absolute shift increases across the
three resolutions. The uncertainty below does not establish a resolution
effect.

``` r
by_size <- function(rows) {
  d0 <- abs(rows$d0); dz <- abs(shift(rows, cutoff))
  c(pairs = nrow(rows),
    mean_abs_d0 = mean(d0), se_d0 = sd(d0) / sqrt(nrow(rows)),
    mean_abs_dz = mean(dz), se_dz = sd(dz) / sqrt(nrow(rows)))
}

at_300 <- subset(results, n_trials == 300)
sizes <- do.call(rbind, lapply(split(at_300, at_300$img_size), by_size))
knitr::kable(signif(sizes, 3))
```

|     | pairs | mean_abs_d0 |   se_d0 | mean_abs_dz |  se_dz |
|:----|------:|------------:|--------:|------------:|-------:|
| 64  |    20 |      0.0449 | 0.00499 |      0.0577 | 0.0074 |
| 128 |    12 |      0.0512 | 0.00943 |      0.0729 | 0.0122 |
| 512 |    10 |      0.0530 | 0.01540 |      0.0812 | 0.0262 |

``` r
gap <- function(small, large, f) {
  x <- abs(f(subset(results, n_trials == 300 & img_size == small)))
  y <- abs(f(subset(results, n_trials == 300 & img_size == large)))
  se <- sqrt(var(x) / length(x) + var(y) / length(y))
  c(ratio = mean(y) / mean(x), difference = mean(y) - mean(x), se = se,
    ses_from_zero = abs(mean(y) - mean(x)) / se)
}

knitr::kable(signif(rbind(
  `d0, 512px vs 64px` = gap(64, 512, function(rows) rows$d0),
  `shift at the cut-off, 512px vs 64px` = gap(64, 512, function(rows) shift(rows, cutoff))
), 3))
```

|                                     | ratio | difference |     se | ses_from_zero |
|:------------------------------------|------:|-----------:|-------:|--------------:|
| d0, 512px vs 64px                   |  1.18 |    0.00819 | 0.0162 |         0.506 |
| shift at the cut-off, 512px vs 64px |  1.41 |    0.02350 | 0.0273 |         0.861 |

At the default resolution the shift at the cut-off is 1.41 times its
value at 64px — but on 10 pairs that is 0.9 standard errors for the
difference of means, so a resolution effect remains unresolved. The
ratio is a point estimate, not a measured resolution multiplier that can
be applied to other configurations.

Read the trial-count rows above as measured at 64px, then, and this
512px row as the one that applies to a study at the package’s own
resolution: around 0.08 in z at 300 trials, against 0.06 at 64px.

The specific combination of 512 pixels and 770 trials is not measured. A
separable model would predict the same trial-count ratio at each
resolution; the next table describes those ratios at 64px and 128px but
does not establish separability at 512px.

``` r
trial_ratio <- function(size) {
  at <- function(n) abs(shift(subset(results, n_trials == n & img_size == size), cutoff))
  x <- at(300); y <- at(770)
  c(pairs_300 = length(x), pairs_770 = length(y),
    ratio_770_over_300 = mean(y) / mean(x),
    se = (mean(y) / mean(x)) * sqrt(var(x) / (length(x) * mean(x)^2) +
                                    var(y) / (length(y) * mean(y)^2)))
}

knitr::kable(signif(rbind(`64px` = trial_ratio(64), `128px` = trial_ratio(128),
                          `sqrt(300/770)` = c(NA, NA, sqrt(300 / 770), NA)), 3))
```

|               | pairs_300 | pairs_770 | ratio_770_over_300 |    se |
|:--------------|----------:|----------:|-------------------:|------:|
| 64px          |        20 |        30 |              0.921 | 0.158 |
| 128px         |        12 |         8 |              0.758 | 0.226 |
| sqrt(300/770) |        NA |        NA |              0.624 |    NA |

The point estimates for the full shift differ from the
inverse-square-root prediction. The shift contains both median and scale
components, each with finite-reference estimation error. Comparing the
tenfold references gives further diagnostics, but neither establishes a
scale trend nor justifies holding it fixed when extrapolating.

``` r
decomposed <- subset(results, !is.na(between_rho))
component <- function(x) c(mean = mean(abs(x)), se = sd(abs(x)) / sqrt(length(x)))

knitr::kable(signif(do.call(rbind, lapply(split(decomposed, decomposed$n_trials), function(rows)
  c(pairs = nrow(rows),
    median_term = component(rows$d0), between_base_median = component(rows$between_d0),
    scale_term = component(rows$rho), between_base_scale = component(rows$between_rho)))), 3))
```

|     | pairs | median_term.mean | median_term.se | between_base_median.mean | between_base_median.se | scale_term.mean | scale_term.se | between_base_scale.mean | between_base_scale.se |
|:----|------:|-----------------:|---------------:|-------------------------:|-----------------------:|----------------:|--------------:|------------------------:|----------------------:|
| 300 |    20 |           0.0449 |        0.00499 |                   0.0380 |                0.00485 |          0.0171 |       0.00253 |                 0.01030 |               0.00209 |
| 770 |    30 |           0.0299 |        0.00367 |                   0.0254 |                0.00436 |          0.0149 |       0.00257 |                 0.00828 |               0.00112 |

Whether either between-base component falls with trials is what the
combined estimate turns on, and twenty and thirty pairs do not settle it
from point estimates alone.

``` r
trend <- function(field) {
  x <- abs(decomposed[[field]][decomposed$n_trials == 300])
  y <- abs(decomposed[[field]][decomposed$n_trials == 770])
  difference <- mean(y) - mean(x)
  se <- sqrt(var(x) / length(x) + var(y) / length(y))
  c(at_300 = mean(x), at_770 = mean(y), difference = difference, se = se,
    ses_from_zero = abs(difference) / se)
}

knitr::kable(signif(rbind(`between-base median` = trend("between_d0"),
                          `between-base scale` = trend("between_rho")), 3))
```

|                     | at_300 |  at_770 | difference |      se | ses_from_zero |
|:--------------------|-------:|--------:|-----------:|--------:|--------------:|
| between-base median | 0.0380 | 0.02540 |   -0.01270 | 0.00653 |         1.940 |
| between-base scale  | 0.0103 | 0.00828 |   -0.00197 | 0.00237 |         0.831 |

A positive mean absolute difference alone cannot establish unequal
population spreads: finite references differ even under equality. Each
pair therefore repeats the between-base estimate with a second response
seed. Conditional on the saved bases, differencing the replicates
removes their common target and any common bias. It probes
response-sampling variability, not bias relative to the population MAD
ratio.

``` r
precise_770 <- decomposed[decomposed$n_trials == 770, ]

# Two shared-draw replicates of the same between-base quantity. Their difference
# carries response-sampling error alone; common finite-sample bias survives it,
# so these are magnitudes to report rather than terms of a test.
estimate <- (precise_770$between_rho + precise_770$between_rho_alt) / 2
error <- precise_770$between_rho - precise_770$between_rho_alt

c(between_base = mean(abs(estimate)), replicate_error_sd = sd(error) / sqrt(2),
  se_of_mean = sd(abs(estimate)) / sqrt(nrow(precise_770)))
      between_base replicate_error_sd         se_of_mean 
       0.007990065        0.004425206        0.001126322 
```

The replicate spread is response-sampling variability at this iteration
count, not a calibration of the magnitude beside it: mean absolute
values and standard deviations are different summaries, and the common
finite-sample bias in a MAD ratio is invisible to a difference of
replicates. Sharing response seeds across base pairs also means these
across-pair summaries do not capture all response-seed uncertainty.

Whether the two bases’ population MAD ratios differ at all, and how much
of the scale term is estimation error, are therefore left open here.
Both need an inferential model for the replicate error rather than a
comparison of summary magnitudes. No conclusion below rests on either.

In these samples the median component’s mean absolute value drops by 1.9
standard errors, close to what the inverse-square-root model predicts —
0.0237 against 0.0254 observed. The scale one has a lower point estimate
at the higher count, but at 0.8 standard errors these samples do not
separate that from no change, so whether it falls at all is open here.

The following extrapolation assumes inverse-square-root scaling for the
median component and no change in the scale component of each sampled
512px pair. Both are modelling assumptions at the unmeasured trial
count. An unresolved scale trend does not support a zero trend, and the
direction of error is unknown because signed components can cancel.

``` r
components <- function(rows) c(median_term = mean(abs(rows$d0)), scale_term = mean(abs(rows$rho)))
knitr::kable(signif(do.call(rbind, lapply(split(at_64px, at_64px$n_trials), components)), 3))
```

|     | median_term | scale_term |
|:----|------------:|-----------:|
| 100 |      0.0788 |     0.0240 |
| 300 |      0.0449 |     0.0171 |
| 770 |      0.0299 |     0.0149 |

``` r

at_512 <- subset(results, img_size == 512)
trial_law <- sqrt(300 / 770)
estimate <- mean(abs(at_512$d0 * trial_law + cutoff * at_512$rho))
# Signed terms can cancel, so the estimate is not one-sided. Adding the two
# magnitudes instead cannot: per pair it bounds that pair, and averaged it bounds
# the mean rather than any individual study.
per_pair_bound <- abs(at_512$d0) * trial_law + cutoff * abs(at_512$rho)
bound <- mean(per_pair_bound)
c(measured_512px_300_trials = mean(abs(shift(at_512, cutoff))),
  estimated_512px_770_trials = estimate, mean_bound = bound,
  largest_pair_bound = max(per_pair_bound))
 measured_512px_300_trials estimated_512px_770_trials                 mean_bound 
                0.08117167                 0.06123524                 0.06616231 
        largest_pair_bound 
                0.22161474 
```

Under these assumptions the estimate scales the median term and holds
the scale term fixed. That is a model rather than a bound on real
studies: both terms are signed and can cancel, so holding the scale term
at its 300-trial value can as easily lower the result as raise it. On
that model a study at 512 pixels and 770 trials carries roughly 0.06 in
z, against 0.08 measured at 512 pixels and 300 trials — a modelled
sample mean rather than a measurement at that configuration. Adding the
two magnitudes rather than the signed terms removes the cancellation and
gives 0.07, which bounds the mean of these modelled shifts rather than
any one of them; the largest single pair among them reaches 0.22.

Extrapolating the 64px inverse-square-root model to 24 trials gives a
mean absolute median component of about 0.16. That is not an estimate of
the full shift in \#299: the issue used 512 reference draws and a
particular synthetic CI, whereas the sweep uses 10,000 draws and models
shifts at the cut-off. No 24-trial scale component is measured here. The
issue’s single reproduction demonstrates the bug, not its typical
impact.

## What it takes to reach a different conclusion

A significance call at `z = 1.96` changes only when the correct InfoVal
lies within the shift of the cut-off *and* the shift points across it.
Applying the measured shifts to a stand-in population of true InfoVals:

``` r
# Each observation is reported as true_z + d0 + true_z * rho, so the scale term
# is evaluated at its own InfoVal rather than at the cut-off.
flip_rate <- function(rows, true_z) {
  reported <- outer(true_z, rows$d0, `+`) + outer(true_z, rows$rho)
  mean((true_z > cutoff) != (reported > cutoff))
}

set.seed(20260906)
spreads <- list("uniform 0-4" = runif(20000, 0, 4),
                "uniform 0-10" = runif(20000, 0, 10),
                "uniform 1.5-2.5" = runif(20000, 1.5, 2.5),
                "uniform 1.94-1.98" = runif(20000, 1.94, 1.98))
flipped <- vapply(spreads, function(z) flip_rate(at_770, z), numeric(1))
knitr::kable(data.frame(true_infoval_spread = names(spreads), flipped = round(flipped, 4)),
             row.names = FALSE)
```

| true_infoval_spread | flipped |
|:--------------------|--------:|
| uniform 0-4         |  0.0139 |
| uniform 0-10        |  0.0054 |
| uniform 1.5-2.5     |  0.0525 |
| uniform 1.94-1.98   |  0.4611 |

These are artificial populations crossed independently with every
sampled reference pair. In these examples the flip rate is 1.4% for the
four-unit uniform spread and 0.5% for the ten-unit spread, rising to 46%
in the near-threshold example. This is not a general monotonic
relationship with distribution width: location, shift direction, and
dependence between InfoVals and reference pairs also matter. None of
these examples estimates an affected-study rate or bounds one beyond the
trivial 0–100%.

And a classification image is only exposed if all of these hold:

1.  InfoVal was computed at all — a feature from 0.4.0 onward, while
    CRAN carried 0.3.4.1 until the 2021 archival.
2.  The stimulus set has two or more base images.
3.  It was generated with the non-default `use_same_parameters = FALSE`.
    In the repository revision examined, this setting appears in no
    vignette, example, test, or InfoVal configuration of the release
    gate — evidence about documented and tested coverage, not about
    researchers’ usage.
4.  The classification image comes from a base image other than the
    first. For k bases, k-1 bases are exposed; the fraction of analyzed
    CIs depends on how many are computed for each base.

## What is not affected

`generateCI()` selects `stimuli_params[[baseimage]]`, so the
classification images and z-maps themselves are built from the right
noise. The bug moves one number, the InfoVal, and only for base images
after the first.

## Limits

How often anyone sets `use_same_parameters = FALSE` is unknown. Issue
\#299 was raised by a code audit and does not provide usage data from
researchers’ studies. Condition 3 above is therefore unquantified rather
than measured to be rare. What this document characterizes is the size
of the damage to an affected classification image, in sampled means and
the largest value among the pairs it drew; it bounds neither that damage
for an individual image nor the number of affected studies.

The references share simulated response draws. Their medians and MADs
retain sampling error, so the measured shift combines differences
between saved noise realizations with finite-reference error under that
coupling. Sharing draws can change the variance of the difference; its
direction and size cannot be inferred from the single-reference control.
The repeated tenfold references diagnose variability at their own
iteration count, not the residual error in the 10,000-draw shift.

Resolution comparisons use 300 trials; trial-count comparisons use 64px
and 128px. The 512px, 770-trial result is an extrapolation with
unresolved resolution dependence and scale trend. The threshold
illustrations use 770 trials at 64px and cannot be transferred to 512px
merely by multiplying by a mean shift ratio. Flip rates depend on the
joint distribution of InfoVals and shifts, not only on mean absolute
shift. The artificial populations assume independence from the sampled
reference pairs; actual study responses need not satisfy that
assumption.

This is not an old-versus-fixed package comparison. A fix must select
the requested base and separate cached references; forwarding
`use_same_parameters` alone is insufficient (#299). Its RNG behavior
depends on the implementation and whether `response_seed` is supplied.
The eventual fix needs its own reproducibility measurement before
quoting a correction in `NEWS.md`.

For orientation, the shared-response mean absolute shift at the cut-off
is 0.05 at 64px and 770 trials, with 0.06 under the 512px, 770-trial
extrapolation. The largest sampled value is 0.35 in the 100-trial
configuration — all of them shifts modelled for an InfoVal sitting at
the cut-off, since this document generates no classification images and
so observes no movement of one. Sample means and maxima bound neither
individual studies nor the population distribution of shifts.
