What the wrong-base InfoVal reference does to a conclusion
================

- [The reference distribution ignores which base you
  score](#the-reference-distribution-ignores-which-base-you-score)
- [Base 2’s own reference](#base-2s-own-reference)
- [The error is scatter, with a directional part too small to
  matter](#the-error-is-scatter-with-a-directional-part-too-small-to-matter)
- [How large the scatter is](#how-large-the-scatter-is)
- [Trials set the size; pixels are not seen
  to](#trials-set-the-size-pixels-are-not-seen-to)
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

This document measures what that costs. It exists because the two
obvious readings are both wrong: “the reference is from the wrong base,
so the InfoVal is invalid” overstates it, and “it is only noise, so it
does not matter” understates it. The answer is that the wrong reference
is a *valid* null built on a different realization of the same noise
process, so the error is scatter — with a directional part some two
orders of magnitude below it — whose size is set by the number of
trials, and a significance call flips only in a narrow window around the
cut-off.

Knitting this runs the whole measurement, which takes tens of minutes:
82 base image pairs, each needing two stimulus sets and three reference
distributions of 10,000 iterations.

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
generate <- function(bases, n_trials, seed, size, same_parameters) {
  path <- tempfile("stim", tmpdir = scratch)
  dir.create(path)
  invisible(capture.output(
    generateStimuli2IFC(bases, n_trials = n_trials, img_size = size, stimulus_path = path,
                        seed = seed, use_same_parameters = same_parameters, nscales = nscales,
                        ncores = 1, save_as_png = FALSE)
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
  vapply(seq_len(nrow(params)),
         function(trial) as.vector(generateNoiseImage(params[trial, ], p)),
         numeric(size^2))
}

reference_norms <- function(noise, iter, seed, block = 500L) {
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

That reimplementation is only trustworthy if it reproduces the package’s
own output, so every configuration below checks it against
`generateReferenceDistribution2IFC()` on the first base before using it
on the second. The same check is what establishes that the *first* base
is scored correctly: the package’s reference matches a null rebuilt from
base 1’s saved parameters.

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
    mc_rho = mad(monte_carlo) / mad(scored_against) - 1
  )
}

sweep <- function(n_trials, seeds, size = img_size) {
  do.call(rbind, Map(function(seed, validate) measure_pair(n_trials, seed, size, validate),
                     seeds, seq_along(seeds) == 1))
}
```

`d0` is the shift in the reference median, expressed in the units
InfoVal is reported in, and `rho` is the error in its scale. An InfoVal
that should have been `z` is reported as `z + d0 + z * rho`, so `d0` is
what a threshold call sees and `rho` matters only for large values.

``` r
results <- rbind(
  sweep(770, 1:30),
  sweep(300, 101:120),
  sweep(100, 201:220),
  sweep(300, 301:312, size = 128)
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

## The error is scatter, with a directional part too small to matter

``` r
at_default <- subset(results, n_trials == 770)

knitr::kable(round(rbind(
  "reference median only (d0)" = centre(at_default$d0),
  "scale only (rho)" = centre(at_default$rho),
  "full shift at the cut-off" = centre(shift(at_default, cutoff)),
  "full shift at z = 10" = centre(shift(at_default, 10))
), 4))
```

|                            |    mean |     se | ses_from_zero |
|:---------------------------|--------:|-------:|--------------:|
| reference median only (d0) | -0.0066 | 0.0065 |        1.0074 |
| scale only (rho)           | -0.0020 | 0.0038 |        0.5409 |
| full shift at the cut-off  | -0.0106 | 0.0114 |        0.9282 |
| full shift at z = 10       | -0.0269 | 0.0403 |        0.6684 |

Base 2’s parameters are an independent draw from the same generator as
base 1’s, so the null it receives is a legitimate one — just built on
another realization of the noise. Across 30 base image pairs at the
package defaults the full shift at the cut-off averages -0.011, which is
0.9 standard errors from zero against a scatter 6 times the size of that
mean.

Neither term is exactly centred, though, and the sign of each is decided
by the same structure. `rho` is the ratio
`mad(own) / mad(scored_against)` less one, and for exchangeable positive
scales `B/A + A/B >= 2` forces `E[B/A] >= 1`, strictly whenever the two
differ at all. `d0` divides by the MAD of one member of the pair only,
so its mean is `-Cov(median, 1 / mad)` rather than zero: a noise
realization with more energy has both a larger median norm and a wider
spread. Both terms therefore carry a directional component, and the
question is its size rather than its existence.

Neither mean can measure it — both sit far below their own standard
errors. Averaging a pair’s two labellings can: had the second base image
been listed first, the bug would have scored the first against *its*
null, and the two orderings share the bulk of the variation, which
cancels.

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

The estimates come out positive wherever they resolve against their own
standard error, which is what the argument above requires of the
expectations, and all of them are minute. At the defaults the
directional part of the shift is 5.3e-04 at the cut-off and 0.002 at an
InfoVal of 10 — around 1% of the scatter at that value, a proportion
that does not grow with the InfoVal because both scale with it. A set of
published InfoVals was therefore scattered rather than moved, to a
precision far finer than anyone reads these numbers to.

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

`monte_carlo` is the same quantity, computed the same way down to the
scale term, for a reference built from the same base image with a
different draw of the simulated responses: the wobble any InfoVal
already carries at `iter` = 10,000. At the default 770 trials the bug is
worth 1.6 times that wobble. It is the same kind of error as the one the
method already tolerates, and not far from the same size.

The scale error `rho` barely exceeds its own Monte Carlo control (0.0149
against 0.0134), so almost all of the damage is the shift in the median,
and the `z * rho` term stays small until the InfoVal itself is large.

``` r
at_z <- function(z, f) vapply(z, function(zz) mean(abs(f(at_default, zz))), numeric(1))

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

The margin over the Monte Carlo wobble narrows as the InfoVal grows,
because both quantities pick up the same scale term. At a large InfoVal
the bug is barely distinguishable in size from re-running the reference
simulation with different random responses.

## Trials set the size; pixels are not seen to

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

The shift falls as one over the square root of the number of trials,
which is what a reference median built from an average over trials
should do.

``` r
compare_sizes <- function(f) {
  sides <- lapply(c(64, 128), function(s) abs(f(subset(results, n_trials == 300 & img_size == s))))
  gap <- diff(vapply(sides, mean, numeric(1)))
  se <- sqrt(sum(vapply(sides, function(x) var(x) / length(x), numeric(1))))
  c(`64px` = mean(sides[[1]]), `128px` = mean(sides[[2]]), difference = gap, se = se,
    ses_from_zero = abs(gap) / se,
    largest_rise_consistent = (gap + 1.96 * se) / mean(sides[[1]]))
}

knitr::kable(signif(rbind(
  d0 = compare_sizes(function(rows) rows$d0),
  `shift at the cut-off` = compare_sizes(function(rows) shift(rows, cutoff))
), 3))
```

|                      |   64px |  128px | difference |     se | ses_from_zero | largest_rise_consistent |
|:---------------------|-------:|-------:|-----------:|-------:|--------------:|------------------------:|
| d0                   | 0.0449 | 0.0512 |     0.0063 | 0.0107 |          0.59 |                   0.607 |
| shift at the cut-off | 0.0577 | 0.0729 |     0.0152 | 0.0143 |          1.06 |                   0.749 |

Quadrupling the pixels moves neither quantity detectably. That is a weak
constraint rather than a demonstration of independence: with 20 and 12
pairs, the comparison is consistent with the shift at the cut-off being
up to 75% larger at the higher resolution, and the default 512px is a
further factor of four away and not measured here at all.

What makes the extrapolation reasonable rather than established is that
the quantity is dimensionless: adding pixels raises the reference median
and its spread together, and InfoVal reports their ratio. Treat the
trial-count rows as the finding and 512px as an inference from them.

This also puts the 24-trial reproduction in the issue in proportion.
Extrapolating the fitted scaling to 24 trials gives a typical shift of
0.16, so the −0.067 observed there is one unremarkable draw, and a real
study with hundreds of trials sits several times below it.

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
flipped <- vapply(spreads, function(z) flip_rate(at_default, z), numeric(1))
knitr::kable(data.frame(true_infoval_spread = names(spreads), flipped = round(flipped, 4)),
             row.names = FALSE)
```

| true_infoval_spread | flipped |
|:--------------------|--------:|
| uniform 0-4         |  0.0139 |
| uniform 0-10        |  0.0054 |
| uniform 1.5-2.5     |  0.0525 |
| uniform 1.94-1.98   |  0.4611 |

Those spreads are stand-ins, not measurements: what InfoVals look like
across real studies is not something this package records, and none of
them bounds the answer. The rate is set by how much of a study’s InfoVal
distribution sits within the shift of the cut-off, so it falls as the
distribution broadens — 1.4% over a four-unit spread, 0.5% over ten —
and rises without a useful bound as it narrows. The last row makes that
concrete: concentrate every InfoVal within a few hundredths of 1.96 and
46% of the calls change, approaching the probability that the shift
points across the threshold at all.

Which row a given study resembles is the whole question, and it is not
one this package can answer: nothing records what InfoVals look like in
practice. A study whose classification images mostly carry strong signal
sits near the top rows; one that reports values hovering at the
threshold sits near the bottom one, and for it the rate is high. The
effect of the bug is bounded only by where a study’s InfoVals actually
fall.

And a classification image is only exposed if all of these hold:

1.  InfoVal was computed at all — a feature from 0.4.0 onward, while
    CRAN carried 0.3.4.1 until the 2021 archival.
2.  The stimulus set has two or more base images.
3.  It was generated with `use_same_parameters = FALSE`, which is not
    the default and appears in no vignette, example, test, or InfoVal
    configuration of the release gate.
4.  The classification image comes from a base image other than the
    first — `(k-1)/k` of the images in such a study.

## What is not affected

`generateCI()` selects `stimuli_params[[baseimage]]`, so the
classification images and z-maps themselves are built from the right
noise. The bug moves one number, the InfoVal, and only for base images
after the first.

## Limits

How often anyone sets `use_same_parameters = FALSE` is unknown. There is
no telemetry, and the only issue naming the setting is \#299 itself,
raised by an audit of the code rather than by anyone reporting it from
their own work — which is an absence of reports, not evidence of an
absence of use. Condition 3 above is therefore unquantified rather than
measured to be rare, and this document bounds the damage per affected
classification image, not the number of affected studies.

The two references compared here share their simulated response draws,
which isolates the noise realization from Monte Carlo error. A fix that
draws base 2’s responses from its own position in the stream adds
roughly the `monte_carlo` column on top, which changes none of the
conclusions above.

The shift measured here is the correction a fix applies, so it is also
the size to quote in a `NEWS.md` “Reproducibility impact” entry: about
0.05 in z at the defaults, scaling as one over the square root of the
trial count, scatter rather than a correction in a direction, and
affecting no base image but the second and later.
