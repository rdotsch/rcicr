How far InfoVal moves when the target CI’s trial design differs from the
reference’s
================

- [Method: norms through the Gram
  matrix](#method-norms-through-the-gram-matrix)
- [Validation against the package](#validation-against-the-package)
- [Stimulus sets](#stimulus-sets)
- [The designs](#the-designs)
- [Measurements](#measurements)
- [Monte Carlo variability](#monte-carlo-variability)
- [What this shows and what it does
  not](#what-this-shows-and-what-it-does-not)

`computeInfoVal2IFC()` scores a classification image against a reference
distribution simulated from **every** saved stimulus, each given one
random +1/-1 response ([issue
\#349](https://github.com/rdotsch/rcicr/issues/349)). `generateCI()`
also accepts a subset of the saved stimuli, repeated presentations, and
participant IDs. The CI it returns records none of that, so the
reference cannot follow the design. Brinkman et al. (2019, Part I, after
Eq. 1) require the reference to use the identical stimulus set,
including its number of stimuli.

This document measures the size of that mismatch at the package’s
default resolution and trial count. It does not propose or evaluate a
remedy. For each design it compares the package’s reference with a
random-responder null simulated for that design, which we call the
*design-matched null*. For a subset of the saved stimuli, this is the
paper’s requirement applied to the stimuli the target CI was actually
computed from. For participant averages the paper gives no such null.
The design-matched null used here, K independent random responders
averaged the way `generateCI()` averages them, is one candidate
definition. It is not validated by the paper and is not proposed as the
correct reference.

``` r
library(rcicr)

img_size <- 512
nscales <- 5
iter <- 10000          # the package's default reference size
design_iter <- 20000   # draws per design-matched null
cutoff <- 1.96
# Two cores: each forked worker needs its own few hundred MB, and this runs in 3 GB.
cores <- if (.Platform$OS.type == "windows") 1L else min(2L, max(1L, parallel::detectCores() - 1L))
```

## Method: norms through the Gram matrix

A classification image is linear in the noise parameters.
`generateCINoise()` averages the weighted parameter rows and renders the
mean through the basis, so any CI `generateCI()` can build is `S w`.
Here `S` holds one column per saved stimulus’s noise image, and `w` is a
weight per saved stimulus that depends on the design. Its Frobenius norm
is therefore `sqrt(t(w) %*% G %*% w)`, where `G = t(S) %*% S` is the
stimulus Gram matrix. `G` comes from the parameters `X` and the basis
`P`, with `G = X t(P) P t(X)`. `P` is sparse: each pixel belongs to one
patch in each of the basis’s layers.

So a norm costs a stimulus-by-stimulus product instead of a 512 by 512
image, and nulls of 20,000 draws per design become affordable. Every
design reduces to a sparse matrix `A` mapping trials to saved stimuli,
with `w = A %*% responses`:

- **Pooled** (no participant IDs): `generateCI()` first averages
  repeated presentations per stimulus, then averages over the distinct
  stimuli. So `A[s, t] = 1 / (count(s) * n_distinct)`.
- **Participants**: one CI per participant from their own rows, repeats
  included, then the mean of those CIs. So `A[s, t]` gets
  `1 / (K * n_k)` for each of participant `k`’s trials on `s`.

``` r
basis_matrix <- function(p) {
  d <- dim(p$patches)
  npix <- d[1] * d[2]
  Matrix::sparseMatrix(i = rep(seq_len(npix), d[3]), j = as.vector(p$patchIdx),
                       x = as.vector(p$patches) / d[3], dims = c(npix, max(p$patchIdx)))
}

gram <- function(params, basis_gram) params %*% basis_gram %*% t(params)

# The same stream the package consumes: one runif() per trial, column by column.
responses <- function(n, k) matrix(((runif(n * k) > 0.5) * 2) - 1, nrow = n, ncol = k)

pooled_design <- function(stimuli, n_saved) {
  counts <- table(stimuli)
  Matrix::sparseMatrix(i = stimuli, j = seq_along(stimuli),
                       x = 1 / (as.vector(counts[as.character(stimuli)]) * length(counts)),
                       dims = c(n_saved, length(stimuli)))
}

participant_design <- function(stimuli, participants, n_saved) {
  n_k <- table(participants)
  Matrix::sparseMatrix(i = stimuli, j = seq_along(stimuli),
                       x = 1 / (length(n_k) * as.vector(n_k[as.character(participants)])),
                       dims = c(n_saved, length(stimuli)))
}

# Norms of the CIs a random responder produces under design A, in blocks so the
# responses matrix stays near 32 MB however many trials the design has.
null_norms <- function(G, A, draws, seed, block = max(100L, min(2000L, as.integer(4e6 / ncol(A))))) {
  set.seed(seed)
  out <- numeric(draws)
  done <- 0L
  while (done < draws) {
    k <- min(block, draws - done)
    W <- as.matrix(A %*% responses(ncol(A), k))
    out[(done + 1):(done + k)] <- sqrt(colSums(W * (G %*% W)))
    done <- done + k
  }
  out
}

observed_norm <- function(G, A, r) {
  w <- as.vector(A %*% r)
  sqrt(sum(w * (G %*% w)))
}

# One core at 512px: each worker holds its own copy of the basis, and four of
# them do not fit beside it in 3 GB.
generate <- function(n_trials, seed, size, ncores = if (size >= 512) 1L else cores) {
  path <- tempfile("stim")
  dir.create(path)
  base <- file.path(path, "base.png")
  set.seed(seed)
  png::writePNG(matrix(runif(size^2), size, size), base)
  invisible(capture.output(
    generateStimuli2IFC(list(base = base), n_trials = n_trials, img_size = size,
                        stimulus_path = path, seed = seed, nscales = nscales, ncores = ncores,
                        save_as_png = FALSE)
  ))
  saved <- new.env()
  load(list.files(path, pattern = "\\.Rdata$", full.names = TRUE)[1], envir = saved)
  saved$rdata <- list.files(path, pattern = "\\.Rdata$", full.names = TRUE)[1]
  saved
}
```

The base image is synthetic noise. Neither the CI norm nor the reference
depends on it.

## Validation against the package

At 64 pixels, where the package’s own reference is cheap, the
reimplementation is checked end to end. The package’s reference is
generated with a `response_seed`. The Gram null consumes the same stream
and must reproduce it. `generateCI()` then builds three targets from one
saved stimulus set: a subset, repeated presentations, and a participant
average. Each target’s norm must match the Gram computation. Last,
`computeInfoVal2IFC()` scores the subset target against the stored
reference, and the returned InfoVal must equal the one computed from the
Gram norms. That last check confirms the issue’s premise by execution,
not from source alone: a subset CI is scored against the full-set null.

``` r
small <- generate(n_trials = 60, seed = 7, size = 64)
G_small <- gram(small$stimuli_params$base, as.matrix(Matrix::crossprod(basis_matrix(small$p))))

invisible(capture.output(package_reference <- suppressWarnings(generateReferenceDistribution2IFC(
  small$rdata, iter = 500, ncores = 1, response_seed = 3, save_rdata = TRUE))))
gram_reference <- null_norms(G_small, pooled_design(1:60, 60), 500, seed = 3)
max(abs(package_reference - gram_reference))
[1] 1.609823e-15

set.seed(11)
targets <- list(
  subset = list(stimuli = sort(sample(60, 25)), participants = NA),
  repeats = list(stimuli = rep(1:60, 2), participants = NA),
  participants = list(stimuli = rep(1:60, 3), participants = rep(c("a", "b", "c"), each = 60))
)
target_ci <- list()
agreement <- sapply(names(targets), function(name) {
  t <- targets[[name]]
  r <- sample(c(-1, 1), length(t$stimuli), replace = TRUE)
  invisible(capture.output(ci <- generateCI(t$stimuli, r, "base", small$rdata,
                                            participants = t$participants, save_as_png = FALSE,
                                            n_cores = 1)))
  target_ci[[name]] <<- ci
  A <- if (all(is.na(t$participants))) pooled_design(t$stimuli, 60) else
    participant_design(t$stimuli, t$participants, 60)
  abs(norm(matrix(ci$ci), "f") - observed_norm(G_small, A, r))
})
agreement
      subset      repeats participants 
6.661338e-16 1.110223e-16 2.775558e-17 

invisible(capture.output(package_z <- computeInfoVal2IFC(target_ci$subset, small$rdata)))
This classification image was built from 25 of the 60 saved stimuli, but the reference is built over 60. Brinkman et al. (2019) require the reference to use the stimuli the CI was built from. To score it that way, pass reference_stimuli = attr(<your CI>, "trial_design")$stimuli.
gram_z <- (norm(matrix(target_ci$subset$ci), "f") - median(gram_reference)) / mad(gram_reference)
c(package = package_z, gram = gram_z, difference = package_z - gram_z)
     package         gram   difference 
1.158067e+01 1.158067e+01 2.575717e-13 
```

At 512 pixels the package’s reference stores a 262,144-by-770 noise
matrix, which this machine does not hold alongside everything else. At
that size the basis matrix is checked against `generateNoiseImage()` on
the saved stimuli themselves instead. The Gram arithmetic above does not
depend on resolution.

## Stimulus sets

Two saved sets at the package defaults, 770 trials and 512 pixels. Two
more at 300 trials check that the picture does not hinge on the trial
count.

``` r
basis <- basis_matrix(generateNoisePattern(img_size, nscales = nscales))
basis_gram <- as.matrix(Matrix::crossprod(basis))

sets <- lapply(list(`770a` = c(770, 21), `770b` = c(770, 22), `300a` = c(300, 23),
                    `300b` = c(300, 24)), function(cfg) {
  saved <- generate(cfg[1], cfg[2], img_size)
  params <- saved$stimuli_params$base
  check <- max(sapply(1:3, function(i)
    max(abs(as.vector(basis %*% params[i, ]) - as.vector(generateNoiseImage(params[i, ], saved$p))))))
  rm(saved)
  invisible(gc())
  list(n = cfg[1], G = gram(params, basis_gram), basis_check = check)
})
# Forked workers below copy whatever the parent still holds as soon as they
# collect garbage, so the parent drops everything but the Gram matrices.
rm(basis, basis_gram, small, G_small)
invisible(gc())
sapply(sets, `[[`, "basis_check")
        770a         770b         300a         300b 
1.942890e-16 1.387779e-16 1.665335e-16 1.942890e-16 
```

## The designs

Every design below uses the saved stimuli of one set, with `N` saved
stimuli.

- **Missing trials**: a random subset of `N`, as when timed-out or
  excluded trials are dropped before `generateCI()`. From half a percent
  to a fifth of the trials missing.
- **Partial trial sets**: half, a quarter, and 100 of 770.
- **Participant averages**: `K` participants who each complete all `N`
  trials, averaged by passing `participants`. A pooled call without IDs
  over the same `K` presentations per stimulus gives the identical
  weights, since per-stimulus averaging and per-participant averaging
  agree when every participant sees every stimulus once.
- **Participant averages with missing trials**: 20 participants, each
  missing a different random 5%.
- **Disjoint equal blocks**: 5 participants, each completing a different
  fifth of the set. Every stimulus gets weight `w = r / N` exactly as in
  the reference, so this is a control. Its design-matched null has the
  same distribution as the reference.

``` r
designs <- function(N, seed) {
  set.seed(seed)
  subset <- function(keep) {
    s <- sort(sample(N, keep))
    pooled_design(s, N)
  }
  everyone <- function(K) participant_design(rep(seq_len(N), K), rep(seq_len(K), each = N), N)
  missing <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20)
  partial <- c(0.5, 0.25)
  out <- c(
    setNames(lapply(missing, function(f) subset(N - round(f * N))),
             sprintf("missing %s%%", 100 * missing)),
    setNames(lapply(partial, function(f) subset(round(f * N))), sprintf("%s%% of trials", 100 * partial)),
    if (N == 770) list(`100 of 770` = subset(100)),
    setNames(lapply(c(2, 5, 10, 30), everyone), sprintf("mean of %d participants", c(2, 5, 10, 30))),
    list(`mean of 20, 5% missing each` = {
      trials <- lapply(1:20, function(k) sort(sample(N, N - round(0.05 * N))))
      participant_design(unlist(trials), rep(1:20, lengths(trials)), N)
    }),
    list(`5 disjoint blocks (control)` = participant_design(seq_len(N), rep(1:5, each = N / 5), N))
  )
  out
}
```

## Measurements

Each set gets the package’s reference, 10,000 draws from all `N` saved
stimuli with one response each. Each design gets its design-matched null
of 20,000 draws. The comparison reports:

- `z_null`: the InfoVal the package reports for the median
  random-responder CI of that design. Under a matched reference it is 0
  by construction.
- `scale`: the design-matched null’s MAD over the reference’s.
- `z_at_1.96`: the InfoVal the package reports for a CI that sits at
  1.96 on the design-matched null. InfoVal is linear in the norm, so
  this is `z_null + 1.96 * scale`.
- `false_pos`: the share of the design-matched null that the package
  scores above 1.96. `matched_pos` is the same share scored against the
  design-matched null itself. InfoVal is not normally distributed, so
  the nominal 2.5% is not the right comparison for `false_pos`.

``` r
measure_set <- function(name, set, design_seed) {
  reference <- null_norms(set$G, pooled_design(seq_len(set$n), set$n), iter, seed = 1)
  ref_med <- median(reference)
  ref_mad <- mad(reference)
  ds <- designs(set$n, design_seed)
  rows <- parallel::mclapply(seq_along(ds), mc.cores = cores, function(i) {
    nulls <- null_norms(set$G, ds[[i]], design_iter, seed = 1000 + i)
    med <- median(nulls)
    s <- mad(nulls)
    data.frame(set = name, design = names(ds)[i],
               z_null = (med - ref_med) / ref_mad, scale = s / ref_mad,
               z_at_1.96 = (med + cutoff * s - ref_med) / ref_mad,
               false_pos = mean((nulls - ref_med) / ref_mad > cutoff),
               matched_pos = mean((nulls - med) / s > cutoff))
  })
  do.call(rbind, rows)
}

results <- do.call(rbind, Map(measure_set, names(sets), sets, seq_along(sets)))
```

Set 770a, the first 770-trial set:

``` r
shown <- results
shown[, 3:5] <- lapply(shown[, 3:5], round, 2)
shown[, 6:7] <- lapply(shown[, 6:7], function(x) sprintf("%.1f%%", 100 * x))
knitr::kable(subset(shown, set == "770a", -set), row.names = FALSE)
```

| design                      | z_null | scale | z_at_1.96 | false_pos | matched_pos |
|:----------------------------|-------:|------:|----------:|:----------|:------------|
| missing 0.5%                |   0.05 |  0.99 |      1.98 | 4.5%      | 4.4%        |
| missing 1%                  |   0.11 |  1.00 |      2.06 | 5.2%      | 4.4%        |
| missing 2%                  |   0.21 |  1.01 |      2.19 | 6.2%      | 4.3%        |
| missing 5%                  |   0.53 |  1.02 |      2.53 | 10.3%     | 4.4%        |
| missing 10%                 |   1.14 |  1.04 |      3.18 | 23.3%     | 4.5%        |
| missing 20%                 |   2.50 |  1.09 |      4.65 | 69.4%     | 4.6%        |
| 50% of trials               |   8.76 |  1.41 |     11.51 | 100.0%    | 4.0%        |
| 25% of trials               |  21.20 |  1.98 |     25.09 | 100.0%    | 4.1%        |
| 100 of 770                  |  37.35 |  2.73 |     42.69 | 100.0%    | 4.0%        |
| mean of 2 participants      |  -6.19 |  0.76 |     -4.71 | 0.0%      | 4.2%        |
| mean of 5 participants      | -11.69 |  0.49 |    -10.72 | 0.0%      | 4.0%        |
| mean of 10 participants     | -14.46 |  0.35 |    -13.77 | 0.0%      | 4.0%        |
| mean of 30 participants     | -17.29 |  0.21 |    -16.88 | 0.0%      | 3.7%        |
| mean of 20, 5% missing each | -16.30 |  0.26 |    -15.80 | 0.0%      | 3.7%        |
| 5 disjoint blocks (control) |   0.01 |  0.98 |      1.94 | 4.0%      | 4.1%        |

The same designs on the other three sets, `z_null` only:

``` r
wide <- reshape(results[, c("set", "design", "z_null")], idvar = "design", timevar = "set",
                direction = "wide")
names(wide) <- sub("z_null.", "", names(wide), fixed = TRUE)
wide[, -1] <- lapply(wide[, -1], round, 2)
knitr::kable(wide, row.names = FALSE)
```

| design                      |   770a |   770b |   300a |   300b |
|:----------------------------|-------:|-------:|-------:|-------:|
| missing 0.5%                |   0.05 |   0.06 |   0.06 |   0.07 |
| missing 1%                  |   0.11 |   0.09 |   0.09 |   0.11 |
| missing 2%                  |   0.21 |   0.21 |   0.22 |   0.23 |
| missing 5%                  |   0.53 |   0.51 |   0.57 |   0.58 |
| missing 10%                 |   1.14 |   1.12 |   1.14 |   1.21 |
| missing 20%                 |   2.50 |   2.45 |   2.48 |   2.52 |
| 50% of trials               |   8.76 |   8.79 |   8.81 |   8.88 |
| 25% of trials               |  21.20 |  21.04 |  20.89 |  21.18 |
| 100 of 770                  |  37.35 |  37.64 |     NA |     NA |
| mean of 2 participants      |  -6.19 |  -6.19 |  -6.19 |  -6.23 |
| mean of 5 participants      | -11.69 | -11.65 | -11.67 | -11.76 |
| mean of 10 participants     | -14.46 | -14.41 | -14.43 | -14.55 |
| mean of 30 participants     | -17.29 | -17.23 | -17.25 | -17.40 |
| mean of 20, 5% missing each | -16.30 | -16.24 | -16.26 | -16.40 |
| 5 disjoint blocks (control) |   0.01 |  -0.01 |  -0.01 |   0.02 |

## Monte Carlo variability

The yardstick is the spread of InfoVal across independent 10,000-draw
references for the same stimulus set and the same CI. It is measured at
the reference median and at 1.96 on the population reference, with a
200,000-draw reference standing in for the population. It is compared
with `z_null` for the smallest departures measured.

``` r
set <- sets[["770a"]]
full <- pooled_design(seq_len(set$n), set$n)
population <- null_norms(set$G, full, 200000, seed = 2)
at <- median(population) + c(0, cutoff) * mad(population)
replicate_z <- do.call(rbind, parallel::mclapply(1:40, mc.cores = cores, function(i) {
  ref <- null_norms(set$G, full, iter, seed = 5000 + i)
  (at - median(ref)) / mad(ref)
}))
mc_sd <- apply(replicate_z, 2, sd)
names(mc_sd) <- c("at median", "at 1.96")
round(mc_sd, 3)
at median   at 1.96 
    0.014     0.032 

small_missing <- subset(results, set == "770a" & grepl("^missing", design))
round(setNames(small_missing$z_null / mc_sd[["at median"]], small_missing$design), 1)
missing 0.5%   missing 1%   missing 2%   missing 5%  missing 10%  missing 20% 
         3.3          7.4         14.8         37.4         80.3        175.7 
```

The 40 replicates estimate a standard deviation, and that estimate
carries its own sampling error of about 11% of its value. The ratios
above are not sharper than that.

## What this shows and what it does not

In set 770a (770 trials, 512 pixels), the median random-responder CI
built from all but a random half percent of the saved stimuli gets an
InfoVal of 0.05, not 0. That is 3.3 times the Monte Carlo standard
deviation of a 10,000-draw reference. With 5% missing it is 0.53, with
20% missing 2.5, and a CI from 100 of the 770 stimuli gets 37.3. Across
all four sets, the median ranges from 0.09 to 0.11 with 1% missing and
from 0.51 to 0.58 with 5%. The direction is upward in every subset
design and every set, because fewer stimuli average away less noise. In
set 770a, against the package’s reference, 10.3% of pure-noise CIs with
5% of trials missing score above 1.96, against 4.4% under their own
null.

Participant averages move the other way and further. In set 770a, two
participants who each saw every stimulus put their median
random-responder average at -6.2. A group CI that would sit at 1.96 on
the two-responder null is reported as -4.7. With thirty participants
that becomes -16.9. Whether the design-matched null is the right
reference for a group CI is the methodological question this document
leaves open. What it establishes is that the full-set, single-responder
reference is not calibrated for this null either.

The disjoint-block control reads 0.01, within the Monte Carlo spread of
zero, as its construction requires. The 300-trial sets and the second
770-trial set show the same pattern.

This document does not measure:

- how often real studies drop trials, subset stimuli, or average
  participants before scoring, so no prevalence of affected InfoVals;
- how much any published conclusion changes;
- which reference is correct for a participant average, or whether one
  reference can serve a subset and a group design at once;
- what a remedy should be, or whether an established full-set result
  would move under one. A full-set, single-response design has the same
  weights as the reference by construction.
