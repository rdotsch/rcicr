#!/usr/bin/env Rscript
#
# What the individual-CI mislabelling did to a second-stage analysis.
#
# The bug (fixed in 1.3.0, see notes/individual-ci-mislabelling.md) permuted
# output filenames: the file named unique(participants)[i] held the
# classification image of sort(unique(participants))[i]. An analysis that
# joins a rating of each image back to the participant its filename names is
# therefore analysing a permuted pairing.
#
# This script measures the consequence, because the advisory makes a claim
# about it that is easy to get wrong in either direction -- "it only adds
# noise" and "it could have produced spurious findings" are both wrong as
# stated. Four designs, at the alpha = 0.05 level:
#
#   A  covariate unrelated to the order participants were labelled in
#      (the ordinary individual-differences design)
#   B  condition assigned in blocks by collection order
#   C  covariate that is collection order itself -- testing date, cohort
#   D  a real trend over collection order with no true condition effect,
#      condition blocked by collection order -- the case for asking whether
#      the permutation can produce a hypothesis-consistent result rather
#      than only destroy one
#
# Usage: Rscript tools/simulate-mislabelling-impact.R [nsim]
# Cost is nsim x cells x two tests per iteration; pass a smaller nsim to
# iterate, and expect the wall clock to scale from there.

nsim_arg <- commandArgs(trailingOnly = TRUE)[1]
NSIM <- if (is.na(nsim_arg)) 20000L else as.integer(nsim_arg)
ALPHA <- 0.05
set.seed(20260819)

# The permutation the bug applied, as an index into collection order: reading
# the file named a[i] hands you the image of s[i].
mislabel_perm <- function(ids) {
  a <- unique(ids)
  s <- sort(unique(ids))
  match(s, a)
}

cat("nsim =", NSIM, " alpha =", ALPHA, "\n\n")

cat("=== How much of the pairing survives ===\n")
for (n in c(9, 12, 20, 30, 50, 100)) {
  p <- mislabel_perm(paste0("p", seq_len(n)))
  cat(sprintf("  p1..p%-3d  %3d of %3d files keep their own participant (%5.1f%%)\n",
              n, sum(p == seq_len(n)), n, 100 * mean(p == seq_len(n))))
}
cat(sprintf("  p01..p12  %3d of  12 (zero-padded IDs sort as collected)\n\n",
            sum(mislabel_perm(sprintf("p%02d", 1:12)) == 1:12)))

# Correlational designs (A and C differ only in how the covariate is drawn).
run_cor <- function(n, rho, xfun, nsim = NSIM) {
  perm <- mislabel_perm(paste0("p", seq_len(n)))
  hit_c <- hit_a <- opp <- logical(nsim)
  r_a <- numeric(nsim)
  for (i in seq_len(nsim)) {
    x <- xfun(n)
    y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    tc <- cor.test(x, y)
    ta <- cor.test(x, y[perm])
    hit_c[i] <- tc$p.value < ALPHA
    hit_a[i] <- ta$p.value < ALPHA
    opp[i] <- ta$p.value < ALPHA && ta$estimate < 0
    r_a[i] <- ta$estimate
  }
  # Both analyses see the same y within an iteration, so the paired difference
  # is the statistic that carries the claim -- its error is far smaller than
  # either rate's, and it is not confounded with the test's own size at small
  # n (Welch is conservative at n = 12, in the correct column too).
  discordant <- sum(hit_c != hit_a)
  data.frame(n = n, rho = rho,
             detected_correct = mean(hit_c), detected_affected = mean(hit_a),
             diff = mean(hit_a) - mean(hit_c),
             diff_se = sqrt(discordant) / nsim,
             mean_r_affected = mean(r_a), sig_opposite = mean(opp))
}

grid <- expand.grid(n = c(12, 20, 30, 50), rho = c(0, 0.3, 0.5))

cat("=== A: covariate unrelated to labelling order ===\n")
resA <- do.call(rbind, Map(function(n, r) run_cor(n, r, function(k) rnorm(k)),
                           grid$n, grid$rho))
print(resA, row.names = FALSE, digits = 3)

cat("\n=== C: covariate is collection order itself (testing date, cohort) ===\n")
resC <- do.call(rbind, Map(function(n, r) run_cor(n, r, function(k) as.numeric(scale(seq_len(k)))),
                           grid$n, grid$rho))
print(resC, row.names = FALSE, digits = 3)


cat("\n=== B: condition assigned in blocks by collection order ===\n")
run_blocked <- function(n, d, nsim = NSIM) {
  perm <- mislabel_perm(paste0("p", seq_len(n)))
  cond <- rep(c(0, 1), each = n / 2)
  hit_c <- hit_a <- sig_p <- sig_o <- logical(nsim)
  for (i in seq_len(nsim)) {
    y <- d * cond + rnorm(n)
    tc <- t.test(y[cond == 1], y[cond == 0])
    ya <- y[perm]
    ta <- t.test(ya[cond == 1], ya[cond == 0])
    est <- ta$estimate[1] - ta$estimate[2]
    hit_c[i] <- tc$p.value < ALPHA
    hit_a[i] <- ta$p.value < ALPHA
    sig_p[i] <- ta$p.value < ALPHA && est > 0
    sig_o[i] <- ta$p.value < ALPHA && est < 0
  }
  discordant <- sum(hit_c != hit_a)
  data.frame(n = n, d = d, detected_correct = mean(hit_c),
             detected_affected = mean(hit_a),
             diff = mean(hit_a) - mean(hit_c), diff_se = sqrt(discordant) / nsim,
             sig_predicted = mean(sig_p),
             sig_opposite = mean(sig_o), crossed_condition = mean(cond[perm] != cond))
}
gridB <- expand.grid(n = c(12, 20, 30, 50), d = c(0, 0.8))
resB <- do.call(rbind, Map(function(n, d) run_blocked(n, d), gridB$n, gridB$d))
print(resB, row.names = FALSE, digits = 3)

# The null cells are the evidence for "cannot create an association", so they
# get their own table. Correct and mislabelled analyses see the same outcome
# vector within an iteration, so the paired difference is the statistic: its
# error is much smaller than either rate's, and it does not confound the
# permutation with the test's own size. That distinction matters here --
# Welch is mildly conservative at n = 12, visible in B's *correct* column.
cat("\n=== Null cells, paired: is the mislabelled rate the same as the correct one? ===\n")
nullcells <- rbind(
  data.frame(design = "A", test = "cor.test", resA[resA$rho == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se")]),
  data.frame(design = "C", test = "cor.test", resC[resC$rho == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se")]),
  data.frame(design = "B", test = "Welch t", resB[resB$d == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se")])
)
nullcells$SEs_from_zero <- abs(nullcells$diff) / nullcells$diff_se
print(nullcells, row.names = FALSE, digits = 3)

# D asks whether the permutation can hand back a hypothesis-consistent result
# rather than only destroy one. Scenario B at d = 0 already answers it for an
# exchangeable outcome, where nothing can be created; the harder case is an
# outcome carrying real structure tied to collection order -- drift, practice,
# season -- with no true condition effect, where the permutation has something
# to move around.
cat("\n=== D: real trend over collection order, no true condition effect ===\n")
run_trend <- function(n, trend, nsim = NSIM) {
  perm <- mislabel_perm(paste0("p", seq_len(n)))
  cond <- rep(c(0, 1), each = n / 2)
  ord <- as.numeric(scale(seq_len(n)))
  cp <- co <- ap <- ao <- logical(nsim)
  for (i in seq_len(nsim)) {
    y <- trend * ord + rnorm(n)
    tc <- t.test(y[cond == 1], y[cond == 0])
    ya <- y[perm]
    ta <- t.test(ya[cond == 1], ya[cond == 0])
    ec <- tc$estimate[1] - tc$estimate[2]
    ea <- ta$estimate[1] - ta$estimate[2]
    cp[i] <- tc$p.value < ALPHA && ec > 0
    co[i] <- tc$p.value < ALPHA && ec < 0
    ap[i] <- ta$p.value < ALPHA && ea > 0
    ao[i] <- ta$p.value < ALPHA && ea < 0
  }
  data.frame(n = n, trend = trend,
             correct_predicted = mean(cp), correct_opposite = mean(co),
             mislabelled_predicted = mean(ap), mislabelled_opposite = mean(ao))
}
gridD <- expand.grid(n = c(12, 20, 30, 50), trend = c(0.3, 0.6, 1.0))
print(do.call(rbind, Map(function(n, t) run_trend(n, t), gridD$n, gridD$trend)),
      row.names = FALSE, digits = 3)

cat("\n=== Summary ===\n")
cat("1. No true association -> the mislabelling cannot create one. Rejection\n")
cat("   stays at alpha in A, B and C, at every n. Permuting one side of a pair\n")
cat("   that is independent of the other leaves it independent, so this is\n")
cat("   structural rather than a property of these particular draws.\n")
cat("2. A true association, ordinary design -> it is destroyed. Detection\n")
cat("   collapses to alpha; the runs that fail to reject -- the other ~95% at\n")
cat("   rho = 0.5, N = 50 -- are the false negatives.\n")
cat("3. A covariate tracking labelling order -> the residual can be either\n")
cat("   sign: an effect can survive attenuated, or come back reversed.\n")
cat("4. Even in D, where the permutation has real structure to move around, it\n")
cat("   never raises hypothesis-consistent rejections above what the correctly\n")
cat("   labelled analysis gives -- it lowers them, and diverts some into the\n")
cat("   opposite direction. It has no systematic pull toward a prediction. A\n")
cat("   single affected study can still land on a significant result in the\n")
cat("   predicted direction, which is why re-running is required either way.\n")
