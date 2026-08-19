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
# about it and the two summaries it invites -- "it only adds noise" and "it
# could have produced spurious findings" -- are each half right: the aggregate
# false positive rate is untouched while individual verdicts flip both ways,
# and in some designs the shuffle raises apparent support rather than only
# destroying it. Five designs, at the alpha = 0.05 level:
#
#   A  covariate unrelated to the order participants were labelled in
#      (the ordinary individual-differences design)
#   B  condition assigned in blocks by collection order
#   C  covariate that is collection order itself -- testing date, cohort
#   D  a real trend over collection order with no true condition effect,
#      condition blocked by collection order -- the case for asking whether
#      the permutation can produce a hypothesis-consistent result rather
#      than only destroy one
#   E  how the damage scales with how far the identifiers departed from
#      sorted order -- A to D fix that at close to the worst case
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
             decisions_changed = discordant / nsim,
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
             decisions_changed = discordant / nsim, sig_predicted = mean(sig_p),
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
  data.frame(design = "A", test = "cor.test", resA[resA$rho == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se", "decisions_changed")]),
  data.frame(design = "C", test = "cor.test", resC[resC$rho == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se", "decisions_changed")]),
  data.frame(design = "B", test = "Welch t", resB[resB$d == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se", "decisions_changed")])
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
# Assignment scheme matters, and picking only one of them rigs the answer.
# Contiguous blocks put every early observation in one condition and every
# late one in the other, so the CORRECT pairing already maximises the trend
# contrast and any permutation can only reduce it -- "mislabelled never
# exceeds correct" would then be true by construction rather than measured.
# Alternating assignment is the opposite case: the correct pairing spreads
# the trend evenly across conditions, leaving almost no contrast for the
# permutation to destroy and plenty of room to create one. At n = 12 the
# expected contrast in position units goes 6.00 -> 0.00 blocked, and
# 1.00 -> 2.33 alternating.
run_trend <- function(n, trend, scheme, nsim = NSIM) {
  perm <- mislabel_perm(paste0("p", seq_len(n)))
  cond <- if (scheme == "blocked") rep(c(0, 1), each = n / 2) else rep(c(0, 1), times = n / 2)
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
  data.frame(n = n, trend = trend, assignment = scheme,
             correct_predicted = mean(cp), correct_opposite = mean(co),
             mislabelled_predicted = mean(ap), mislabelled_opposite = mean(ao),
             amplified = mean(ap) > mean(cp))
}
gridD <- expand.grid(n = c(12, 20, 30, 50), trend = c(0.3, 0.6, 1.0),
                     scheme = c("blocked", "alternating"), stringsAsFactors = FALSE)
resD <- do.call(rbind, Map(function(n, t, s) run_trend(n, t, s),
                           gridD$n, gridD$trend, gridD$scheme))
print(resD, row.names = FALSE, digits = 3)


# E: how much damage depends on how far out of order the identifiers were.
# A and C hard-code "p1 ... pN" in collection order, which past ten
# participants is close to the worst case -- exactly one file keeps its own
# participant. An affected study can be far milder: zero-padded identifiers
# entered with a couple of participants out of sequence swap only those, and
# attenuate an association rather than destroying it. E sweeps that severity
# using the real mechanism -- build an appearance order, let mislabel_perm()
# derive the permutation from it -- rather than by constructing permutations
# directly.
cat("\n=== E: severity as a function of how out-of-order the IDs were ===\n")
run_severity <- function(n, rho, swaps, nsim = NSIM) {
  ids <- sprintf("p%02d", seq_len(n))          # zero-padded: sorts as collected
  appearance <- ids
  if (swaps > 0) {                             # displace `swaps` adjacent pairs
    at <- seq(1, by = 2, length.out = swaps)
    appearance[c(rbind(at, at + 1))] <- appearance[c(rbind(at + 1, at))]
  }
  perm <- mislabel_perm(appearance)
  kept <- mean(perm == seq_len(n))
  hit_c <- hit_a <- logical(nsim)
  for (i in seq_len(nsim)) {
    x <- rnorm(n)
    y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
    hit_c[i] <- cor.test(x, y)$p.value < ALPHA
    hit_a[i] <- cor.test(x, y[perm])$p.value < ALPHA
  }
  data.frame(n = n, rho = rho, pairs_swapped = swaps,
             share_correctly_paired = kept,
             detected_correct = mean(hit_c), detected_affected = mean(hit_a))
}
gridE <- expand.grid(n = 50, rho = 0.5, swaps = c(1, 2, 5, 10, 25))
resE <- do.call(rbind, Map(function(n, r, s) run_severity(n, r, s),
                           gridE$n, gridE$rho, gridE$swaps))
print(resE, row.names = FALSE, digits = 3)
cat("  for comparison, the p1..pN scheme used above keeps",
    sprintf("%.2f", mean(mislabel_perm(paste0("p", 1:50)) == 1:50)),
    "of 50 correctly paired\n")

# Every figure below is computed from this run, so it stays true when nsim
# is overridden rather than quoting the seeded 20,000-iteration defaults.
nullrows <- nullcells
dc <- range(nullrows$decisions_changed)
alt <- resD[resD$assignment == "alternating", ]
amp <- alt[alt$amplified, ]
worst <- alt[which.max(alt$mislabelled_predicted - alt$correct_predicted), ]
ord <- resA[resA$rho == 0.5 & resA$n == 50, ]

cat("\n=== Summary ===\n")
cat(sprintf(paste0(
  "1. No true association -> the false positive RATE is unchanged. Rejection\n",
  "   stays at alpha in A, B and C, at every n; permuting one side of a pair\n",
  "   independent of the other leaves it independent, so this is structural.\n",
  "   Individual verdicts still flip: decisions_changed runs %.3f to %.3f in\n",
  "   the null cells, about half of them a rejection the correct labels would\n",
  "   not have given. An unchanged rate is not a clean bill for one dataset.\n"),
  dc[1], dc[2]))
mild <- resE[which.max(resE$share_correctly_paired), ]
cat(sprintf(paste0(
  "2. A true association -> weakened in proportion to how scrambled the\n",
  "   ordering was. With p1..pN in collection order detection collapses to\n",
  "   alpha (%.3f at rho = 0.5, N = 50, against %.3f correctly labelled) and\n",
  "   the %.0f%% of runs that fail to reject are false negatives -- but that\n",
  "   scheme is close to the worst case. E: with %.0f%% of pairs still correct,\n",
  "   detection is %.3f. Severity is the reader's own ordering, not a constant.\n"),
  ord$detected_affected, ord$detected_correct, 100 * (1 - ord$detected_affected),
  100 * mild$share_correctly_paired, mild$detected_affected))
cat("3. A covariate tracking labelling order -> the residual can be either\n")
cat("   sign: an effect can survive attenuated, or come back reversed.\n")
cat(sprintf(paste0(
  "4. D shows the direction depends on the design, and cuts both ways. With\n",
  "   contiguous condition blocks the correct pairing already maximises the\n",
  "   trend contrast, so the permutation can only reduce it. With alternating\n",
  "   assignment the correct pairing balances the trend across conditions and\n",
  "   the permutation unbalances it, RAISING hypothesis-consistent rejections\n",
  "   above the correctly labelled analysis in %d of %d cells -- %.3f to %.3f\n",
  "   at n = %d, trend = %.1f. So the shuffle can hand back apparent support.\n"),
  nrow(amp), nrow(alt), worst$correct_predicted, worst$mislabelled_predicted,
  worst$n, worst$trend))
cat("   The sign it pushes is set by the ID and assignment schemes, and nothing\n")
cat("   in lexical sorting refers to a hypothesis -- but how often that sign\n")
cat("   agrees with a prediction is NOT measured here: every cell fixes trend\n")
cat("   and predicted direction as positive rather than sampling them. Either\n")
cat("   way an affected analysis has to be re-run before it means anything.\n")
