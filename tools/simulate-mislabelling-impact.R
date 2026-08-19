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
# stated. Three designs, at the alpha = 0.05 level:
#
#   A  covariate unrelated to the order participants were labelled in
#      (the ordinary individual-differences design)
#   B  condition assigned in blocks by collection order
#   C  covariate that is collection order itself -- testing date, cohort
#
# Usage: Rscript tools/simulate-mislabelling-impact.R [nsim]
# Runs in well under a minute at the default nsim.

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
    tc <- suppressWarnings(cor.test(x, y))
    ta <- suppressWarnings(cor.test(x, y[perm]))
    hit_c[i] <- tc$p.value < ALPHA
    hit_a[i] <- ta$p.value < ALPHA
    opp[i] <- ta$p.value < ALPHA && ta$estimate < 0
    r_a[i] <- ta$estimate
  }
  data.frame(n = n, rho = rho,
             detected_correct = mean(hit_c), detected_affected = mean(hit_a),
             mean_r_affected = mean(r_a), sig_opposite = mean(opp))
}

grid <- expand.grid(n = c(12, 20, 30, 50), rho = c(0, 0.3, 0.5))

cat("=== A: covariate unrelated to labelling order ===\n")
print(do.call(rbind, Map(function(n, r) run_cor(n, r, function(k) rnorm(k)),
                         grid$n, grid$rho)),
      row.names = FALSE, digits = 3)

cat("\n=== C: covariate is collection order itself (testing date, cohort) ===\n")
print(do.call(rbind, Map(function(n, r) run_cor(n, r, function(k) as.numeric(scale(seq_len(k)))),
                         grid$n, grid$rho)),
      row.names = FALSE, digits = 3)

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
  data.frame(n = n, d = d, detected_correct = mean(hit_c),
             detected_affected = mean(hit_a), sig_predicted = mean(sig_p),
             sig_opposite = mean(sig_o), crossed_condition = mean(cond[perm] != cond))
}
gridB <- expand.grid(n = c(12, 20, 30, 50), d = c(0, 0.8))
print(do.call(rbind, Map(function(n, d) run_blocked(n, d), gridB$n, gridB$d)),
      row.names = FALSE, digits = 3)

cat("\n=== Summary ===\n")
cat("1. No true association -> the mislabelling cannot create one. Rejection\n")
cat("   stays at alpha in all three designs, at every n. Permuting one side of\n")
cat("   a pair that is independent of the other leaves it independent, so this\n")
cat("   is structural rather than a property of these particular draws.\n")
cat("2. A true association, ordinary design -> it is destroyed. Detection\n")
cat("   collapses to alpha, and every one of those is a false negative.\n")
cat("3. A covariate tracking labelling order -> the residual can be either\n")
cat("   sign: an effect can survive attenuated, or come back reversed.\n")
