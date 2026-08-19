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
#
# To measure one specific analysis rather than read a row off scenario E:
#   Rscript tools/simulate-mislabelling-impact.R --ids p1,p2,p10,p3 [rho] [nsim]
#   Rscript tools/simulate-mislabelling-impact.R --ids ids.txt [rho] [nsim]
# All-numeric identifiers need a fifth argument, numeric or character, saying
# how the participants vector stored them -- the sort, and so the answer,
# differs between the two.
#
# This mode cannot represent a FACTOR participants vector with a custom level
# order. The bug ordered participants by factor(participants), which keeps a
# factor's own levels, so factor(c("p2","p1"), levels = c("p2","p1")) was
# correctly labelled while these identifiers passed as text look swapped.
# Reconstructing from strings loses the levels and there is nowhere to put
# them, so that case has to be checked in R against the original object, where
# sort() follows the levels for you:
#   appearance <- unique(participants)
#   mean(appearance == sort(unique(participants)))
# The identifiers go in the order they were passed to generateCI(). This runs
# the caller's own permutation, which E's rows cannot substitute for: the share
# of correctly paired files fixes the average attenuation but not the whole
# distribution, and two orderings sharing that share can differ by around 0.03
# in detection depending on their cycle structure.

args <- commandArgs(trailingOnly = TRUE)
ids_mode <- length(args) >= 2 && args[1] == "--ids"
nsim_arg <- if (ids_mode) args[4] else args[1]
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

# Reader-facing mode: measure the permutation one real analysis actually got,
# rather than approximating it from the share of correctly paired files.
if (ids_mode) {
  raw <- args[2]
  # A path that does not exist must be an error, never data: silently treating
  # "ids.txt" as a single identifier yields one participant, which reports
  # "correctly paired" and hands an affected analysis a false all-clear.
  # A comma settles it: no filename this mode accepts contains one, while an
  # identifier list may well end in something filename-shaped ("subject.dat")
  # or contain a slash, so the path heuristic must not see the raw argument
  # before the list form has been ruled out.
  looks_like_path <- !grepl(",", raw) &&
    (grepl("[/\\\\]", raw) || grepl("\\.(txt|csv|tsv|dat)$", raw))
  if (looks_like_path && !file.exists(raw)) {
    stop("no such file: ", raw, call. = FALSE)
  }
  # Whitespace is only noise in the comma-separated form, where "p1, p2" is
  # the natural way to type it. In a file it may be data: c("p1 ", "p1") are
  # two participants to factor(), and trimming would merge them and could
  # report the rest as correctly paired. So file input keeps the bytes and
  # says when it found any, since stray spaces from an editor would otherwise
  # split one participant into two and raise a false alarm instead.
  from_file <- file.exists(raw)
  ids <- if (from_file) readLines(raw) else trimws(strsplit(raw, ",")[[1]])
  # An empty line is ambiguous and both silent readings are wrong. Dropping it
  # changes the permutation if the original vector really held an empty
  # identifier; keeping it invents a participant on any file whose last line is
  # blank, which is most of them -- writeLines() and most editors leave one. So
  # neither: say what was found and let the caller resolve it.
  if (any(!nzchar(ids))) {
    where <- paste(which(!nzchar(ids)), collapse = ", ")
    stop(if (from_file) paste0("blank line(s) in ", raw, ": lines ", where)
         else paste0("empty field(s) in the identifier list: position ", where),
         ".\n",
         if (from_file) {
           "  If they are formatting -- a trailing newline, say -- remove them and re-run.\n"
         } else {
           "  If that is a stray comma, remove it and re-run.\n"
         },
         "  If the participants vector really held an empty identifier, this ",
         "mode cannot\n  represent it; check that case in R against the ",
         "original object.", call. = FALSE)
  }
  if (from_file && any(ids != trimws(ids))) {
    cat("NOTE: some identifiers in this file carry leading or trailing space,\n")
    cat("  kept as given because factor() would have treated them as distinct.\n")
    cat("  If they are an artefact of how the file was written rather than of\n")
    cat("  the original participants vector, strip them and re-run.\n")
  }
  if (length(unique(ids)) < 2) {
    stop("need at least two distinct participant identifiers, got ",
         length(unique(ids)), call. = FALSE)
  }
  rho <- if (length(args) >= 3 && !is.na(args[3])) as.numeric(args[3]) else 0.5

  # Type decides the sort, so it cannot be guessed. R sorts a numeric
  # participants vector numerically, leaving 1:12 unaffected, while the
  # character vector "1" ... "12" sorts lexically and is affected. Everything
  # arriving on a command line is a string, so an all-numeric list is ambiguous
  # and the caller has to say which they had: guessing "numeric" would tell a
  # character-ID reader their filenames were fine, and as.numeric() would also
  # collapse "01" and "1" into one participant after uniqueness was checked.
  numeric_looking <- !any(is.na(suppressWarnings(as.numeric(ids))))
  if (numeric_looking) {
    type <- if (length(args) >= 5) args[5] else NA_character_
    if (is.na(type) || !type %in% c("numeric", "character")) {
      stop("these identifiers are all numbers, and the answer depends on how ",
           "they were stored.\n",
           "  participants = c(1, 2, ...)       -> sorted numerically, ",
           "usually unaffected\n",
           "  participants = c(\"1\", \"2\", ...) -> sorted lexically, ",
           "affected from the tenth on\n",
           "  Re-run with a fifth argument, numeric or character, e.g.\n",
           "    --ids ", args[2], " ", if (is.na(args[3])) "0.5" else args[3],
           " ", if (is.na(args[4])) "20000" else args[4], " character",
           call. = FALSE)
    }
    if (type == "numeric") {
      converted <- as.numeric(ids)
      if (length(unique(converted)) != length(unique(ids))) {
        stop("reading these as numbers merges identifiers that differ as text ",
             "(\"01\" and \"1\", say); they must have been character",
             call. = FALSE)
      }
      ids <- converted
      cat("identifiers read as numbers, matching participants = c(1, 2, ...)\n")
    } else {
      cat("identifiers read as character strings, matching ",
          "participants = c(\"1\", \"2\", ...)\n", sep = "")
    }
  }
  n <- length(unique(ids))
  perm <- mislabel_perm(ids)
  kept <- mean(perm == seq_len(n))
  cat(sprintf("%d participants, %.3f of files correctly paired\n", n, kept))
  cat("  (assumes participants was a plain vector; a factor with a custom level\n")
  cat("   order sorts by its levels, which this mode cannot represent -- check\n")
  cat("   that case in R against the original object)\n")
  # This mode re-derives the ordering with sort(), so it inherits the running
  # session's collation rather than the one the affected analysis ran under.
  # Only identifiers outside plain lowercase ASCII can order differently, so
  # the warning fires for those rather than on every run.
  if (any(grepl("[^0-9a-z]", as.character(ids)))) {
    cat("  WARNING: these identifiers contain characters whose sort order is\n")
    cat("   locale-dependent (upper case, accents or punctuation), and this run\n")
    cat(sprintf("   used LC_COLLATE=%s. If the original analysis ran under a\n",
                Sys.getlocale("LC_COLLATE")))
    cat("   different collation, the permutation below is not the one it got.\n")
  }
  if (kept == 1) {
    cat("Nothing to do: this ordering was already sorted, so the filenames were correct.\n")
  } else if (n < 3) {
    # cor.test() needs three finite pairs. The pairing share above is still
    # exact and is the part that describes the mislabelling; only the detection
    # rates need a sample to estimate.
    cat("Two participants: the share above is exact, but detection rates need at\n")
    cat("  least three, so they are not estimated here. Both files carry the\n")
    cat("  other participant's image; the recovery in the advisory still applies.\n")
  } else {
    hit_c <- hit_a <- logical(NSIM)
    for (i in seq_len(NSIM)) {
      x <- rnorm(n)
      y <- rho * x + sqrt(1 - rho^2) * rnorm(n)
      hit_c[i] <- cor.test(x, y)$p.value < ALPHA
      hit_a[i] <- cor.test(x, y[perm])$p.value < ALPHA
    }
    cat(sprintf(paste0(
      "At rho = %.2f, a correctly labelled analysis of this size detects the\n",
      "association %.3f of the time; with this ordering's mislabelling, %.3f.\n",
      "Verdicts differ on %.3f of datasets.\n"),
      rho, mean(hit_c), mean(hit_a), mean(hit_c != hit_a)))
  }
  quit(save = "no")
}


# Not monotone in n: the count turns on where the identifiers change width,
# which is why 99-101 keep eleven and 120 is back to one. Readers are given
# mean(sort(unique(participants)) == unique(participants)) rather than a rule.
cat("=== How much of the pairing survives ===\n")
for (n in c(9, 12, 20, 30, 50, 99, 100, 120)) {
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
             flip_to_sig = sum(!hit_c & hit_a) / nsim,
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
  data.frame(n = n, d = d, flip_to_sig = sum(!hit_c & hit_a) / nsim,
             detected_correct = mean(hit_c),
             detected_affected = mean(hit_a),
             diff = mean(hit_a) - mean(hit_c), diff_se = sqrt(discordant) / nsim,
             decisions_changed = discordant / nsim, sig_predicted = mean(sig_p),
             sig_opposite = mean(sig_o), crossed_condition = mean(cond[perm] != cond))
}
gridB <- expand.grid(n = c(12, 20, 30, 50), d = c(0, 0.8))
resB <- do.call(rbind, Map(function(n, d) run_blocked(n, d), gridB$n, gridB$d))
print(resB, row.names = FALSE, digits = 3)

# The null cells are the evidence for the unchanged false positive rate, so they
# get their own table. Correct and mislabelled analyses see the same outcome
# vector within an iteration, so the paired difference is the statistic: its
# error is much smaller than either rate's, and it does not confound the
# permutation with the test's own size. That distinction matters here --
# Welch is mildly conservative at n = 12, visible in B's *correct* column.
cat("\n=== Null cells, paired: is the mislabelled rate the same as the correct one? ===\n")
nullcells <- rbind(
  data.frame(design = "A", test = "cor.test", resA[resA$rho == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se", "decisions_changed", "flip_to_sig")]),
  data.frame(design = "C", test = "cor.test", resC[resC$rho == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se", "decisions_changed", "flip_to_sig")]),
  data.frame(design = "B", test = "Welch t", resB[resB$d == 0, c("n", "detected_correct", "detected_affected", "diff", "diff_se", "decisions_changed", "flip_to_sig")])
)
nullcells$SEs_from_zero <- abs(nullcells$diff) / nullcells$diff_se
print(nullcells, row.names = FALSE, digits = 3)

# D asks whether the permutation can hand back a hypothesis-consistent result
# rather than only destroy one. Scenario B at d = 0 already answers it for an
# exchangeable outcome, where the rate is unchanged although individual
# verdicts still flip both ways. The harder case is an outcome carrying real
# structure tied to collection order -- drift, practice, season -- with no true
# condition effect, where the permutation has something to move around.
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
# A cell with no discordant decisions -- reachable at small nsim, which the
# header encourages -- would give 0/0 and print NA% for the whole range.
observed <- nullrows[nullrows$decisions_changed > 0, ]
share_to_sig <- observed$flip_to_sig / observed$decisions_changed
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
  "   the null cells%s.\n",
  "   An unchanged rate is not a clean bill for one dataset.\n"),
  dc[1], dc[2],
  if (length(share_to_sig) == 0) {
    ", none of which changed any decision at this nsim"
  } else {
    sprintf(", of which %.0f%% to %.0f%% are rejections\n   the correct labels would not have given",
            100 * min(share_to_sig), 100 * max(share_to_sig))
  }))
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
