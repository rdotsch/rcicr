# Temporary diagnostic for #425: does this BLAS return A %*% B correctly when
# called repeatedly with the same shapes, as rcicr's gramNorms() does?
# Compares each BLAS product with R's internal (BLAS-free) matrix product.
set.seed(1)
report <- function(label, n, k, reps = 200) {
  a <- crossprod(matrix(rnorm(n * n), n))
  bad_cols <- 0L
  ratios <- numeric(0)
  for (r in seq_len(reps)) {
    b <- matrix(((runif(n * k) > 0.5) * 2) - 1, n, k)
    options(matprod = "blas")
    got <- a %*% b
    options(matprod = "internal")
    want <- a %*% b
    options(matprod = "default")
    off <- colSums(abs(got - want)) > 1e-8 * colSums(abs(want))
    if (any(off)) {
      bad_cols <- bad_cols + sum(off)
      ratios <- c(ratios, colSums(got[, off, drop = FALSE] * b[, off, drop = FALSE]) /
                    colSums(want[, off, drop = FALSE] * b[, off, drop = FALSE]))
    }
  }
  cat(sprintf("%-28s n=%4d k=%4d reps=%d  wrong columns: %d of %d", label, n, k, reps,
              bad_cols, reps * k))
  if (length(ratios)) cat(sprintf("  r'Ar ratio: min %.4f median %.4f max %.4f",
                                  min(ratios), median(ratios), max(ratios)))
  cat("\n")
}
cat("BLAS:", extSoftVersion()[["BLAS"]], "\n")
for (n in c(8L, 20L, 30L, 40L, 100L)) for (k in c(3L, 5L, 20L, 1000L)) report("gram %*% responses", n, k)

# tcrossprod(), as in gramFromParameters(), against the internal product.
x <- matrix(rnorm(40 * 500), 40)
options(matprod = "blas"); g1 <- tcrossprod(x)
options(matprod = "internal"); g2 <- tcrossprod(x)
options(matprod = "default")
cat("tcrossprod max relative difference:", max(abs(g1 - g2)) / max(abs(g2)), "\n")
