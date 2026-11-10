# Temporary diagnostic for #425: why did the attach check stay silent under
# the broken threading setup? Relative error of 10 successive products, for
# the check's fixed matrices and for the old probe's random ones, under
# matprod "default" and "blas".
rel <- function(a, b, mode) {
  options(matprod = mode); got <- a %*% b; options(matprod = "default")
  want <- vapply(seq_len(ncol(b)), function(j) colSums(a * b[, j]), numeric(nrow(a)))
  max(abs(got - want)) / max(abs(want))
}
n <- 40L; k <- 20L
x <- matrix(sin(seq_len(n * n)), n); fixed_a <- x + t(x)
fixed_b <- matrix(ifelse(cos(seq_len(n * k) * 1.7) > 0, 1, -1), n, k)
set.seed(1)
rand_a <- crossprod(matrix(rnorm(n * n), n))
rand_b <- matrix(((runif(n * k) > 0.5) * 2) - 1, n, k)
cat("BLAS:", extSoftVersion()[["BLAS"]], " MKL_THREADING_LAYER:", Sys.getenv("MKL_THREADING_LAYER", "<unset>"), "\n")
for (mode in c("default", "blas")) for (what in c("fixed", "random")) {
  a <- if (what == "fixed") fixed_a else rand_a
  b <- if (what == "fixed") fixed_b else rand_b
  cat(sprintf("%-7s %-6s %s\n", mode, what,
              paste(formatC(vapply(1:10, function(i) rel(a, b, mode), 0), format = "e", digits = 1), collapse = " ")))
}
