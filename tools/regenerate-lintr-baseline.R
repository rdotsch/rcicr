# Regenerates the `exclusions:` block in .lintr from a fresh lint_package()
# run, ignoring whatever exclusions are already there -- so this baselines
# lints, it never accumulates them. Rerun after:
#   - fixing a lint (the entry drops out on its own)
#   - deliberately accepting a new one (it gets added)
#   - changing the `linters:` block in .lintr
#
# Needs the package installed, not just load_all()-ed -- object_usage_linter
# requires it (?lintr::executing_linters).
#
# Usage: Rscript tools/regenerate-lintr-baseline.R

library(lintr)

options(lintr.exclusions = list())
lints <- lint_package(".")

df <- as.data.frame(lints)

format_file_block <- function(filename, rows) {
  by_linter <- split(rows$line_number, rows$linter)
  by_linter <- by_linter[sort(names(by_linter))]
  linter_entries <- vapply(names(by_linter), function(linter) {
    line_numbers <- paste(sort(unique(by_linter[[linter]])), collapse = ", ")
    sprintf("      %s = c(%s)", linter, line_numbers)
  }, character(1))
  sprintf('    "%s" = list(\n%s\n    )', filename, paste(linter_entries, collapse = ",\n"))
}

by_file <- split(df, df$filename)
by_file <- by_file[sort(names(by_file))]
file_blocks <- vapply(names(by_file), function(f) format_file_block(f, by_file[[f]]), character(1))

exclusions_block <- paste0(
  "exclusions: list(\n",
  paste(file_blocks, collapse = ",\n"),
  "\n  )\n"
)

lintr_lines <- readLines(".lintr")
exclusions_start <- grep("^exclusions:", lintr_lines)
stopifnot(length(exclusions_start) == 1)
new_lintr <- c(lintr_lines[seq_len(exclusions_start - 1L)], exclusions_block)
writeLines(new_lintr, ".lintr")

cat(length(lints), "lints baselined across", length(by_file), "files\n")
