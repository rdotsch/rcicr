# Regenerates the `exclusions:` block in .lintr from a fresh lint_package()
# run, ignoring whatever exclusions are already there -- so this baselines
# lints, it never accumulates them. Rerun after:
#   - fixing a lint (the entry drops out on its own)
#   - deliberately accepting a new one (it gets added)
#   - changing the `linters:` block in .lintr
#
# Exclusions are written per file *and linter* as `Inf`, never as line numbers.
# Line numbers do not survive an edit above them: adding four lines to
# R/generateCI.R once shifted seven of them and failed CI in a function ninety
# lines below the change (issue #250).
#
# Needs the package installed, not just load_all()-ed -- object_usage_linter
# requires it (?lintr::executing_linters).
#
# Usage: Rscript tools/regenerate-lintr-baseline.R

library(lintr)

# Linters whose hits mean something is wrong rather than merely unfashionable:
# an undefined name, code left commented out. Excluding one of these for a
# whole file would hide the next real instance in it, so they are annotated at
# the offending line with `# nolint: <linter>.` instead, and this script
# refuses to baseline them. An annotated line never reaches lint_package()
# output, so in normal use this list is inert -- it only fires when a new
# semantic lint appears and someone reruns the script instead of annotating it.
must_annotate <- c(
  "commented_code_linter",
  "object_length_linter",
  "object_name_linter",
  "object_usage_linter",
  "seq_linter"
)

options(lintr.exclusions = list())
lints <- lint_package(".")

df <- as.data.frame(lints)

offending <- df[df$linter %in% must_annotate, ]
if (nrow(offending) > 0) {
  lines <- sprintf(
    "  %s:%d  add  # nolint: %s.",
    offending$filename, offending$line_number, offending$linter
  )
  stop(paste0(
    "These lints must be annotated inline, not baselined -- a whole-file\n",
    "exclusion for them would hide the next real one in the same file:\n",
    paste(lines, collapse = "\n"),
    "\nAnnotate each line above, then rerun this script."
  ), call. = FALSE)
}

format_file_block <- function(filename, rows) {
  linters <- sort(unique(rows$linter))
  linter_entries <- sprintf("      %s = Inf", linters)
  sprintf('    "%s" = list(\n%s\n    )', filename, paste(linter_entries, collapse = ",\n"))
}

by_file <- split(df, df$filename)
by_file <- by_file[sort(names(by_file))]
file_blocks <- vapply(names(by_file), function(f) format_file_block(f, by_file[[f]]), character(1))

exclusions_block <- if (length(file_blocks) == 0L) {
  # The empty case needs the one-line form. A multi-line `list(\n\n  )` does not
  # parse, and lint_package() then aborts on the config instead of reporting
  # zero -- a .lintr that errors is worse than a stale one.
  "exclusions: list()"
} else {
  paste0(
    "exclusions: list(\n",
    paste(file_blocks, collapse = ",\n"),
    # No trailing newline: writeLines() adds one per element, and a second left
    # a blank line at EOF that pre-commit.ci then stripped in its own commit.
    "\n  )"
  )
}

lintr_lines <- readLines(".lintr")
exclusions_start <- grep("^exclusions:", lintr_lines)
stopifnot(length(exclusions_start) == 1)
new_lintr <- c(lintr_lines[seq_len(exclusions_start - 1L)], exclusions_block)
writeLines(new_lintr, ".lintr")

cat(nrow(df), "lints baselined across", length(by_file), "files\n")
