#!/usr/bin/env Rscript

# CITATION.cff duplicates what DESCRIPTION and inst/CITATION already say, and
# this package is cited in published work -- two citation records that disagree
# are worse than one. So the duplication is enforced rather than trusted: every
# field of CITATION.cff is derived from a source of truth here, and the
# `ubuntu-latest (release)` job fails when one drifts.
#
# `version` and `date-released` come from the newest *release* -- DESCRIPTION's
# version between releases carries a `.9000` suffix, and a citation names a
# released artifact -- so they are checked against NEWS.md's dated heading.

fail <- character()
check <- function(label, want, have) {
  ok <- identical(as.character(want), as.character(have))
  cat(if (ok) "  ok   " else "  FAIL ", label, ": ", sQuote(have),
    if (!ok) paste0(" (expected ", sQuote(want), ")") else "", "\n",
    sep = ""
  )
  if (!ok) fail <<- c(fail, label)
}

cff <- yaml::read_yaml("CITATION.cff")

# GitHub renders the "Cite this repository" button from this file and silently
# declines when it cannot parse it, so a structural break has no other symptom.
# The comparisons below cover every derived field; these are the ones nothing
# derives, so a hand edit to one would otherwise pass unnoticed. Checked by value
# rather than against the CFF schema, which would mean vendoring the schema and
# depending on jsonvalidate for a file this size.
invalid <- character()
note <- function(...) invalid <<- c(invalid, paste0(...))
str1 <- function(x) is.character(x) && length(x) == 1L && nzchar(trimws(x))

for (field in setdiff(c("cff-version", "message", "title", "authors"), names(cff))) {
  note("missing ", field)
}

constants <- c(`cff-version` = "1.2.0", type = "software")
for (field in names(constants)) {
  have <- if (is.null(cff[[field]])) "(absent)" else as.character(cff[[field]])
  if (!identical(have, constants[[field]])) {
    note(field, " must be ", sQuote(constants[[field]]), ", not ", sQuote(have))
  }
}

if (!str1(cff$message)) note("message must be a non-empty string")
if (!length(cff$authors)) {
  note("authors must be a non-empty list")
} else if (!all(vapply(cff$authors, function(a) str1(a$`family-names`) || str1(a$name), TRUE))) {
  note("every author needs family-names (a person) or name (an entity)")
}

if (length(invalid)) {
  cat("CITATION.cff is not valid CFF 1.2.0:\n", paste0("  ", invalid, "\n"), sep = "")
  quit(status = 1L)
}

desc <- read.dcf("DESCRIPTION")[1, ]
cit <- utils::readCitationFile("inst/CITATION", meta = as.list(desc))[[1]]

release <- sub("\\.9000$", "", desc[["Version"]])
check("version (DESCRIPTION, development suffix dropped)", release, cff$version)

# The heading is what makes a version a release; R's own news parser reads these.
news <- readLines("NEWS.md", warn = FALSE)
heading <- grep(paste0("^# rcicr ", release, " \\([0-9]{4}-[0-9]{2}-[0-9]{2}\\)$"),
  news,
  value = TRUE
)
if (length(heading) != 1L) {
  cat("  FAIL date-released: no single '# rcicr ", release,
    " (YYYY-MM-DD)' heading in NEWS.md\n",
    sep = ""
  )
  fail <- c(fail, "date-released")
} else {
  check("date-released (NEWS.md heading)", sub(".*\\((.*)\\)$", "\\1", heading), cff$`date-released`)
}

check("title (inst/CITATION)", cit$title, cff$title)

authors <- cit$author
check(
  "authors (inst/CITATION)",
  paste(vapply(authors, function(a) paste(a$given, a$family), ""), collapse = "; "),
  paste(vapply(cff$authors, function(a) paste(a$`given-names`, a$`family-names`), ""),
    collapse = "; "
  )
)

# DESCRIPTION carries the SPDX-less R spelling; CFF wants SPDX.
spdx <- c("GPL-2" = "GPL-2.0-only", "GPL-3" = "GPL-3.0-only", "MIT + file LICENSE" = "MIT")
check("license (DESCRIPTION, as SPDX)", spdx[[desc[["License"]]]], cff$license)

urls <- sub("/$", "", trimws(strsplit(desc[["URL"]], ",")[[1]]))
for (field in c("repository-code", "url")) {
  have <- sub("/$", "", cff[[field]])
  if (!have %in% urls) {
    cat("  FAIL ", field, ": ", sQuote(have), " is not in DESCRIPTION's URL field\n", sep = "")
    fail <- c(fail, field)
  } else {
    cat("  ok   ", field, ": ", sQuote(have), "\n", sep = "")
  }
}

if (length(fail)) {
  cat("\nCITATION.cff disagrees with the package metadata in: ",
    paste(fail, collapse = ", "), "\n",
    "Fix CITATION.cff, or the source it should have matched.\n",
    sep = ""
  )
  quit(status = 1L)
}
cat("\nCITATION.cff agrees with DESCRIPTION, inst/CITATION and NEWS.md.\n")
