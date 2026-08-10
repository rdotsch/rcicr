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

# CFF 1.2.0's required keys. GitHub renders the "Cite this repository" button
# from this file and silently declines when it cannot parse it, so a structural
# break has no other symptom.
required <- c("cff-version", "message", "title", "authors")
missing <- required[!required %in% names(cff)]
if (length(missing) || !length(cff$authors)) {
  cat("CITATION.cff is not valid CFF: missing ",
    paste(c(missing, if (!length(cff$authors)) "a non-empty authors list"), collapse = ", "),
    "\n",
    sep = ""
  )
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
