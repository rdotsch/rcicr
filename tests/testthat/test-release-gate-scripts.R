# The release gate (tools/compare-release-output.R and the battery it drives)
# runs at most a handful of times per release, and its failure mode is a syntax
# error discovered at the worst possible moment. Parsing it costs nothing.
#
# tools/ is excluded from the built package via .Rbuildignore, so this can only
# run from a source checkout.

gate_scripts <- c("compare-harness.R", "compare-release-output.R")

test_that("the release-gate scripts parse", {
  tools_dir <- test_path("..", "..", "tools")
  skip_if_not(dir.exists(tools_dir), "tools/ is not present (built package)")

  for (script in gate_scripts) {
    path <- file.path(tools_dir, script)
    expect_true(file.exists(path))
    expect_no_error(parse(path))
  }
})
