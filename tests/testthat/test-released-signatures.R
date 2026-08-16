# The guiding constraint, made mechanical: researchers re-run stored scripts,
# and some of those scripts pass arguments positionally. A new formal inserted
# anywhere but the end rebinds every argument after it -- silently, with no
# error, at whatever the caller's next value happens to be.
#
# This is not hypothetical. `generateCI(zmappointsize=)` was first added before
# `sigma` (#212), where a positional call would have handed sigma's value to it
# and blurred the z-map at the default instead.
#
# Two baselines, because arguments have been appended across several releases
# and they catch different things:
#
#   1.2.3  the latest release, and the one doing the work. Its argument list
#          contains v1.0.1's as a prefix -- measured, for all 17 exported
#          functions -- so it pins everything the older baseline does *and* the
#          arguments added in 1.1.0-1.2.3. Checking only v1.0.1 would let a new
#          formal be inserted after its arguments but before the later ones,
#          passing while breaking every script written against 1.1.0 or newer.
#   1.0.1  a floor, and never regenerated. It is the CRAN-archived release most
#          stored scripts were written against, and the fixture nobody is
#          tempted to refresh when a check goes red.
#
# Refresh the 1.2.3 fixture at release time, renaming it for the new version;
# leave the 1.0.1 one alone. 1.2.2's breaking change removed argument
# *defaults*, which does not move a positional binding and so does not appear
# here -- names are what a positional call depends on.
#
# To regenerate after a deliberate change:
#
#   mkdir -p /tmp/rel101 # nolint: commented_code_linter.
#   for f in $(git ls-tree --name-only v1.0.1 R/); do
#     git show v1.0.1:$f > /tmp/rel101/$(basename $f); done
#   # source them into an environment, then
#   # saveRDS(lapply(fns, function(f) names(formals(f))), fixture_path) # nolint: commented_code_linter.
#
# Regenerating is the deliberate act. Doing it to make this test pass is how a
# reordered signature reaches a user.

test_that("every exported function still starts with its released arguments", {
  baselines <- c("1.0.1", "1.2.3")
  exports <- sort(getNamespaceExports("rcicr"))
  compared <- 0

  for (version in baselines) {
  released <- readRDS(test_path("fixtures", paste0("released-formals-", version, ".rds")))
  for (nm in exports) {
    old <- released[[nm]]
    if (is.null(old)) next  # added after this release; nothing to preserve
    new <- names(formals(get(nm, envir = asNamespace("rcicr"))))
    compared <- compared + 1

    # A prefix, not an equality: appending an optional argument is allowed, and
    # is how every argument added since 1.0.1 went in.
    #
    expect_equal(
      new[seq_along(old)], old,
      info = paste0(
        nm, "() no longer starts with the arguments it shipped with in ",
        version, ". ",
        "Released: ", paste(old, collapse = ", "), ". ",
        "Now: ", paste(new, collapse = ", "), ". ",
        "Append new arguments instead of inserting them."
      )
    )
  }
  }

  # Guards the guard: a fixture that failed to load, or an export list that came
  # back empty, would otherwise pass this file in silence.
  expect_gt(compared, 10)
})
