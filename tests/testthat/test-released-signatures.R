# The guiding constraint, made mechanical: researchers re-run stored scripts,
# and some of those scripts pass arguments positionally. A new formal inserted
# anywhere but the end rebinds every argument after it -- silently, with no
# error, at whatever the caller's next value happens to be.
#
# This is not hypothetical. `generateCI(zmappointsize=)` was first added before
# `sigma` (#212), where a positional call would have handed sigma's value to it
# and blurred the z-map at the default instead.
#
# The baseline is v1.0.1, the same version the release gate pins for numbers:
# the CRAN-archived release most stored scripts were written against. Measured
# before choosing it -- all 17 exported functions still begin with their v1.0.1
# arguments, so the stricter baseline costs no exception list. 1.2.2's breaking
# change removed argument *defaults*, which does not move positional bindings
# and so does not show up here; names are what a positional call depends on.
#
# To regenerate after a deliberate change:
#
#   mkdir -p /tmp/rel101
#   for f in $(git ls-tree --name-only v1.0.1 R/); do
#     git show v1.0.1:$f > /tmp/rel101/$(basename $f); done
#   # source them into an environment, then
#   # saveRDS(lapply(fns, function(f) names(formals(f))), fixture_path)
#
# Regenerating is the deliberate act. Doing it to make this test pass is how a
# reordered signature reaches a user.

test_that("every exported function still starts with its released arguments", {
  released <- readRDS(test_path("fixtures", "released-formals-1.0.1.rds"))
  exports <- sort(getNamespaceExports("rcicr"))
  compared <- 0

  for (nm in exports) {
    old <- released[[nm]]
    if (is.null(old)) next  # added after 1.0.1; nothing released to preserve
    new <- names(formals(get(nm, envir = asNamespace("rcicr"))))
    compared <- compared + 1

    # A prefix, not an equality: appending an optional argument is allowed, and
    # is how every argument added since 1.0.1 went in.
    expect_equal(
      new[seq_along(old)], old,
      info = paste0(
        nm, "() no longer starts with the arguments it shipped with in 1.0.1. ",
        "Released: ", paste(old, collapse = ", "), ". ",
        "Now: ", paste(new, collapse = ", "), ". ",
        "Append new arguments instead of inserting them."
      )
    )
  }

  # Guards the guard: a fixture that failed to load, or an export list that came
  # back empty, would otherwise pass this file in silence.
  expect_gt(compared, 10)
})
