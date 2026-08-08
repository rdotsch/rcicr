# CRAN review of the 1.2.1 submission

**Received** 2026-08-07, from Konstanze Lauseker (she/her), CRAN.
**Submission** `rcicr_1.2.1.tar.gz`, sent 2026-07-29. A request for changes, not a
rejection on the merits — the auto-check had already passed clean.

**Answered by** 1.2.2 (tagged, never submitted) and then 1.2.3. The point-by-point response
is `cran-comments.md`, which is the text pasted into the "Optional comment" field at
<https://cran.r-project.org/submit.html>; it is `.Rbuildignore`d and never travels inside
the tarball.

**Why this file exists.** Point 4 below lists *two* `.Rd` files. Three successive drafts of
our response addressed only the second and told the reviewer we could not find what she had
pointed at — while the fix for the first sat in the same commit, made under a different
point. Working from a summary of the review rather than its text cost that, three times.
Keep the full text so the next round is answered against what was written.

---

## The review, verbatim

Thanks,

Please omit the redundant "Functions to" at the start of your description.

If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
authors (year) <doi:...>
authors (year, ISBN:...)
or if those are not available: <https:...>
with no space after 'doi:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")
For more details: <https://contributor.r-project.org/cran-cookbook/description_issues.html#references>

Please write TRUE and FALSE instead of T and F. Please don't use "T" or "F" as vector names.
For more details: <https://contributor.r-project.org/cran-cookbook/code_issues.html#tf-instead-of-truefalse>

-> Warning: 'T' and 'F' instead of TRUE and FALSE:
  man/generateCI.Rd:
    generateCI(
      stimuli,
      responses,
      baseimage,
      rdata,
      participants = NA,
      save_individual_cis = FALSE,
      save_as_png = TRUE,
      filename = "",
      targetpath = "./cis",
      antiCI = FALSE,
      scaling = "independent",
      scaling_constant = 0.1,
      individual_scaling = "independent",
      individual_scaling_constant = 0.1,
      zmap = F,
      zmapmethod = "quick",
      zmapdecoration = T,
      sigma = 3,
      threshold = 3,
      zmaptargetpath = "./zmaps",
      n_cores = default_ncores(),
      mask = NA
    )
  man/plotZmap.Rd:
    plotZmap(
      zmap,
      bgimage = "",
      sigma,
      threshold = 3,
      mask = NULL,
      decoration = T,
      targetpath = "zmaps",
      filename = "zmap",
      size = 512,
      ...
    )

Some code lines in examples are commented out. Please never do that. Ideally find toy examples that can be regularly executed and checked. Lengthy examples (> 5 sec), can be wrapped in \donttest{}.

-> Examples in comments in:
      generateNoiseImage.Rd
      generateNoisePattern.Rd

\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Does not seem necessary. Please replace \dontrun with \donttest.

Please unwrap the examples if they are executable in < 5 sec, or replace dontrun{} with \donttest{}.

For more details: <https://contributor.r-project.org/cran-cookbook/general_issues.html#structuring-of-examples>

Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). This is not allowed by CRAN policies. Please omit any default path in writing functions. In your examples/vignettes/tests you can write to tempdir().
For more details: <https://contributor.r-project.org/cran-cookbook/code_issues.html#writing-files-and-directories-to-the-home-filespace>

Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.
e.g.:
...
oldpar <- par(no.readonly = TRUE) # code line i
on.exit(par(oldpar)) # code line i + 1
...
par(mfrow=c(2,2)) # somewhere after
...

...
oldwd <- getwd() # code line i
on.exit(setwd(oldwd)) # code line i+1
...
setwd(...) # somewhere after
...
e.g.:
If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.
For more details: <https://contributor.r-project.org/cran-cookbook/code_issues.html#change-of-options-graphical-parameters-and-working-directory>

-> R/plotZmap.R

Please always make sure to reset to user's options(), working directory or par() after you changed it in examples and vignettes and demos.
e.g.:
oldpar <- par(mfrow = c(1,2))
...
par(oldpar)

old <- options(digits = 3)
...
options(old)
For more details: <https://contributor.r-project.org/cran-cookbook/code_issues.html#change-of-options-graphical-parameters-and-working-directory>

-> inst/doc/reverse-correlation-walkthrough.R, please reset the par() like shown above, on.exit() should only be used within functions.

Please fix and resubmit.

Best,
Konstanze Lauseker (she/her)

---

## Reading notes

Things in the text above that a summary loses:

- **Point 4 names two files**, `generateNoiseImage.Rd` *and* `generateNoisePattern.Rd`. Only
  the first had commented-out examples; the second already held a live
  `generateNoisePattern(256)`. Do not read it as naming one.
- **Point 3 is two instructions**: replace `T`/`F`, *and* do not use `T` or `F` as variable
  names. The second is easy to answer and easy to skip.
- **The two `par()` paragraphs ask for different things.** Inside functions she wants an
  immediate `on.exit()`; for examples and vignettes she explicitly says *"on.exit() should
  only be used within functions"* and shows a plain `oldpar <- par(...)` / `par(oldpar)`
  pair. `plotZmap()` and the vignette differing is what she asked for, not an inconsistency.
- **Her `on.exit()` example uses `par(no.readonly = TRUE)`.** We deliberately do not, in
  package code, because it also captures derived parameters such as `pin` that a later
  `plot.window()` invalidates, so the restore itself errors. Both forms satisfy the request;
  `cran-comments.md` no longer spends a paragraph explaining the choice.
- **She quotes the generated `\usage` blocks**, not the roxygen sources — which is why a
  sweep of `R/` alone can report clean while the `.Rd` files she reads still show the
  problem.
