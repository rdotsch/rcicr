# A closure that puts the caller's random stream back as it was. Assigning the
# saved .Random.seed back restores kind and position together; an RNGkind()
# round trip does not, since RNGkind() reseeds. Without a seed to assign back,
# the kind lives only in the session, so it is re-applied before the seed that
# call creates is removed.
captureRandomStream <- function() {
  had <- exists('.Random.seed', envir = globalenv(), inherits = FALSE)
  before <- if (had) get('.Random.seed', envir = globalenv(), inherits = FALSE)
  kind <- RNGkind()
  function() {
    if (had) {
      # nolint next: object_name_linter. R owns this name, not this package.
      assign('.Random.seed', before, envir = globalenv())
      return(invisible(NULL))
    }
    # sample.kind = "Rounding" warns on every call that sets it.
    if (!identical(RNGkind(), kind)) suppressWarnings(do.call(RNGkind, as.list(kind)))
    if (exists('.Random.seed', envir = globalenv(), inherits = FALSE)) {
      rm('.Random.seed', envir = globalenv())
    }
    invisible(NULL)
  }
}

# A cache hit consumed no random numbers before automatic migration existed.
preserveRandomStream <- function(expr) {
  restore <- captureRandomStream()
  on.exit(restore(), add = TRUE)
  expr
}
