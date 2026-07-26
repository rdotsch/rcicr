test_that("simulateNoiseIntensities currently errors (known bug, not fixed here)", {
  # R/simulateNoiseIntensities.R references `data[,by]` to size its progress
  # bar, but neither `data` nor `by` are parameters of this function (its only
  # parameters are `nrep` and `img_size`). `data` resolves to the imported
  # `utils::data` function object, so `data[,by]` fails immediately with
  # "object of type 'closure' is not subsettable". Documenting current
  # (broken) behavior; not fixed as part of this task - flag as a follow-up
  # issue instead.
  suppressWarnings(expect_error(
    simulateNoiseIntensities(nrep = 2, img_size = 32),
    "not subsettable"
  ))
})
