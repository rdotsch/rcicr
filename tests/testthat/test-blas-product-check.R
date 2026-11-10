test_that("the matrix-product check passes under this session's BLAS", {
  expect_true(blas_product_ok())
})

test_that("attaching reports a wrong matrix product, and only then", {
  local_mocked_bindings(blas_product_ok = function() FALSE)
  expect_message(rcicr:::.onAttach("", "rcicr"), "MKL_THREADING_LAYER=GNU")
  local_mocked_bindings(blas_product_ok = function() TRUE)
  expect_silent(rcicr:::.onAttach("", "rcicr"))
})
