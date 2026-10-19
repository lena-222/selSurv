test_that("target tables are available", {
  expect_true(is.data.frame(selSurv_targets()))
  expect_true(is.data.frame(selSurv_methods()))
})
