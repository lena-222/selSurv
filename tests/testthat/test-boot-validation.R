test_that("B is validated", {
  expect_identical(.selsurv_validate_B(500), 500L)

  expect_error(.selsurv_validate_B(0), "positive integer")
  expect_error(.selsurv_validate_B(2.5), "positive integer")
  expect_error(.selsurv_validate_B(NA_real_), "positive integer")
})


test_that("confidence level is validated", {
  expect_identical(.selsurv_validate_conf_level(0.95), 0.95)

  expect_error(.selsurv_validate_conf_level(0), "between 0 and 1")
  expect_error(.selsurv_validate_conf_level(1), "between 0 and 1")
  expect_error(.selsurv_validate_conf_level(NA_real_), "between 0 and 1")
})


test_that("seed is validated", {
  expect_null(.selsurv_validate_seed(NULL))
  expect_identical(.selsurv_validate_seed(2026), 2026L)

  expect_error(.selsurv_validate_seed(-1), "non-negative integer")
  expect_error(.selsurv_validate_seed(1.5), "non-negative integer")
})


test_that("dots must be named and cannot duplicate public arguments", {
  expect_invisible(.selsurv_validate_dots(list(engine_option = TRUE)))

  expect_error(
    .selsurv_validate_dots(list(TRUE)),
    "must be named"
  )

  expect_error(
    .selsurv_validate_dots(list(B = 100)),
    "reserved"
  )
})

