test_that("fully implicit specification uses debiased wild bootstrap", {
  specification <- .selsurv_resolve_boot_spec()

  expect_identical(specification$inference, "debiased")
  expect_identical(specification$type, "wild")
})


test_that("bootstrap type determines inference family", {
  wild <- .selsurv_resolve_boot_spec(type = "wild")
  both <- .selsurv_resolve_boot_spec(type = c("wild", "efron"))
  model_based <- .selsurv_resolve_boot_spec(type = "model_based")

  expect_identical(wild$inference, "debiased")
  expect_identical(both$inference, "debiased")
  expect_identical(both$type, c("wild", "efron"))
  expect_identical(model_based$inference, "conditional")
})


test_that("inference family determines default bootstrap type", {
  debiased <- .selsurv_resolve_boot_spec(inference = "debiased")
  conditional <- .selsurv_resolve_boot_spec(inference = "conditional")

  expect_identical(debiased$type, "wild")
  expect_identical(conditional$type, "model_based")
})


test_that("type all expands within an explicit inference family", {
  debiased <- .selsurv_resolve_boot_spec(
    inference = "debiased",
    type = "all"
  )
  conditional <- .selsurv_resolve_boot_spec(
    inference = "conditional",
    type = "all"
  )

  expect_identical(debiased$type, c("wild", "efron"))
  expect_identical(conditional$type, "model_based")
})


test_that("mixed inferential families are rejected", {
  expect_error(
    .selsurv_resolve_boot_spec(
      type = c("wild", "model_based")
    ),
    "different inference targets"
  )
})


test_that("explicit incompatible arguments are rejected", {
  expect_error(
    .selsurv_resolve_boot_spec(
      inference = "conditional",
      type = "wild"
    ),
    "not available"
  )

  expect_error(
    .selsurv_resolve_boot_spec(
      inference = "debiased",
      type = "model_based"
    ),
    "not available"
  )
})


test_that("unknown bootstrap types are rejected", {
  expect_error(
    .selsurv_resolve_boot_spec(type = "pairs"),
    "Unknown bootstrap type"
  )
})


test_that("duplicate types are removed while preserving order", {
  specification <- .selsurv_resolve_boot_spec(
    type = c("wild", "efron", "wild")
  )

  expect_identical(specification$type, c("wild", "efron"))
})
