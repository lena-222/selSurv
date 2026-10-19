test_that("debiased bootstrap prepares the debiased object automatically", {
  local_mocked_bindings(
    coxlasso_debias = function(object) {
      structure(
        list(source = object),
        class = "selSurv_debias"
      )
    },
    .package = "selSurv"
  )

  fit <- structure(list(id = "fit"), class = "selsurv_fit")
  result <- .selsurv_prepare_boot_object(fit, inference = "debiased")

  expect_s3_class(result, "selSurv_debias")
  expect_identical(attr(result, "selsurv_fit"), fit)
})


test_that("an existing debiased object is reused", {
  debiased <- structure(
    list(id = "debiased"),
    class = "selSurv_debias"
  )

  result <- .selsurv_prepare_boot_object(
    debiased,
    inference = "debiased"
  )

  expect_identical(result, debiased)
})


test_that("conditional bootstrap recovers the original fit when available", {
  fit <- structure(list(id = "fit"), class = "selsurv_fit")
  debiased <- structure(
    list(id = "debiased"),
    class = "selSurv_debias"
  )
  attr(debiased, "selsurv_fit") <- fit

  result <- .selsurv_prepare_boot_object(
    debiased,
    inference = "conditional"
  )

  expect_identical(result, fit)
})


test_that("conditional bootstrap rejects a detached debiased object", {
  debiased <- structure(
    list(id = "debiased"),
    class = "selSurv_debias"
  )

  expect_error(
    .selsurv_prepare_boot_object(
      debiased,
      inference = "conditional"
    ),
    "requires the original fitted object"
  )
})


test_that("v12 debiased objects are reused by unified dispatch", {
  fit <- structure(list(id = "fit"), class = "selSurv_fit")
  debiased <- structure(list(fit = fit), class = "selSurv_debias")

  expect_identical(
    .selsurv_prepare_boot_object(debiased, inference = "debiased"),
    debiased
  )
  expect_identical(
    .selsurv_prepare_boot_object(debiased, inference = "conditional"),
    fit
  )
})

test_that("conditional unified API cannot silently invoke reselection", {
  fit <- structure(list(id = "fit"), class = "selSurv_fit")

  expect_error(
    .selsurv_run_boot_engine(
      object = fit, inference = "conditional", type = "model_based",
      B = 2L, conf_level = 0.90, seed = 1L,
      dots = list(workflow = "reselect_fixed_lambda")
    ),
    "fixed-support"
  )
})
