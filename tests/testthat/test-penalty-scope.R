make_penalty_scope_data <- function(seed = 9016L) {
  set.seed(seed)

  n <- 140L
  X <- matrix(rnorm(n * 4L), nrow = n, ncol = 4L)
  colnames(X) <- paste0("x", seq_len(ncol(X)))

  eta <- 0.7 * X[, 1] - 0.5 * X[, 2]
  event_time <- rexp(n, rate = 0.08 * exp(eta))
  censor_time <- rexp(n, rate = 0.05)

  list(
    X = X,
    time = pmin(event_time, censor_time),
    status = as.integer(event_time <= censor_time)
  )
}


capture_scope_warnings <- function(code) {
  messages <- character()

  value <- withCallingHandlers(
    code,
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, messages = messages)
}


test_that("debiased inference rejects non-Lasso pilots", {
  dat <- make_penalty_scope_data()

  fit_en <- coxlasso_fit(
    dat$X,
    dat$time,
    dat$status,
    nfolds = 5L,
    alpha = 0.5
  )

  expect_error(
    coxlasso_debias(fit_en),
    "alpha = 1",
    fixed = TRUE
  )

  expect_error(
    coxlasso_boot(
      fit_en,
      inference = "debiased",
      B = 1L,
      seed = 9016L
    ),
    "alpha = 1",
    fixed = TRUE
  )
})


test_that("score bootstraps defensively reject non-Lasso debiased objects", {
  fake_fit <- list(alpha = 0.5)
  fake_db <- structure(
    list(fit = fake_fit),
    class = "selSurv_debias"
  )

  expect_error(
    coxlasso_wildboot(fake_db, B = 1L),
    "alpha = 1",
    fixed = TRUE
  )

  expect_error(
    coxlasso_efron_events(fake_db, B = 1L),
    "alpha = 1",
    fixed = TRUE
  )
})


test_that("Lasso debiasing remains available", {
  dat <- make_penalty_scope_data(9017L)

  fit_lasso <- coxlasso_fit(
    dat$X,
    dat$time,
    dat$status,
    nfolds = 5L,
    alpha = 1
  )

  db <- coxlasso_debias(fit_lasso)

  expect_s3_class(db, "selSurv_debias")
  expect_identical(fit_lasso$alpha, 1)
})


test_that("split-sample non-Lasso penalty warns but computes", {
  dat <- make_penalty_scope_data(9018L)

  captured <- capture_scope_warnings(
    coxlasso_split(
      dat$X,
      dat$time,
      dat$status,
      seed = 9018L,
      nfolds = 5L,
      penalty_alpha = 0.5
    )
  )

  expect_true(
    any(grepl(
      "inferential validation currently concerns Lasso",
      captured$messages,
      fixed = TRUE
    ))
  )

  expect_s3_class(captured$value, "selSurv_split")
  expect_equal(captured$value$penalty_alpha, 0.5)
  expect_equal(captured$value$fit_selection$alpha, 0.5)
})


test_that("fixed-support model-based bootstrap warns for non-Lasso pilot", {
  dat <- make_penalty_scope_data(9019L)

  fit_en <- coxlasso_fit(
    dat$X,
    dat$time,
    dat$status,
    nfolds = 5L,
    alpha = 0.5
  )

  captured <- capture_scope_warnings(
    coxlasso_boot_parametric(
      fit_en,
      B = 1L,
      workflow = "fixed_support",
      seed = 9019L
    )
  )

  expect_true(
    any(grepl(
      "inferential validation currently concerns Lasso",
      captured$messages,
      fixed = TRUE
    ))
  )

  expect_s3_class(captured$value, "selSurv_boot")
})


test_that("procedure-level re-selection warns for non-Lasso penalty", {
  dat <- make_penalty_scope_data(9020L)

  captured <- capture_scope_warnings(
    coxlasso_boot_reselect(
      dat$X,
      dat$time,
      dat$status,
      B = 1L,
      methods = "fixed_support",
      seed = 9020L,
      nfolds = 5L,
      alpha = 0.5,
      progress = FALSE
    )
  )

  expect_true(
    any(grepl(
      "inferential validation currently concerns Lasso",
      captured$messages,
      fixed = TRUE
    ))
  )

  expect_s3_class(captured$value, "selsurv_reselection")
})


test_that("Lasso split-sample route does not emit penalty-scope warning", {
  dat <- make_penalty_scope_data(9021L)

  captured <- capture_scope_warnings(
    coxlasso_split(
      dat$X,
      dat$time,
      dat$status,
      seed = 9021L,
      nfolds = 5L,
      penalty_alpha = 1
    )
  )

  expect_false(
    any(grepl(
      "inferential validation currently concerns Lasso",
      captured$messages,
      fixed = TRUE
    ))
  )
})
