test_that("Cox glmnet defaults are explicit and version-compatible", {
  old_cv <- .selsurv_prepare_cox_glmnet_args(
    list(),
    selection_ties = "breslow",
    cv = TRUE,
    glmnet_formals = character()
  )

  expect_identical(old_cv$type.measure, "deviance")
  expect_identical(old_cv$grouped, TRUE)
  expect_false("cox.ties" %in% names(old_cv))

  expect_error(
    .selsurv_prepare_cox_glmnet_args(
      list(),
      selection_ties = "efron",
      cv = FALSE,
      glmnet_formals = character()
    ),
    "does not support Efron"
  )

  new_breslow <- .selsurv_prepare_cox_glmnet_args(
    list(),
    selection_ties = "breslow",
    cv = TRUE,
    glmnet_formals = "cox.ties"
  )

  expect_identical(new_breslow$cox.ties, "breslow")
  expect_identical(new_breslow$type.measure, "deviance")
  expect_identical(new_breslow$grouped, TRUE)

  new_efron <- .selsurv_prepare_cox_glmnet_args(
    list(),
    selection_ties = "efron",
    cv = FALSE,
    glmnet_formals = "cox.ties"
  )

  expect_identical(new_efron$cox.ties, "efron")
})


test_that("explicit user CV settings remain overridable", {
  args <- .selsurv_prepare_cox_glmnet_args(
    list(
      type.measure = "C",
      grouped = FALSE,
      cox.ties = "efron"
    ),
    selection_ties = "breslow",
    cv = TRUE,
    glmnet_formals = "cox.ties"
  )

  expect_identical(args$type.measure, "C")
  expect_identical(args$grouped, FALSE)
  expect_identical(args$cox.ties, "efron")
})


test_that("coxlasso_fit records the Cox selection and refit conventions", {
  set.seed(9015)

  n <- 140L
  X <- matrix(rnorm(n * 4L), nrow = n)
  colnames(X) <- paste0("x", 1:4)

  eta <- 1.1 * X[, 1] - 0.9 * X[, 2]
  event <- rexp(n, rate = exp(eta))
  censor <- rexp(n, rate = 0.25)

  time <- pmin(event, censor)
  status <- as.integer(event <= censor)

  fit <- coxlasso_fit(
    X,
    time,
    status,
    lambda_choice = "min",
    nfolds = 5
  )

  expect_identical(fit$cv_type_measure, "deviance")
  expect_identical(fit$cv_grouped, TRUE)
  expect_identical(fit$selection_ties, "breslow")
  expect_identical(fit$refit_ties, "efron")

  if (length(fit$support) > 0L) {
    expect_identical(fit$refit$fit$method, "efron")
  }
})


test_that("Breslow bootstrap baseline is explicit for Efron refits", {
  set.seed(90152)

  n <- 160L
  X <- matrix(rnorm(n * 3L), nrow = n)
  colnames(X) <- paste0("x", 1:3)

  eta <- 0.8 * X[, 1] - 0.6 * X[, 2] + 0.25 * X[, 3]
  Te <- rexp(n, rate = exp(eta))
  Tc <- rexp(n, rate = 0.35)

  time <- pmax(1, round(pmin(Te, Tc) * 5))
  status <- as.integer(Te <= Tc)

  dat <- data.frame(
    time = time,
    status = status,
    X
  )

  fit <- survival::coxph(
    survival::Surv(time, status) ~ x1 + x2 + x3,
    data = dat,
    ties = "efron",
    x = TRUE,
    model = TRUE
  )

  got <- .selsurv_breslow_basehaz(
    fit,
    centered = FALSE
  )

  sf1 <- survival::survfit(
    fit,
    se.fit = FALSE,
    ctype = 1
  )

  beta <- stats::coef(fit)
  beta <- ifelse(is.na(beta), 0, beta)
  offset <- sum(fit$means * beta)

  expected <- as.numeric(sf1$cumhaz) * exp(-offset)

  expect_equal(
    got$hazard,
    expected,
    tolerance = 1e-12
  )

  default_bh <- survival::basehaz(
    fit,
    centered = FALSE
  )

  expect_gt(
    max(abs(default_bh$hazard - got$hazard)),
    1e-6
  )
})


test_that("clean re-selection refits default to Efron and legacy forces Breslow", {
  expect_identical(
    eval(formals(fit_unpenalized_cox)$ties),
    "efron"
  )

  expect_identical(
    eval(formals(fit_original_workflow)$ties),
    "efron"
  )

  expect_identical(
    eval(formals(run_reselection_bootstrap)$ties),
    "efron"
  )

  expect_identical(
    eval(formals(coxlasso_boot_reselect)$ties),
    "efron"
  )

  expect_identical(
    eval(formals(fit_support_population_target)$ties),
    "efron"
  )

  expect_identical(
    .reselection_resolve_refit_ties(
      "efron",
      "independent_streams"
    ),
    "efron"
  )

  expect_identical(
    .reselection_resolve_refit_ties(
      "efron",
      "legacy_nki70"
    ),
    "breslow"
  )
})


test_that("re-selection Cox-Lasso calls use explicit CV conventions", {
  args <- .reselection_prepare_cox_glmnet_args(
    list(),
    cv = TRUE
  )

  expect_identical(args$type.measure, "deviance")
  expect_identical(args$grouped, TRUE)

  supports_cox_ties <- "cox.ties" %in%
    names(formals(glmnet::glmnet))

  if (supports_cox_ties) {
    expect_identical(args$cox.ties, "breslow")
  } else {
    expect_false("cox.ties" %in% names(args))
  }
})


test_that("reference re-selection generator uses Breslow baseline from Efron refit", {
  set.seed(90153)

  n <- 150L
  X <- matrix(rnorm(n * 3L), nrow = n)
  colnames(X) <- paste0("x", 1:3)

  eta <- 0.9 * X[, 1] - 0.7 * X[, 2]
  Te <- rexp(n, rate = exp(eta))
  Tc <- rexp(n, rate = 0.30)

  time <- pmax(1, round(pmin(Te, Tc) * 5))
  status <- as.integer(Te <= Tc)

  ref <- fit_unpenalized_cox(
    X,
    time,
    status,
    support = 1:3
  )

  expect_true(ref$success)
  expect_identical(ref$fit$method, "efron")

  gen <- fit_modelbased_generator(
    X,
    time,
    status,
    original_refit = ref,
    generator_mode = "reference"
  )

  expected <- .selsurv_breslow_basehaz(
    ref$fit,
    centered = FALSE
  )

  expect_equal(
    gen$basehaz$time,
    expected$time
  )

  expect_equal(
    gen$basehaz$hazard,
    expected$hazard,
    tolerance = 1e-12
  )
})


test_that("exact Cox refits remain available outside legacy mode", {
  expect_identical(
    .reselection_resolve_refit_ties(
      "exact",
      "independent_streams"
    ),
    "exact"
  )

  expect_identical(
    .reselection_resolve_refit_ties(
      "exact",
      "legacy_nki70"
    ),
    "breslow"
  )
})


test_that("parametric bootstrap defaults old fit objects to Breslow selection", {
  fake_fit <- list(selection_ties = NULL)

  selection_ties <- if (!is.null(fake_fit$selection_ties)) {
    fake_fit$selection_ties
  } else {
    "breslow"
  }

  expect_identical(selection_ties, "breslow")
})


test_that("re-selection package defaults use general-purpose policy", {
  expect_identical(
    eval(formals(coxlasso_boot_reselect)$B),
    500L
  )
  expect_identical(
    eval(formals(coxlasso_boot_reselect)$conf_level),
    0.95
  )
  expect_null(
    eval(formals(coxlasso_boot_reselect)$seed)
  )
  expect_identical(
    eval(formals(coxlasso_boot_reselect)$standardize),
    TRUE
  )

  expect_identical(
    eval(formals(run_reselection_bootstrap)$B),
    500L
  )
  expect_identical(
    eval(formals(run_reselection_bootstrap)$conf_level),
    0.95
  )
  expect_null(
    eval(formals(run_reselection_bootstrap)$seed)
  )
  expect_identical(
    eval(formals(run_reselection_bootstrap)$standardize),
    TRUE
  )
})


test_that("NULL re-selection seed generates a reproducible root seed", {
  set.seed(90156)

  s1 <- .reselection_resolve_root_seed(
    NULL,
    "independent_streams"
  )

  set.seed(90156)

  s2 <- .reselection_resolve_root_seed(
    NULL,
    "independent_streams"
  )

  expect_identical(s1, s2)
  expect_type(s1, "integer")
  expect_gte(s1, 0L)
  expect_lte(s1, .Machine$integer.max)

  expect_identical(
    .reselection_resolve_root_seed(
      2026L,
      "independent_streams"
    ),
    2026L
  )

  expect_error(
    .reselection_resolve_root_seed(
      NULL,
      "legacy_nki70"
    ),
    "requires an explicit seed"
  )
})


test_that("split-sample defaults separate interval and penalty alpha", {
  expect_identical(
    eval(formals(coxlasso_split)$conf_level),
    0.95
  )
  expect_identical(
    eval(formals(coxlasso_split)$penalty_alpha),
    1
  )
  expect_false("alpha" %in% names(formals(coxlasso_split)))
})


test_that("historical direct bootstrap defaults remain unchanged", {
  expect_identical(
    eval(formals(coxlasso_boot_parametric)$B),
    200
  )
  expect_identical(
    eval(formals(coxlasso_boot_parametric)$alpha),
    0.1
  )

  expect_identical(
    eval(formals(coxlasso_wildboot)$B),
    500
  )
  expect_identical(
    eval(formals(coxlasso_wildboot)$alpha),
    0.1
  )

  expect_identical(
    eval(formals(coxlasso_efron_events)$B),
    500
  )
  expect_identical(
    eval(formals(coxlasso_efron_events)$alpha),
    0.1
  )

  expect_identical(
    eval(formals(coxlasso_debias)$ridge),
    1e-4
  )
})
