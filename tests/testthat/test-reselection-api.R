test_that("re-selection API returns a stable object on a small synthetic example", {
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")

  set.seed(1401)
  n <- 80L
  p <- 5L
  x <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("x", seq_len(p))
  lp <- drop(x %*% c(0.7, -0.6, 0.3, 0, 0))
  t_event <- rexp(n, rate = 0.08 * exp(lp))
  t_cens <- rexp(n, rate = 0.05)
  time <- pmin(t_event, t_cens)
  status <- as.integer(t_event <= t_cens)

  original <- fit_original_fixed_lambda(
    x, time, status,
    lambda = 0.05,
    standardize = FALSE,
    ties = "breslow"
  )

  fit <- coxlasso_boot_reselect(
    x, time, status,
    B = 2L,
    methods = c("fixed_support", "reselect_fixed_lambda"),
    seed = 11L,
    nfolds = 5L,
    original = original,
    standardize = FALSE,
    ties = "breslow",
    rng_mode = "independent_streams",
    progress = FALSE
  )

  expect_s3_class(fit, "selsurv_reselection")
  expect_s3_class(fit, "cox_reselection_bootstrap")

  d <- method_level_diagnostics(fit)
  expect_true(all(c("fixed_support", "reselect_fixed_lambda") %in% d$method))
  expect_true(all(d$selection_success_rate >= 0 & d$selection_success_rate <= 1))

  s <- summary(fit)
  expect_s3_class(s, "summary.selsurv_reselection")
  expect_true(all(c("intervals", "diagnostics", "settings") %in% names(s)))
})

test_that("re-selection aliases normalize to canonical method names", {
  expect_equal(
    .normalize_reselection_methods(c("reselect_fixed", "reselect_cvmin")),
    c("reselect_fixed_lambda", "reselect_cv_min")
  )
})
