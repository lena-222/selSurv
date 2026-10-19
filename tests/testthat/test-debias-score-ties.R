# Numerical regression tests for the corrected Cox partial-likelihood kernel.
# These tests are placed in an isolated 9017 source snapshot.
testthat::test_that("Cox score and observed information handle tied risk sets", {
  f <- getFromNamespace("cox_score_info", "selSurv")
  X <- cbind(x1 = c(0, 1, 2, 3), x2 = c(1, -1, 0, 2))
  time <- c(1, 1, 1, 2)
  status <- c(1L, 1L, 0L, 1L)
  beta <- c(0.2, -0.1)
  ord <- c(3L, 1L, 4L, 2L)
  for (method in c("breslow", "efron")) {
    actual <- f(X, time, status, beta, ties = method)
    reordered <- f(X[ord, , drop = FALSE], time[ord],
                   status[ord], beta, ties = method)
    testthat::expect_equal(actual$score, reordered$score, tolerance = 1e-12)
    testthat::expect_equal(actual$info, reordered$info, tolerance = 1e-12)
    testthat::expect_equal(actual$contributions[ord, , drop = FALSE],
                           reordered$contributions, tolerance = 1e-12)
    testthat::expect_equal(unname(colSums(actual$contributions) / nrow(X)),
                           unname(actual$score), tolerance = 1e-12)
    testthat::expect_equal(actual$info, t(actual$info), tolerance = 1e-12)
    testthat::expect_true(all(is.finite(actual$info)))
  }
  testthat::expect_gt(max(abs(f(X, time, status, beta, "breslow")$score -
                               f(X, time, status, beta, "efron")$score)), 1e-6)
})

testthat::test_that("Debiasing respects penalized selection ties and older fits", {
  f <- getFromNamespace("cox_score_info", "selSurv")
  X <- cbind(x1 = c(0, 1, 2), x2 = c(-1, 0, 1))
  time <- c(1, 1, 1)
  status <- c(1L, 1L, 0L)
  beta <- c(x1 = 0.1, x2 = -0.2)
  fit <- list(X = X, time = time, status = status,
              beta_lasso = beta, alpha = 1,
              selection_ties = "efron", refit_ties = "breslow")
  db <- selSurv::coxlasso_debias(fit)
  testthat::expect_equal(unname(db$score),
                         unname(f(X, time, status, beta, "efron")$score),
                         tolerance = 1e-12)
  fit$selection_ties <- NULL
  db_old <- selSurv::coxlasso_debias(fit)
  testthat::expect_equal(unname(db_old$score),
                         unname(f(X, time, status, beta, "breslow")$score),
                         tolerance = 1e-12)
})
