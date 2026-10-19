testthat::test_that(
  "Cox score and information agree with an independent coxph reference", {

  score_info <- getFromNamespace("cox_score_info", "selSurv")

  X <- cbind(
    x1 = c(0, 1, 2, 3, 0.3, -0.5, 1.5, 2.2),
    x2 = c(1, -1, 0, 2, -0.2, 1.3, -0.7, 0.5)
  )
  status <- c(1L, 1L, 0L, 1L, 0L, 1L, 0L, 1L)
  beta <- c(0.3, -0.1)
  E <- diag(2)
  n <- nrow(X)

  scenarios <- list(
    untied = seq_len(n),
    tied = c(1, 1, 1, 2, 2, 3, 4, 5)
  )

  for (time in scenarios) {
    dat <- data.frame(
      time = time,
      status = status,
      x1 = X[, 1],
      x2 = X[, 2]
    )

    for (ties in c("breslow", "efron")) {
      ll <- function(b) {
        fit <- suppressWarnings(survival::coxph(
          survival::Surv(time, status) ~ x1 + x2,
          data = dat,
          ties = ties,
          init = as.numeric(b),
          control = survival::coxph.control(iter.max = 0)
        ))
        as.numeric(tail(fit$loglik, 1L))
      }

      h1 <- 1e-5
      h2 <- 2e-4

      score_ref <- vapply(seq_len(2), function(j) {
        (ll(beta + h1 * E[, j]) -
           ll(beta - h1 * E[, j])) / (2 * h1 * n)
      }, numeric(1))

      info_ref <- matrix(0, 2, 2)
      f0 <- ll(beta)

      for (j in seq_len(2)) {
        info_ref[j, j] <- -(
          ll(beta + h2 * E[, j]) - 2 * f0 +
            ll(beta - h2 * E[, j])
        ) / (h2^2 * n)
      }

      e1 <- E[, 1]
      e2 <- E[, 2]

      info_ref[1, 2] <- -(
        ll(beta + h2 * (e1 + e2)) -
          ll(beta + h2 * (e1 - e2)) -
          ll(beta + h2 * (-e1 + e2)) +
          ll(beta - h2 * (e1 + e2))
      ) / (4 * h2^2 * n)

      info_ref[2, 1] <- info_ref[1, 2]

      actual <- score_info(X, time, status, beta, ties = ties)

      testthat::expect_lt(
        max(abs(as.numeric(actual$score) - score_ref)),
        1e-6
      )
      testthat::expect_lt(
        max(abs(actual$info - info_ref)),
        2e-5
      )
    }
  }
})
