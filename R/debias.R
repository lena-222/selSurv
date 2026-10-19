cox_score_info <- function(X, time, status, beta, ties = c("breslow", "efron")) {
  # Work on full risk sets at each distinct event time. Censoring at an
  # event time therefore remains at risk at that time.
  ties <- match.arg(ties)
  n <- nrow(X)
  p <- ncol(X)
  eta <- as.numeric(X %*% beta)
  if (any(!is.finite(eta))) {
    stop("Non-finite linear predictors in Cox score calculation.", call. = FALSE)
  }

  score <- numeric(p)
  info <- matrix(0, p, p)
  contributions <- matrix(0, n, p)
  event_times <- sort(unique(time[status == 1L]))

  for (tt in event_times) {
    risk <- which(time >= tt)
    failures <- which(time == tt & status == 1L)
    d <- length(failures)
    XR <- X[risk, , drop = FALSE]
    XD <- X[failures, , drop = FALSE]

    # A time-specific log-sum-exp shift avoids clipping the linear predictor.
    # The common shift cancels from all ratios of risk-set sums.
    shift <- max(eta[risk])
    wR <- exp(eta[risk] - shift)
    wD <- exp(eta[failures] - shift)
    S0 <- sum(wR)
    S1 <- colSums(sweep(XR, 1L, wR, `*`))
    S2 <- crossprod(sweep(XR, 1L, sqrt(wR), `*`))
    E0 <- sum(wD)
    E1 <- colSums(sweep(XD, 1L, wD, `*`))
    E2 <- crossprod(sweep(XD, 1L, sqrt(wD), `*`))

    adjustment <- numeric(p)
    for (k in seq_len(d) - 1L) {
      fraction <- if (ties == "efron") k / d else 0
      denom <- S0 - fraction * E0
      if (!is.finite(denom) || denom <= 0) {
        stop("Invalid denominator in Cox score calculation.", call. = FALSE)
      }
      mean_risk <- (S1 - fraction * E1) / denom
      adjustment <- adjustment + mean_risk
      info <- info + (S2 - fraction * E2) / denom -
        tcrossprod(mean_risk)
    }

    # Allocate the average risk adjustment symmetrically among tied failures.
    # This guarantees permutation invariance of the score contributions.
    contributions[failures, ] <- sweep(XD, 2L, adjustment / d, `-`)
    score <- score + colSums(XD) - adjustment
  }
  list(score = score / n, info = info / n,
       contributions = contributions, n_events = sum(status))
}

#' Approximate debiasing for a Cox-Lasso fit
#'
#' Applies a one-step Cox score correction using a ridge-stabilized
#' inverse information matrix. This research implementation does not
#' supply a general high-dimensional validity guarantee. It requires
#' a Lasso pilot (`fit$alpha = 1`); other penalties are rejected.
#'
#' The score calculation warrants additional scrutiny for tied
#' event times. Results from such data should be treated cautiously.
#'
#' @param fit A `selSurv_fit` object.
#' @param ridge Ridge-stabilization term added to the
#'   information matrix before inversion.
#' @return A `selSurv_debias` object.
#' @seealso \code{\link{selSurv_workflows}}
#' @export
coxlasso_debias <- function(fit, ridge = 1e-4) {
  .selsurv_require_lasso_fit(fit, "coxlasso_debias()")

  beta <- fit$beta_lasso
  ties <- if (is.null(fit$selection_ties)) "breslow" else fit$selection_ties
  si <- cox_score_info(fit$X, fit$time, fit$status, beta, ties = ties)
  p <- length(beta)
  info_ridge <- si$info + diag(ridge, p)
  theta <- tryCatch(solve(info_ridge), error = function(e) MASS::ginv(info_ridge))
  correction <- as.numeric(theta %*% si$score)
  beta_debiased <- beta + correction
  names(beta_debiased) <- names(beta)
  vcov <- theta %*% si$info %*% t(theta) / nrow(fit$X)
  se <- sqrt(pmax(diag(vcov), 0))
  names(se) <- names(beta)
  out <- list(
    fit = fit,
    beta_lasso = beta,
    beta_debiased = beta_debiased,
    se = se,
    theta = theta,
    score = si$score,
    info = si$info,
    score_contributions = si$contributions,
    ridge = ridge
  )
  class(out) <- "selSurv_debias"
  out
}

#' @export
print.selSurv_debias <- function(x, ...) {
  cat("selSurv debiased Cox-Lasso object\n")
  cat("p =", length(x$beta_debiased), " ridge =", x$ridge, "\n")
  invisible(x)
}
