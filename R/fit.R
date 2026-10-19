#' Fit a Cox-Lasso workflow
#'
#' Fits a penalized Cox model with cross-validated tuning and stores
#' the selected model. The default mixing parameter (`alpha = 1`)
#' yields the Lasso.
#'
#' @usage
#' coxlasso_fit(
#'   X,
#'   time,
#'   status,
#'   lambda_choice = c("min", "1se"),
#'   nfolds = 10,
#'   alpha = 1,
#'   standardize = TRUE,
#'   ...
#' )
#' @param X Numeric matrix or data frame of covariates. In data frames,
#'   numeric and integer columns remain numeric, while logical, factor,
#'   ordered-factor, and character columns are expanded using treatment
#'   coding with the first level as the reference. Ordered factors are treated
#'   as ordinary categorical predictors. Categorical predictors contribute
#'   \eqn{k-1} separately penalized design-matrix coordinates; they are not
#'   selected as groups. To control the reference category, supply an explicit
#'   factor level order or use `relevel()` before fitting.
#' @param time Follow-up time.
#' @param status Event indicator coded 0/1.
#' @param lambda_choice Either `"min"` or `"1se"`.
#' @param nfolds Number of cross-validation folds.
#' @param alpha Elastic-net mixing parameter passed to `glmnet`.
#'   `alpha = 1` gives the Lasso. The downstream inference
#'   workflows have been evaluated primarily for this setting.
#' @param standardize Passed to `glmnet`.
#' @param ... Additional arguments passed to `glmnet::cv.glmnet()`.
#' @return A `selSurv_fit` object.
#' @examples
#' set.seed(2026)
#' n <- 120L
#' X <- matrix(stats::rnorm(n * 4L), nrow = n)
#' colnames(X) <- paste0("x", seq_len(ncol(X)))
#' event_time <- stats::rexp(n, rate = 0.12 * exp(X[, 1] - 0.7 * X[, 2]))
#' censor_time <- stats::rexp(n, rate = 0.06)
#' time <- pmin(event_time, censor_time)
#' status <- as.integer(event_time <= censor_time)
#' fit <- coxlasso_fit(X, time, status, nfolds = 3, lambda_choice = "min")
#' fit$selected
#' @export
coxlasso_fit <- function(X,
                         time,
                         status,
                         lambda_choice = c("min", "1se"),
                         nfolds = 10,
                         alpha = 1,
                         standardize = TRUE,
                         ...) {
  lambda_choice <- match.arg(lambda_choice)
  inp <- check_survival_inputs(X, time, status)
  X <- inp$X
  y <- survival::Surv(inp$time, inp$status)
  dots <- list(...)

  selection_ties <- if ("cox.ties" %in% names(dots)) {
    match.arg(
      as.character(dots[["cox.ties"]]),
      c("breslow", "efron")
    )
  } else {
    "breslow"
  }

  cv_args <- c(
    list(
      x = X,
      y = y,
      family = "cox",
      alpha = alpha,
      nfolds = nfolds,
      standardize = standardize
    ),
    dots
  )

  cv_args <- .selsurv_prepare_cox_glmnet_args(
    cv_args,
    selection_ties = selection_ties,
    cv = TRUE
  )

  cvfit <- do.call(glmnet::cv.glmnet, cv_args)
  lambda <- if (lambda_choice == "min") cvfit$lambda.min else cvfit$lambda.1se
  beta <- as.numeric(stats::coef(cvfit, s = lambda))
  names(beta) <- colnames(X)
  support <- selected_from_beta(beta)
  refit <- fit_cox_support(
    X,
    inp$time,
    inp$status,
    support,
    ties = "efron"
  )
  theta_refit <- zero_pad(refit$coef, support, ncol(X), colnames(X))
  se_refit <- zero_pad(refit$se, support, ncol(X), colnames(X))
  out <- list(
    X = X,
    time = inp$time,
    status = inp$status,
    y = y,
    cvfit = cvfit,
    lambda = lambda,
    lambda_choice = lambda_choice,
    alpha = alpha,
    standardize = standardize,
    cv_type_measure = cv_args[["type.measure"]],
    cv_grouped = cv_args[["grouped"]],
    selection_ties = selection_ties,
    refit_ties = "efron",
    beta_lasso = beta,
    support = support,
    selected = colnames(X)[support],
    refit = refit,
    theta_refit = theta_refit,
    coef = theta_refit,
    se = se_refit,
    call = match.call()
  )
  class(out) <- "selSurv_fit"
  out
}

#' Refit Cox model on selected support
#'
#' @param fit A `selSurv_fit` object.
#' @param support Optional integer support. Defaults to `fit$support`.
#' @param X Optional numeric matrix or data frame. After input normalization,
#'   its design-matrix columns must match `fit$X` exactly in names and order.
#' @param time Optional time vector for bootstrap data.
#' @param status Optional status vector for bootstrap data.
#' @return A list with support-specific and zero-padded estimates.
#' @export
coxlasso_refit <- function(fit, support = fit$support, X = fit$X, time = fit$time, status = fit$status) {
  inp <- check_survival_inputs(X, time, status)

  if (is.null(fit$X) ||
      !identical(colnames(inp$X), colnames(fit$X))) {
    stop(
      "X must expand to exactly the same design-matrix columns, in the same ",
      "order, as fit$X.",
      call. = FALSE
    )
  }

  refit <- fit_cox_support(inp$X, inp$time, inp$status, support)
  theta <- zero_pad(refit$coef, support, ncol(inp$X), colnames(inp$X))
  se <- zero_pad(refit$se, support, ncol(inp$X), colnames(inp$X))
  list(
    support = support,
    selected = colnames(inp$X)[support],
    refit = refit,
    theta = theta,
    coef = theta,
    se = se,
    coef_support = refit$coef,
    se_support = refit$se
  )
}

#' @export
print.selSurv_fit <- function(x, ...) {
  cat("selSurv Cox-Lasso fit\n")
  cat("n =", nrow(x$X), " p =", ncol(x$X), " events =", sum(x$status), "\n")
  cat("lambda_choice =", x$lambda_choice, " lambda =", signif(x$lambda, 4), "\n")
  cat("selected variables:", if (length(x$selected)) paste(x$selected, collapse = ", ") else "none", "\n")
  invisible(x)
}
