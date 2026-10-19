#' Split-sample inference after Cox-Lasso selection
#'
#' Uses one subsample to choose predictors and a disjoint
#' subsample to estimate coefficients and form Wald intervals.
#'
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
#' @param split_fraction Fraction of observations used for selection.
#' @param seed Optional random seed.
#' @param lambda_choice Cox-Lasso tuning rule used in the selection split.
#' @param nfolds Number of folds for the selection split.
#' @param conf_level Confidence level for the split-sample Wald intervals.
#'   Must be strictly between zero and one.
#' @param penalty_alpha Elastic-net mixing parameter used for penalized Cox
#'   screening in the selection split. Passed as `alpha` to `coxlasso_fit()`.
#'   Values other than 1 are computationally supported but trigger a warning
#'   because other penalties have not been assessed separately
#'   in the inference workflows.
#' @param ... Additional arguments passed to `coxlasso_fit()`. For backward
#'   compatibility, a named `alpha` argument is interpreted as the deprecated
#'   confidence-interval significance level; use `conf_level` instead.
#' @return A `selSurv_split` object.
#' @export
coxlasso_split <- function(
    X,
    time,
    status,
    split_fraction = 0.5,
    seed = NULL,
    lambda_choice = c("min", "1se"),
    nfolds = 10,
    conf_level = 0.95,
    penalty_alpha = 1,
    ...
) {
  lambda_choice <- match.arg(lambda_choice)

  conf_level_supplied <- !missing(conf_level)
  dots <- list(...)

  # Before 0.1.0.9015, `alpha` denoted the CI significance level in this
  # function, which prevented callers from controlling the Cox-Lasso
  # elastic-net mixing parameter. Preserve that old interpretation with a
  # warning rather than silently changing the meaning of existing calls.
  if ("alpha" %in% names(dots)) {
    legacy_alpha <- dots[["alpha"]]

    if (isTRUE(conf_level_supplied)) {
      stop(
        "Do not supply both `conf_level` and the deprecated split-sample ",
        "`alpha` argument.",
        call. = FALSE
      )
    }

    if (length(legacy_alpha) != 1L ||
        is.na(legacy_alpha) ||
        !is.numeric(legacy_alpha) ||
        !is.finite(legacy_alpha) ||
        legacy_alpha <= 0 ||
        legacy_alpha >= 1) {
      stop(
        "Deprecated split-sample `alpha` must be strictly between 0 and 1.",
        call. = FALSE
      )
    }

    warning(
      "`alpha` as a confidence-interval significance level in ",
      "`coxlasso_split()` is deprecated. Use `conf_level` instead. ",
      "Use `penalty_alpha` for the Cox-Lasso elastic-net mixing parameter.",
      call. = FALSE
    )

    conf_level <- 1 - legacy_alpha
    dots[["alpha"]] <- NULL
  }

  if (length(conf_level) != 1L ||
      is.na(conf_level) ||
      !is.numeric(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1) {
    stop(
      "`conf_level` must be one number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  if (length(penalty_alpha) != 1L ||
      is.na(penalty_alpha) ||
      !is.numeric(penalty_alpha) ||
      !is.finite(penalty_alpha) ||
      penalty_alpha < 0 ||
      penalty_alpha > 1) {
    stop(
      "`penalty_alpha` must be one number between 0 and 1.",
      call. = FALSE
    )
  }

  .selsurv_warn_nonlasso(
    penalty_alpha,
    "coxlasso_split()",
    "penalty_alpha"
  )

  if (!is.null(seed)) set.seed(seed)
  inp <- check_survival_inputs(X, time, status)
  n <- nrow(inp$X)
  idx_train <- sort(sample(seq_len(n), size = floor(split_fraction * n)))
  idx_test <- setdiff(seq_len(n), idx_train)
  fit_args <- c(
    list(
      X = inp$X[idx_train, , drop = FALSE],
      time = inp$time[idx_train],
      status = inp$status[idx_train],
      lambda_choice = lambda_choice,
      nfolds = nfolds,
      alpha = penalty_alpha
    ),
    dots
  )

  fit_sel <- do.call(
    coxlasso_fit,
    fit_args
  )
  support <- fit_sel$support
  ref <- fit_cox_support(
    inp$X[idx_test, , drop = FALSE],
    inp$time[idx_test],
    inp$status[idx_test],
    support,
    ties = "efron"
  )
  theta <- zero_pad(ref$coef, support, ncol(inp$X), colnames(inp$X))
  ci <- data.frame(term = colnames(inp$X), estimate = theta, lower = NA_real_, upper = NA_real_, row.names = NULL)
  if (length(support) > 0 && !is.null(ref$fit)) {
    cf <- stats::coef(ref$fit)
    se <- sqrt(diag(stats::vcov(ref$fit)))
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    ci$lower[support] <- cf - z * se
    ci$upper[support] <- cf + z * se
  }
  out <- list(
    fit_selection = fit_sel,
    support = support,
    selected = colnames(inp$X)[support],
    refit = ref,
    theta = theta,
    ci = ci,
    train_index = idx_train,
    test_index = idx_test,
    conf_level = conf_level,
    alpha = 1 - conf_level,
    penalty_alpha = penalty_alpha
  )
  class(out) <- "selSurv_split"
  out
}

#' @export
print.selSurv_split <- function(x, ...) {
  cat("selSurv split-sample object\n")
  cat("selected variables:", if (length(x$selected)) paste(x$selected, collapse = ", ") else "none", "\n")
  print(utils::head(x$ci))
  invisible(x)
}
