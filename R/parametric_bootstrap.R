generate_cox_bootstrap <- function(fit, generator_refit = fit$refit) {
  X <- fit$X
  n <- nrow(X)
  support <- generator_refit$support
  if (length(support) == 0 || is.null(generator_refit$fit)) {
    base_fit <- survival::coxph(survival::Surv(fit$time, fit$status) ~ 1)
    bh <- .selsurv_breslow_basehaz(
      base_fit,
      centered = FALSE
    )
    lp <- rep(0, n)
  } else {
    bh <- .selsurv_breslow_basehaz(
      generator_refit$fit,
      centered = FALSE
    )
    beta <- generator_refit$coef
    lp <- as.numeric(X[, support, drop = FALSE] %*% beta)
  }
  target <- -log(stats::runif(n)) / exp(pmin(lp, 30))
  event_time <- inverse_baseline_time(bh, target)
  censor_sampler <- km_censor_sampler(fit$time, fit$status)
  censor_time <- censor_sampler(n)
  time_star <- pmax(pmin(event_time, censor_time), 1e-8)
  status_star <- as.integer(event_time <= censor_time)
  list(X = X, time = time_star, status = status_star)
}

fit_bootstrap_workflow <- function(fit, boot_data, workflow, lambda_choice, nfolds, alpha_glm, standardize) {
  Xb <- boot_data$X

  selection_ties <- if (!is.null(fit$selection_ties)) {
    fit$selection_ties
  } else {
    "breslow"
  }
  tb <- boot_data$time
  sb <- boot_data$status
  p <- ncol(Xb)
  if (workflow == "fixed_support") {
    ref <- fit_cox_support(
      Xb, tb, sb, fit$support,
      ties = "efron"
    )
    theta <- zero_pad(ref$coef, fit$support, p, colnames(Xb))
    return(list(theta = theta, support = fit$support, selected = colnames(Xb)[fit$support], ok = ref$converged))
  }
  if (workflow == "reselect_fixed_lambda") {
    glmnet_args <- list(
      x = Xb,
      y = survival::Surv(tb, sb),
      family = "cox",
      alpha = alpha_glm,
      standardize = standardize,
      lambda = fit$lambda
    )
    glmnet_args <- .selsurv_prepare_cox_glmnet_args(
      glmnet_args,
      selection_ties = selection_ties,
      cv = FALSE
    )
    gf <- tryCatch(
      do.call(glmnet::glmnet, glmnet_args),
      error = function(e) e
    )
    if (inherits(gf, "error")) return(list(theta = rep(NA_real_, p), support = integer(0), selected = character(0), ok = FALSE))
    beta <- as.numeric(stats::coef(gf, s = fit$lambda))
    support <- selected_from_beta(beta)
    ref <- fit_cox_support(
      Xb, tb, sb, support,
      ties = "efron"
    )
    theta <- zero_pad(ref$coef, support, p, colnames(Xb))
    return(list(theta = theta, support = support, selected = colnames(Xb)[support], ok = ref$converged))
  }
  fb <- tryCatch(
    coxlasso_fit(
      Xb, tb, sb,
      lambda_choice = lambda_choice,
      nfolds = nfolds,
      alpha = alpha_glm,
      standardize = standardize,
      cox.ties = selection_ties
    ),
    error = function(e) e
  )
  if (inherits(fb, "error")) return(list(theta = rep(NA_real_, p), support = integer(0), selected = character(0), ok = FALSE))
  list(theta = fb$theta_refit, support = fb$support, selected = fb$selected, ok = TRUE, fit = fb)
}

#' Model-based bootstrap for Cox-Lasso post-selection workflows
#'
#' Runs a model-based Cox bootstrap, either keeping the observed
#' support fixed or repeating selection within each resample. Under
#' re-selection, unselected coefficients are stored as zero. The
#' resulting distribution describes the fitted procedure rather
#' than an exact selection-event-conditional distribution.
#'
#' @usage
#' coxlasso_boot_parametric(
#'   fit,
#'   B = 200,
#'   workflow = c(
#'     "fixed_support", "reselect_fixed_lambda", "reselect_retuned"
#'   ),
#'   lambda_choice = fit$lambda_choice,
#'   nfolds = 10,
#'   alpha = 0.1,
#'   ci_type = c("percentile", "basic"),
#'   seed = NULL
#' )
#' @param fit A `selSurv_fit` object.
#' @param B Number of bootstrap samples.
#' @param workflow Bootstrap workflow: fixed support, fixed-penalty re-selection,
#'   or retuned re-selection.
#' @param lambda_choice Tuning rule for `"reselect_retuned"`.
#' @param nfolds Number of folds for retuning.
#' @param alpha Interval significance level, so `1 - alpha` is
#'   the nominal confidence level. This is *not* the mixing parameter.
#' @param ci_type `"percentile"` or `"basic"`.
#' @param seed Optional random seed.
#' @return A `selSurv_boot` object.
#' @section Choosing between bootstrap interfaces:
#' For fixed-support bootstrap, use [coxlasso_boot()].
#' Set `inference = "conditional"`.
#'
#' To repeat selection, use [coxlasso_boot_reselect()].
#' This interface has its own tuning and return conventions.
#' Fits with a non-Lasso mixing parameter can be processed, but
#' a warning is issued because these settings have not been
#' assessed separately for inference. Here `alpha` is the
#' significance level, not the penalty mixing parameter.
#' @seealso \code{\link{selSurv_workflows}}
#' @export
coxlasso_boot_parametric <- function(fit, B = 200, workflow = c("fixed_support", "reselect_fixed_lambda", "reselect_retuned"), lambda_choice = fit$lambda_choice, nfolds = 10, alpha = 0.1, ci_type = c("percentile", "basic"), seed = NULL) {
  workflow <- match.arg(workflow)
  ci_type <- match.arg(ci_type)

  .selsurv_warn_nonlasso(
    fit$alpha,
    "coxlasso_boot_parametric()",
    "fit$alpha"
  )

  if (!is.null(seed)) set.seed(seed)
  p <- ncol(fit$X)
  theta_star <- matrix(NA_real_, B, p)
  colnames(theta_star) <- colnames(fit$X)
  selected_list <- vector("list", B)
  ok <- logical(B)
  for (b in seq_len(B)) {
    bd <- generate_cox_bootstrap(fit)
    fb <- fit_bootstrap_workflow(fit, bd, workflow = workflow, lambda_choice = lambda_choice, nfolds = nfolds, alpha_glm = fit$alpha, standardize = fit$standardize)
    theta_star[b, ] <- fb$theta
    selected_list[[b]] <- fb$selected
    ok[b] <- isTRUE(fb$ok)
  }
  theta_used <- theta_star[ok, , drop = FALSE]
  ci <- make_interval(fit$theta_refit, theta_used, alpha = alpha, type = ci_type)
  reselect_freq <- colMeans(theta_star != 0, na.rm = TRUE)
  zero_freq <- colMeans(theta_star == 0, na.rm = TRUE)
  diagnostics <- data.frame(term = colnames(fit$X), reselect_frequency = reselect_freq, zero_frequency = zero_freq, row.names = NULL)
  out <- list(
    method = paste0("model_based_", workflow),
    workflow = workflow,
    theta_hat = fit$theta_refit,
    theta_star = theta_star,
    ci = ci,
    diagnostics = diagnostics,
    selected = selected_list,
    ok = ok,
    success_rate = mean(ok),
    alpha = alpha,
    B = B,
    ci_type = ci_type
  )
  class(out) <- "selSurv_boot"
  out
}
