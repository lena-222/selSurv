#' Multiplier wild bootstrap for debiased Cox-Lasso estimator
#'
#' @param db A `selSurv_debias` object.
#' @param B Number of bootstrap draws.
#' @param alpha Significance level.
#' @param type Interval type, currently percentile on the bootstrap debiased scale.
#' @param seed Optional random seed.
#' @return A `selSurv_boot` object.
#' @details The debiased estimator must originate from a Cox-Lasso fit with
#'   penalty mixing parameter `alpha = 1`.
#' @export
coxlasso_wildboot <- function(db, B = 500, alpha = 0.1, type = c("percentile"), seed = NULL) {
  .selsurv_require_lasso_fit(db$fit, "coxlasso_wildboot()")

  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(db$score_contributions)
  p <- ncol(db$score_contributions)
  theta_star <- matrix(NA_real_, B, p)
  colnames(theta_star) <- names(db$beta_debiased)
  for (b in seq_len(B)) {
    g <- stats::rnorm(n)
    boot_score <- colSums(db$score_contributions * g) / n
    theta_star[b, ] <- db$beta_debiased + as.numeric(db$theta %*% boot_score)
  }
  ci <- make_interval(db$beta_debiased, theta_star, alpha = alpha, type = "percentile")
  out <- list(method = "wild_score", theta_hat = db$beta_debiased, theta_star = theta_star, ci = ci, alpha = alpha, B = B)
  class(out) <- "selSurv_boot"
  out
}

#' Event-level bootstrap for debiased Cox-Lasso score contributions
#'
#' @param db A `selSurv_debias` object.
#' @param B Number of bootstrap draws.
#' @param alpha Significance level.
#' @param seed Optional random seed.
#' @return A `selSurv_boot` object.
#' @details The debiased estimator must originate from a Cox-Lasso fit with
#'   penalty mixing parameter `alpha = 1`.
#' @export
coxlasso_efron_events <- function(db, B = 500, alpha = 0.1, seed = NULL) {
  .selsurv_require_lasso_fit(db$fit, "coxlasso_efron_events()")

  if (!is.null(seed)) set.seed(seed)
  event_rows <- which(db$fit$status == 1L)
  if (!length(event_rows)) stop("No events available for event-level bootstrap.", call. = FALSE)
  contrib <- db$score_contributions[event_rows, , drop = FALSE]
  p <- ncol(contrib)
  n <- nrow(db$fit$X)
  theta_star <- matrix(NA_real_, B, p)
  colnames(theta_star) <- names(db$beta_debiased)
  centered <- scale(contrib, center = TRUE, scale = FALSE)
  for (b in seq_len(B)) {
    idx <- sample(seq_len(nrow(centered)), nrow(centered), replace = TRUE)
    boot_score <- colSums(centered[idx, , drop = FALSE]) / n
    theta_star[b, ] <- db$beta_debiased + as.numeric(db$theta %*% boot_score)
  }
  ci <- make_interval(db$beta_debiased, theta_star, alpha = alpha, type = "percentile")
  out <- list(method = "efron_events", theta_hat = db$beta_debiased, theta_star = theta_star, ci = ci, alpha = alpha, B = B)
  class(out) <- "selSurv_boot"
  out
}

#' @export
print.selSurv_boot <- function(x, ...) {
  cat("selSurv bootstrap object\n")
  cat("method =", x$method, " B =", x$B, " alpha =", x$alpha, "\n")
  print(utils::head(x$ci))
  invisible(x)
}
