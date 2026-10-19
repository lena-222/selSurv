# Experimental user-facing API for procedure-level Cox-Lasso re-selection.
#
# This API is intentionally separate from coxlasso_boot(inference = "conditional"),
# which remains the fixed-support model-based bootstrap in selSurv 0.1.0.

#' Model-based Cox bootstrap with re-selection
#'
#' Uses a model-based Cox bootstrap to compare fixed-support
#' refitting with procedures that repeat variable selection.
#' An unselected coefficient is recorded as zero; failed Cox
#' refits are recorded as missing, not as zero.
#'
#' @usage
#' coxlasso_boot_reselect(
#'   x, time, status,
#'   B = 500L,
#'   methods = c(
#'     "fixed_support", "reselect_fixed_lambda", "reselect_cv_min"
#'   ),
#'   conf_level = 0.95,
#'   seed = NULL,
#'   nfolds = 10L,
#'   original = NULL,
#'   original_support = NULL,
#'   original_lambda = NULL,
#'   standardize = TRUE,
#'   alpha = 1,
#'   type_measure = "deviance",
#'   ties = "efron",
#'   censoring = "reverse_km",
#'   rng_mode = c("independent_streams", "legacy_nki70"),
#'   glmnet_args = list(),
#'   progress = interactive()
#' )
#' @param x Numeric matrix or data frame of covariates. In data frames,
#'   numeric and integer columns remain numeric, while logical, factor,
#'   ordered-factor, and character columns are expanded using treatment
#'   coding with the first level as the reference. Ordered factors are treated
#'   as ordinary categorical predictors. Categorical predictors contribute
#'   \eqn{k-1} separately penalized design-matrix coordinates; they are not
#'   selected as groups. To control the reference category, supply an explicit
#'   factor level order or use `relevel()` before fitting.
#' @param time Follow-up times.
#' @param status Event indicator coded 0/1.
#' @param B Number of bootstrap replicates.
#' @param methods Character vector of bootstrap procedures.
#'   Supported names are listed in \code{\link{selSurv_workflows}}.
#' @param conf_level Confidence level.
#' @param seed Optional non-negative integer seed. If `NULL`, one random
#'   root seed is generated and stored in the returned object; all internal
#'   bootstrap and CV streams are then derived deterministically from it.
#' @param nfolds Number of cross-validation folds.
#' @param original Optional precomputed observed-data selection and refit
#'   object. This is an advanced reproducibility interface; most users should
#'   leave it as `NULL` and let the function fit the observed-data workflow.
#' @param original_support Observed-data support; provide it together
#'   with `original_lambda`.
#' @param original_lambda Optional observed-data Cox-Lasso penalty. Must be
#'   supplied together with `original_support`.
#' @param standardize Passed to glmnet. The default is `TRUE`; set to
#'   `FALSE` when `x` has already been standardized or when reproducing an
#'   analysis that used externally standardized covariates.
#' @param alpha Elastic-net mixing parameter for selection.
#'   Lasso (`alpha = 1`) is the primary evaluation setting.
#'   Other values are accepted with an inferential-scope warning.
#' @param type_measure Cross-validation loss passed to `cv.glmnet()`.
#' @param ties Ties method for the unpenalized Cox refits. The normal
#'   default is `"efron"`; exact `legacy_nki70` reproduction forces Breslow.
#' @param censoring Bootstrap censoring model; currently `"reverse_km"`.
#' @param rng_mode `"independent_streams"` for method-order-invariant random
#'   streams, or `"legacy_nki70"` only when reproducing the archived NKI70 run.
#' @param glmnet_args Additional arguments forwarded to glmnet fits.
#' @param progress Show bootstrap progress.
#'
#' @return An object of class `selsurv_reselection` and
#'   `cox_reselection_bootstrap`.
#'
#' @details
#' Re-selection characterizes variability in the selection-and-refit
#' procedure. It does not condition on the complete selection
#' event, and its intervals should not automatically be interpreted
#' as confidence intervals for a single fixed coefficient.
#'
#' The default RNG mode uses separate deterministic streams for bootstrap-data
#' generation and CV folds so that results for a named method do not depend on
#' which other methods are requested. If `seed = NULL`, a root seed is drawn
#' once at the start of the call and stored in the returned object, after which
#' the same deterministic stream construction is used.
#' `rng_mode = "legacy_nki70"` is retained
#' solely for reproducibility of the archived NKI70 B=500 analysis.
#'
#' @examples
#' \donttest{
#' set.seed(2026)
#' n <- 120L
#' X <- matrix(stats::rnorm(n * 4L), nrow = n)
#' colnames(X) <- paste0("x", seq_len(ncol(X)))
#' event_time <- stats::rexp(n, rate = 0.12 * exp(X[, 1] - 0.7 * X[, 2]))
#' censor_time <- stats::rexp(n, rate = 0.06)
#' time <- pmin(event_time, censor_time)
#' status <- as.integer(event_time <= censor_time)
#' # Illustrative B only; increase it before interpreting intervals.
#' rs <- coxlasso_boot_reselect(
#'   x = X, time = time, status = status,
#'   B = 8, nfolds = 3,
#'   methods = c(
#'     "fixed_support", "reselect_fixed_lambda",
#'     "reselect_cv_min"
#'   ),
#'   conf_level = 0.95, seed = 2026, progress = FALSE
#' )
#' summary(rs)$diagnostics
#' }
#' @seealso \code{\link{selSurv_workflows}}
#' @export
coxlasso_boot_reselect <- function(
    x, time, status,
    B = 500L,
    methods = c(
      "fixed_support",
      "reselect_fixed_lambda",
      "reselect_cv_min"
    ),
    conf_level = 0.95,
    seed = NULL,
    nfolds = 10L,
    original = NULL,
    original_support = NULL,
    original_lambda = NULL,
    standardize = TRUE,
    alpha = 1,
    type_measure = "deviance",
    ties = "efron",
    censoring = "reverse_km",
    rng_mode = c("independent_streams", "legacy_nki70"),
    glmnet_args = list(),
    progress = interactive()
) {
  ans <- run_reselection_bootstrap(
    x = x, time = time, status = status,
    B = B, methods = methods,
    conf_level = conf_level, seed = seed,
    nfolds = nfolds,
    original = original,
    original_support = original_support,
    original_lambda = original_lambda,
    standardize = standardize,
    alpha = alpha,
    type_measure = type_measure,
    ties = ties,
    censoring = censoring,
    rng_mode = match.arg(rng_mode),
    glmnet_args = glmnet_args,
    progress = progress
  )
  class(ans) <- c("selsurv_reselection", class(ans))
  ans
}

#' @export
summary.selsurv_reselection <- function(object, ...) {
  structure(
    list(
      call = object$call,
      intervals = summarise_reselection(object),
      diagnostics = method_level_diagnostics(object),
      settings = object$settings
    ),
    class = "summary.selsurv_reselection"
  )
}

#' @export
print.selsurv_reselection <- function(x, ...) {
  NextMethod("print", x, ...)
}

#' @export
print.summary.selsurv_reselection <- function(x, ...) {
  cat("selSurv procedure-level Cox-Lasso re-selection bootstrap\n")
  cat("RNG mode:", x$settings$rng_mode, "\n\n")
  cat("Method-level diagnostics:\n")
  print(x$diagnostics, row.names = FALSE)
  cat("\nCoefficient intervals for originally selected coordinates:\n")
  keep <- c(
    "method", "variable", "estimate",
    "bootstrap_se_lower", "bootstrap_se_upper",
    "reselection_frequency", "boot_zero_frequency"
  )
  print(x$intervals[, intersect(keep, names(x$intervals)), drop = FALSE],
        row.names = FALSE)
  invisible(x)
}
