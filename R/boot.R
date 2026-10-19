#' Unified bootstrap inference for Cox-Lasso models
#'
#' Computes bootstrap intervals for either an approximately debiased
#' coefficient vector or a fixed-support Cox refit. Select the
#' procedure with `inference` or `type`; the two routes have
#' different inferential interpretations.
#'
#' @param object A fitted Cox-Lasso object or an object returned by
#'   [coxlasso_debias()].
#' @param inference Optional inferential family. One of `"debiased"` or
#'   `"conditional"`.
#' @param type Bootstrap type: `"wild"` or `"efron"` for debiasing;
#'   `"model_based"` for fixed-support resampling. Compatible types may
#'   be combined. With an explicit `inference`, `"all"` selects every
#'   available type in that family.
#' @param B Number of bootstrap repetitions. Must be a positive integer.
#' @param conf_level Confidence level strictly between zero and one.
#' @param seed Optional non-negative integer seed. When multiple bootstrap
#'   types are requested, consecutive seeds are used in the supplied order.
#' @param ... Additional arguments forwarded to the selected backend engine.
#'
#' @details
#' Supplying `type = "wild"` or `"efron"` automatically selects
#' `inference = "debiased"`. Supplying `type = "model_based"` automatically
#' selects `inference = "conditional"`.
#'
#' Conversely, `inference = "debiased"` defaults to `type = "wild"`, while
#' `inference = "conditional"` defaults to `type = "model_based"`.
#'
#' For debiased bootstrap inference, a fitted object is automatically passed
#' through [coxlasso_debias()] before the bootstrap engine is run. Debiasing
#' remains independently available without bootstrap by calling
#' `coxlasso_debias(fit)`.
#'
#' Under `inference = "conditional"`, the support selected in the
#' observed data is held fixed during bootstrap refitting. This
#' does not condition on the complete selection event.
#'
#' For `inference = "debiased"`, the score bootstrap quantifies
#' variation around the implemented one-step correction; it does
#' not independently establish high-dimensional coverage validity.
#'
#' @return An object of class `selsurv_boot`. Its `results` element contains
#'   one `selsurv_boot_method` object per requested bootstrap type.
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
#' fit <- coxlasso_fit(X, time, status, nfolds = 3, lambda_choice = "min")
#'
#' # Illustrative B only; use many more draws for inference.
#' db <- coxlasso_debias(fit)
#' boot_score <- coxlasso_boot(fit, type = "wild", B = 8, seed = 2026)
#' boot_fixed <- coxlasso_boot(
#'   fit, inference = "conditional", B = 8, seed = 2026
#' )
#' summary(boot_fixed)
#' }
#'
#' @seealso \code{\link{selSurv_workflows}}
#' @export
coxlasso_boot <- function(
    object,
    inference = NULL,
    type = NULL,
    B = 500L,
    conf_level = 0.95,
    seed = NULL,
    ...
) {
  call <- match.call()
  dots <- list(...)
  .selsurv_validate_dots(dots)

  specification <- .selsurv_resolve_boot_spec(
    inference = inference,
    type = type
  )

  B <- .selsurv_validate_B(B)
  conf_level <- .selsurv_validate_conf_level(conf_level)
  seed <- .selsurv_validate_seed(seed)

  inference <- specification$inference
  types <- specification$type

  analysis_object <- .selsurv_prepare_boot_object(
    object = object,
    inference = inference
  )

  results <- vector("list", length(types))
  names(results) <- types

  for (i in seq_along(types)) {
    current_seed <- if (is.null(seed)) NULL else seed + i - 1L

    results[[i]] <- .selsurv_run_boot_engine(
      object = analysis_object,
      inference = inference,
      type = types[[i]],
      B = B,
      conf_level = conf_level,
      seed = current_seed,
      dots = dots
    )
  }

  structure(
    list(
      call = call,
      inference = inference,
      type = types,
      results = results,
      object = analysis_object,
      settings = list(
        B = B,
        conf_level = conf_level,
        seed = seed
      )
    ),
    class = c("selsurv_boot", "selsurv_result")
  )
}


#' @export
print.selsurv_boot <- function(x, ...) {
  cat("<selsurv_boot>\n")
  cat("Inference: ", x$inference, "\n", sep = "")
  cat("Types:     ", paste(x$type, collapse = ", "), "\n", sep = "")
  cat("B:         ", x$settings$B, "\n", sep = "")
  cat("Level:     ", x$settings$conf_level, "\n", sep = "")
  invisible(x)
}


#' @export
summary.selsurv_boot <- function(object, ...) {
  method_summaries <- lapply(
    object$results,
    function(x) {
      tryCatch(
        summary(x$result, ...),
        error = function(e) x$result
      )
    }
  )

  structure(
    list(
      inference = object$inference,
      type = object$type,
      settings = object$settings,
      methods = method_summaries,
      call = object$call
    ),
    class = "summary.selsurv_boot"
  )
}


#' @export
print.summary.selsurv_boot <- function(x, ...) {
  cat("selSurv bootstrap summary\n")
  cat("Inference: ", x$inference, "\n", sep = "")
  cat("Types:     ", paste(x$type, collapse = ", "), "\n", sep = "")
  cat("B:         ", x$settings$B, "\n\n", sep = "")

  for (type in x$type) {
    cat("--- ", type, " ---\n", sep = "")
    print(x$methods[[type]], ...)
    cat("\n")
  }

  invisible(x)
}
