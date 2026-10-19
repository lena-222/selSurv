#' Compare selection diagnostics across bootstrap workflows
#'
#' Combines the `diagnostics` tables returned by `coxlasso_boot_parametric()`
#' for fixed-support, fixed-lambda re-selection, and retuned re-selection
#' workflows.
#'
#' @param fixed Optional result from `coxlasso_boot_parametric(..., workflow = "fixed_support")`.
#' @param reselect_fixed_lambda Optional result from `coxlasso_boot_parametric(..., workflow = "reselect_fixed_lambda")`.
#' @param reselect_retuned Optional result from `coxlasso_boot_parametric(..., workflow = "reselect_retuned")`.
#' @param digits Number of digits used for rounding. Use `NULL` for no rounding.
#'
#' @return A data frame with one row per term and workflow-specific re-selection and zero-frequency columns.
#' @export
compare_selection_diagnostics <- function(fixed = NULL,
                                          reselect_fixed_lambda = NULL,
                                          reselect_retuned = NULL,
                                          digits = 2) {
  pieces <- list()

  add_piece <- function(object, prefix) {
    if (is.null(object)) {
      return(NULL)
    }
    if (is.null(object$diagnostics)) {
      stop("Object for prefix '", prefix, "' has no $diagnostics component.", call. = FALSE)
    }

    d <- object$diagnostics
    required <- c("term", "reselect_frequency", "zero_frequency")
    missing <- setdiff(required, names(d))
    if (length(missing) > 0) {
      stop(
        "Diagnostics for prefix '", prefix, "' are missing columns: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }

    out <- d[, required, drop = FALSE]
    names(out) <- c("term", paste0(prefix, "_reselect"), paste0(prefix, "_zero"))
    out
  }

  pieces$fixed <- add_piece(fixed, "fixed")
  pieces$t2 <- add_piece(reselect_fixed_lambda, "t2")
  pieces$retuned <- add_piece(reselect_retuned, "retuned")
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]

  if (length(pieces) == 0) {
    stop("At least one bootstrap result must be supplied.", call. = FALSE)
  }

  ans <- Reduce(function(a, b) merge(a, b, by = "term", all = TRUE), pieces)
  ans <- ans[order(ans$term), , drop = FALSE]
  rownames(ans) <- NULL

  if (!is.null(digits)) {
    num_cols <- vapply(ans, is.numeric, logical(1))
    ans[num_cols] <- lapply(ans[num_cols], round, digits = digits)
  }

  ans
}
