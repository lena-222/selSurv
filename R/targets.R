#' Explain the implemented inferential targets
#'
#' @return A data frame describing the package workflows and targets.
#' @export
selSurv_targets <- function() {
  data.frame(
    target = c("full_model", "fixed_support", "procedure_level"),
    estimand = c(
      "Full-model Cox coefficient or debiased approximation.",
      "Cox projection coefficient on the support selected in the observed data.",
      "Zero-padded coefficient produced by the complete selection-refit workflow."
    ),
    non_selection = c("Not part of target", "Conditioned away", "Encoded as zero"),
    main_workflows = c("coxlasso_debias(), score bootstraps", "coxlasso_refit(), fixed_support bootstrap, split sample", "reselect_fixed_lambda or reselect_retuned model bootstrap"),
    stringsAsFactors = FALSE
  )
}

#' List implemented package methods
#'
#' @return A data frame with methods, targets, and status.
#' @export
selSurv_methods <- function() {
  data.frame(
    function_name = c("coxlasso_fit", "coxlasso_refit", "coxlasso_debias", "coxlasso_wildboot", "coxlasso_efron_events", "coxlasso_boot_parametric", "coxlasso_split"),
    target = c("selection", "fixed_support", "full_model", "debiased_estimator", "debiased_estimator", "fixed_support_or_procedure_level", "split_selected_model"),
    status = c("implemented", "implemented", "research implementation", "research implementation", "research implementation", "research implementation", "implemented"),
    stringsAsFactors = FALSE
  )
}
