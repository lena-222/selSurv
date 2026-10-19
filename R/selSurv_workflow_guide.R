#' Choosing an inference workflow in selSurv
#'
#' @description
#' A guide to inference after Cox-Lasso variable selection.
#'
#' Start with [coxlasso_fit()]. Subsequent methods differ in the
#' coefficient or procedure they describe and in whether variable
#' selection is repeated during resampling.
#'
#' @details
#' These workflows address different inferential questions.
#' An interval's interpretation depends on the target and the
#' procedure used to construct it.
#'
#' @section Which function should I use?:
#' * **Approximate full-model-coordinate inference.** Use
#'   [coxlasso_debias()] for a ridge-stabilized one-step correction.
#'   Optionally use [coxlasso_boot()] with
#'   `inference = "debiased"` and either `type = "wild"`
#'   or `type = "efron"` for a score-level bootstrap.
#'
#'   This is an approximate, research-stage method. Neither the
#'   one-step correction nor the bootstrap alone establishes
#'   general high-dimensional confidence-interval validity.
#'
#' * **Fixed-support model-based bootstrap.** Use [coxlasso_boot()] with
#'   `inference = "conditional"`. The default bootstrap type is
#'   `"model_based"`. The original selected support is held fixed
#'   throughout bootstrap refitting.
#'
#'   In this interface, "conditional" means that the *observed*
#'   support is kept fixed. No conditioning on the complete
#'   selection event is performed.
#'
#' * **Procedure-level re-selection.**
#'   Use [coxlasso_boot_reselect()] to compare fixed-support refitting
#'   with repeated Cox-Lasso selection. Available workflows include:
#'   - `"fixed_support"`: retain the original support;
#'   - `"reselect_fixed_lambda"`: reselect with a fixed penalty;
#'   - `"reselect_cv_min"`: reselect using CV-min;
#'   - `"reselect_cv_1se"`: reselect using the one-SE rule;
#'   - `"reselect_aic"`: reselect by AIC.
#'
#'   The reported bootstrap vectors assign zero to unselected
#'   coordinates; failures remain missing. Selection frequencies
#'   and model overlap describe stability. Resulting intervals
#'   need not have standard fixed-parameter coverage interpretations.
#'
#' @section Relation to the separate parametric interface:
#' [coxlasso_boot_parametric()] provides a separate model-based
#' interface, including an option to retune the penalty. Its
#' arguments, interval construction, and results differ from
#' those of the recommended interfaces above. Use it when its
#' specific workflow is needed; the functions are not aliases.
#'
#' @section Meaning of alpha and confidence-level arguments:
#' * In [coxlasso_fit()], `alpha` is the elastic-net mixing parameter.
#'   `alpha = 1` corresponds to Lasso.
#' * In [coxlasso_boot_reselect()], `alpha` is also the mixing parameter.
#' * In [coxlasso_split()], `penalty_alpha` is the mixing parameter;
#'   `conf_level` specifies the confidence level.
#' * In [coxlasso_boot_parametric()], `alpha` denotes the
#'   significance level, rather than the mixing parameter.
#' * The lower-level score bootstrap functions also use `alpha`
#'   for the significance level.
#'
#' Use named arguments when moving between these interfaces.
#'
#' @section Scope and tied event times:
#' The methods have been evaluated primarily for Lasso selection
#' (`alpha = 1`). [coxlasso_fit()] permits other mixing parameters.
#' The debiasing and score-bootstrap procedures require a Lasso
#' pilot. Other inference routes may run with non-Lasso penalties,
#' but those settings have not been evaluated separately.
#'
#' Selection and unpenalized refitting can use different methods
#' for handling tied event times. In particular, the current
#' one-step score calculation has not been established to agree
#' with the intended Cox score when failures are tied. Interpret
#' debiased results for tied times with caution.
#'
#' The Efron *event-score bootstrap* resamples event contributions;
#' it should not be confused with the Efron partial-likelihood
#' approximation for tied event times.
#'
#' @seealso
#' * [coxlasso_fit()]
#' * [coxlasso_debias()]
#' * [coxlasso_boot()]
#' * [coxlasso_boot_reselect()]
#' * [coxlasso_boot_parametric()]
#' * [coxlasso_split()]
#' @name selSurv_workflows
#' @keywords package
NULL
