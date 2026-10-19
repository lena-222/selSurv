.selsurv_validate_design_names <- function(nm, what = "design-matrix column names") {
  if (is.null(nm)) {
    return(invisible(NULL))
  }

  if (anyNA(nm) || any(!nzchar(nm))) {
    stop(what, " must be non-missing and non-empty.", call. = FALSE)
  }

  if (anyDuplicated(nm)) {
    dup <- unique(nm[duplicated(nm)])
    stop(
      what, " must be unique. Duplicated name(s): ",
      paste(dup, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  bad_syntax <- grepl("`|[[:cntrl:]]", nm)

  if (any(bad_syntax)) {
    stop(
      what,
      " must not contain backticks or control characters. Problematic name(s): ",
      paste(nm[bad_syntax], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  reserved <- nm %in% c(".time", ".status")

  if (any(reserved)) {
    stop(
      what, " must not use the reserved name(s) ",
      paste(nm[reserved], collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(nm)
}


.selsurv_as_design_matrix <- function(X, arg = "X") {
  if (is.data.frame(X)) {
    if (nrow(X) < 1L) {
      stop(arg, " must contain at least one row.", call. = FALSE)
    }

    if (ncol(X) < 1L) {
      stop(arg, " must contain at least one covariate.", call. = FALSE)
    }

    if (anyNA(X)) {
      stop(
        arg, " must contain only finite, non-missing covariate values.",
        call. = FALSE
      )
    }

    original_names <- names(X)

    .selsurv_validate_design_names(
      original_names,
      paste0(arg, " column names")
    )

    dat <- lapply(seq_along(X), function(j) {
      x <- X[[j]]
      nm <- original_names[[j]]

      if (is.logical(x) &&
          is.null(dim(x)) &&
          is.null(attr(x, "class"))) {
        # Logical predictors are binary categorical predictors.
        # Constant logical columns are rejected below.
        x <- factor(x, levels = c(FALSE, TRUE))
      } else if (is.ordered(x)) {
        # Ordered factors are treated as categorical predictors rather than
        # receiving the default polynomial contrasts.
        x <- factor(
          as.character(x),
          levels = levels(x),
          ordered = FALSE
        )
      } else if (is.factor(x)) {
        # Keep the user-supplied factor level order.
        x <- x
      } else if (is.character(x) &&
                 is.null(dim(x)) &&
                 is.null(attr(x, "class"))) {
        # Character predictors use ordinary unordered-factor coding.
        x <- factor(x)
      } else {
        plain_numeric <- is.numeric(x) &&
          is.null(dim(x)) &&
          is.null(attr(x, "class"))

        if (!plain_numeric) {
          stop(
            arg, " column `", nm, "` has unsupported class `",
            paste(class(x), collapse = "/"),
            "`. Supported data-frame columns are plain numeric/integer, ",
            "logical, factor, ordered factor, and character.",
            call. = FALSE
          )
        }
      }

      if (is.factor(x)) {
        # Unused factor levels do not define estimable coordinates.
        x <- droplevels(x)

        if (nlevels(x) < 2L) {
          stop(
            arg, " column `", nm,
            "` must contain at least two observed categorical levels.",
            call. = FALSE
          )
        }
      }

      x
    })

    dat <- as.data.frame(
      dat,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    names(dat) <- original_names

    factor_names <- names(dat)[
      vapply(dat, is.factor, logical(1))
    ]

    contrasts_arg <- if (length(factor_names)) {
      stats::setNames(
        lapply(
          factor_names,
          function(nm) {
            stats::contr.treatment(
              levels(dat[[nm]]),
              base = 1L
            )
          }
        ),
        factor_names
      )
    } else {
      NULL
    }

    # Build ordinary treatment coding with an intercept and remove only the
    # intercept afterwards. Every categorical predictor therefore contributes
    # k - 1 coordinates, independent of its position in the data frame.
    X <- stats::model.matrix(
      ~ .,
      data = dat,
      contrasts.arg = contrasts_arg,
      na.action = stats::na.fail
    )

    X <- X[
      ,
      colnames(X) != "(Intercept)",
      drop = FALSE
    ]
  } else {
    X <- as.matrix(X)
  }

  storage.mode(X) <- "double"

  if (nrow(X) < 1L) {
    stop(arg, " must contain at least one row.", call. = FALSE)
  }

  if (ncol(X) < 1L) {
    stop(arg, " must contain at least one usable covariate.", call. = FALSE)
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }

  .selsurv_validate_design_names(
    colnames(X),
    paste0(arg, " design-matrix column names")
  )

  if (anyNA(X) || any(!is.finite(X))) {
    stop(
      arg, " must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }

  X
}

check_survival_inputs <- function(X, time, status) {
  X <- .selsurv_as_design_matrix(X, arg = "X")
  time <- as.numeric(time)
  status <- as.integer(status)

  if (nrow(X) != length(time) || length(time) != length(status)) {
    stop("X, time, and status must have matching row/length dimensions.", call. = FALSE)
  }
  if (any(!is.finite(time)) || any(time <= 0)) {
    stop("time must contain positive finite values.", call. = FALSE)
  }
  if (!all(status %in% c(0L, 1L))) {
    stop("status must be coded as 0/1.", call. = FALSE)
  }

  list(X = X, time = time, status = status)
}

selected_from_beta <- function(beta, tol = 1e-8) {
  which(abs(as.numeric(beta)) > tol)
}

zero_pad <- function(values, support, p, names = NULL) {
  out <- rep(0, p)
  if (length(support) > 0) {
    out[support] <- values
  }
  if (!is.null(names)) {
    names(out) <- names
  }
  out
}

safe_match_choice <- function(x, choices) {
  x <- match.arg(x, choices)
  x
}

quantile_names <- function(alpha) {
  c(lower = alpha / 2, upper = 1 - alpha / 2)
}

make_interval <- function(theta_hat, theta_star, alpha = 0.1, type = c("percentile", "basic")) {
  type <- match.arg(type)
  q <- apply(theta_star, 2, stats::quantile, probs = quantile_names(alpha), na.rm = TRUE, names = FALSE)
  if (is.null(dim(q))) {
    q <- matrix(q, nrow = 2)
  }
  if (type == "percentile") {
    lower <- q[1, ]
    upper <- q[2, ]
  } else {
    lower <- 2 * theta_hat - q[2, ]
    upper <- 2 * theta_hat - q[1, ]
  }
  data.frame(term = names(theta_hat), estimate = as.numeric(theta_hat), lower = as.numeric(lower), upper = as.numeric(upper), row.names = NULL)
}

as_term_dataframe <- function(theta, lower = NULL, upper = NULL) {
  out <- data.frame(term = names(theta), estimate = as.numeric(theta), row.names = NULL)
  if (!is.null(lower)) out$lower <- as.numeric(lower)
  if (!is.null(upper)) out$upper <- as.numeric(upper)
  out
}

cox_formula_from_matrix <- function(time, status, X, support) {
  dat <- data.frame(.time = time, .status = status)
  if (length(support) > 0) {
    for (j in support) dat[[colnames(X)[j]]] <- X[, j]
    rhs <- paste(sprintf("`%s`", colnames(X)[support]), collapse = " + ")
  } else {
    rhs <- "1"
  }
  list(formula = stats::as.formula(paste0("survival::Surv(.time, .status) ~ ", rhs)), data = dat)
}

fit_cox_support <- function(X, time, status, support, ties = "efron") {
  fdat <- cox_formula_from_matrix(time, status, X, support)
  if (length(support) == 0) {
    return(list(coef = numeric(0), fit = NULL, support = integer(0), converged = TRUE))
  }
  fit <- tryCatch(
    survival::coxph(fdat$formula, data = fdat$data, ties = ties, x = TRUE, model = TRUE),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(coef = rep(NA_real_, length(support)), fit = NULL, support = support, converged = FALSE, error = conditionMessage(fit)))
  }
  cf <- stats::coef(fit)
  names(cf) <- colnames(X)[support]
  se <- tryCatch(sqrt(diag(stats::vcov(fit))), error = function(e) rep(NA_real_, length(cf)))
  names(se) <- names(cf)
  list(coef = cf, se = se, fit = fit, support = support, converged = TRUE)
}

km_censor_sampler <- function(time, status) {
  censor_status <- as.integer(status == 0)
  if (sum(censor_status) < 2) {
    rate <- 1 / mean(time)
    return(function(n) stats::rexp(n, rate = rate))
  }
  sf <- survival::survfit(survival::Surv(time, censor_status) ~ 1)
  times <- sf$time
  surv <- sf$surv
  cdf <- 1 - surv
  cdf <- pmin(pmax(cdf, 0), 1)
  keep <- !duplicated(cdf) & is.finite(cdf)
  times <- times[keep]
  cdf <- cdf[keep]
  if (length(times) < 2 || max(cdf) <= 0) {
    return(function(n) sample(time, n, replace = TRUE))
  }
  function(n) {
    u <- stats::runif(n, min = 0, max = max(cdf))
    stats::approx(x = cdf, y = times, xout = u, rule = 2, ties = "ordered")$y
  }
}

inverse_baseline_time <- function(basehaz_df, target_hazard) {
  bh <- basehaz_df
  bh <- bh[is.finite(bh$hazard) & is.finite(bh$time), , drop = FALSE]
  if (nrow(bh) < 2 || max(bh$hazard) <= 0) {
    rate <- 1 / max(1, max(bh$time, na.rm = TRUE))
    return(stats::rexp(length(target_hazard), rate = rate))
  }
  hazard <- pmax(bh$hazard, 0)
  time <- pmax(bh$time, 0)
  keep <- !duplicated(hazard)
  hazard <- hazard[keep]
  time <- time[keep]
  max_h <- max(hazard)
  max_t <- max(time)
  out <- stats::approx(x = hazard, y = time, xout = pmin(target_hazard, max_h), rule = 2, ties = "ordered")$y
  beyond <- target_hazard > max_h
  if (any(beyond)) {
    tail_rate <- max(1e-8, max_h / max(max_t, 1e-8))
    out[beyond] <- max_t + stats::rexp(sum(beyond), rate = tail_rate)
  }
  out
}

# Compute the Breslow cumulative baseline hazard at the coefficients of a
# fitted Cox model.  This is deliberately independent of the coefficient-fit
# tie convention: an Efron Cox refit may still be used as the plug-in
# coefficient estimator while the bootstrap generator uses a Breslow baseline.
.selsurv_breslow_basehaz <- function(fit, centered = FALSE) {
  if (!inherits(fit, "coxph")) {
    stop("fit must inherit from `coxph`.", call. = FALSE)
  }

  sf <- survival::survfit(
    fit,
    se.fit = FALSE,
    ctype = 1
  )

  hazard <- as.numeric(sf$cumhaz)

  if (!isTRUE(centered)) {
    beta <- stats::coef(fit)
    beta <- ifelse(is.na(beta), 0, beta)

    means <- fit$means
    if (is.null(means)) {
      means <- numeric(0)
    }

    offset <- sum(means * beta)
    hazard <- hazard * exp(-offset)
  }

  out <- data.frame(
    hazard = hazard,
    time = as.numeric(sf$time)
  )

  if (!is.null(sf$strata)) {
    out$strata <- factor(
      rep(names(sf$strata), sf$strata),
      levels = names(sf$strata)
    )
  }

  out
}


# Prepare glmnet arguments for Cox models using selSurv's explicit
# selection convention.  glmnet < 5 does not expose cox.ties and uses
# Breslow internally.  Newer versions can receive cox.ties explicitly.
.selsurv_prepare_cox_glmnet_args <- function(
    args,
    selection_ties = c("breslow", "efron"),
    cv = FALSE,
    glmnet_formals = names(formals(glmnet::glmnet))
) {
  if (!is.list(args)) {
    stop("args must be a list.", call. = FALSE)
  }

  selection_ties <- match.arg(selection_ties)
  supports_cox_ties <- "cox.ties" %in% glmnet_formals

  if (!supports_cox_ties &&
      identical(selection_ties, "efron")) {
    stop(
      "This glmnet version does not support Efron tie handling for Cox ",
      "fits. Use Breslow tie handling or a glmnet version exposing ",
      "`cox.ties`.",
      call. = FALSE
    )
  }

  if ("cox.ties" %in% names(args)) {
    requested_ties <- match.arg(
      as.character(args[["cox.ties"]]),
      c("breslow", "efron")
    )

    if (!supports_cox_ties &&
        !identical(requested_ties, "breslow")) {
      stop(
        "This glmnet version does not support Efron tie handling for Cox ",
        "fits. Use Breslow tie handling or a glmnet version exposing ",
        "`cox.ties`.",
        call. = FALSE
      )
    }

    # glmnet < 5 silently accepts cox.ties through ... without using it.
    # Remove an explicit Breslow request there rather than forwarding an
    # ineffective argument.
    if (!supports_cox_ties) {
      args[["cox.ties"]] <- NULL
    }
  } else if (supports_cox_ties) {
    args[["cox.ties"]] <- selection_ties
  }

  if (isTRUE(cv)) {
    if (!"type.measure" %in% names(args)) {
      args[["type.measure"]] <- "deviance"
    }
    if (!"grouped" %in% names(args)) {
      args[["grouped"]] <- TRUE
    }
  }

  args
}

# Penalization-scope safeguards --------------------------------------------

.selsurv_is_lasso_alpha <- function(alpha, tol = 1e-12) {
  is.numeric(alpha) &&
    length(alpha) == 1L &&
    !is.na(alpha) &&
    is.finite(alpha) &&
    abs(alpha - 1) <= tol
}

.selsurv_require_lasso_fit <- function(fit, caller) {
  alpha <- fit$alpha

  if (is.null(alpha) ||
      length(alpha) != 1L ||
      !is.numeric(alpha) ||
      is.na(alpha) ||
      !is.finite(alpha)) {
    stop(
      caller,
      " requires a fitted object with a recorded Cox penalty mixing ",
      "parameter. A Lasso fit with `alpha = 1` is required.",
      call. = FALSE
    )
  }

  if (!.selsurv_is_lasso_alpha(alpha)) {
    stop(
      caller,
      " currently supports Cox-Lasso fits with `alpha = 1` only. ",
      "Debiased inference for Elastic-Net or Ridge pilots has not been ",
      "validated in selSurv.",
      call. = FALSE
    )
  }

  invisible(alpha)
}

.selsurv_warn_nonlasso <- function(
    alpha,
    caller,
    argument = "alpha"
) {
  if (is.null(alpha) ||
      length(alpha) != 1L ||
      !is.numeric(alpha) ||
      is.na(alpha) ||
      !is.finite(alpha)) {
    return(invisible(FALSE))
  }

  if (!.selsurv_is_lasso_alpha(alpha)) {
    warning(
      caller,
      " is proceeding with `",
      argument,
      " = ",
      format(alpha, digits = 8),
      "`. selSurv's inferential validation currently concerns Lasso ",
      "screening (`alpha = 1`). Other mixing parameters are available ",
      "for fitting flexibility but are not separately validated ",
      "inferentially.",
      call. = FALSE
    )
    return(invisible(TRUE))
  }

  invisible(FALSE)
}
