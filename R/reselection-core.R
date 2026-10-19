
# Cox-Lasso model-based bootstrap with re-selection
# Reconstructed to match the current reselection manuscript.
#
# Core principles:
#   * original analysis: CV-min Cox-Lasso -> support Mhat -> unpenalized Cox refit
#   * bootstrap generator fitted ONCE from the original unpenalized refit
#   * observed X kept fixed
#   * same generated bootstrap dataset is used across compared bootstrap methods
#   * methods:
#       - fixed_support
#       - reselect_fixed_lambda
#       - reselect_cv_min
#       - reselect_cv_1se (tuning-rule sensitivity)
#       - reselect_aic (tuning-rule sensitivity)
#   * re-selection coefficients are recorded on the original p-dimensional scale:
#       non-selected coordinate = 0
#   * failed non-empty Cox refits are NA, NOT zero
#   * main comparison = bootstrap-SE Wald on the zero-padded distribution
#   * the observed workflow can be supplied explicitly, avoiding accidental
#       re-selection/re-tuning before the bootstrap starts
#   * two RNG modes are available: a method-order-invariant default and an
#       archived NKI70 reproduction mode
#   * reverse-percentile basic, percentile, and studentized intervals = sensitivities
#
# Packages: survival, glmnet

.required_pkgs <- c("survival", "glmnet")

check_reselection_packages <- function() {
  miss <- .required_pkgs[!vapply(.required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("Missing required package(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.reselection_prepare_cox_glmnet_args <- function(
    args,
    cv = FALSE
) {
  # Normal package execution delegates to the shared helper.
  if (exists(
        ".selsurv_prepare_cox_glmnet_args",
        mode = "function",
        inherits = TRUE
      )) {
    return(
      .selsurv_prepare_cox_glmnet_args(
        args,
        selection_ties = "breslow",
        cv = cv
      )
    )
  }

  # Standalone compatibility for archived replication scripts that source
  # reselection-core.R directly.
  supports_cox_ties <- "cox.ties" %in%
    names(formals(glmnet::glmnet))

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

    if (!supports_cox_ties) {
      args[["cox.ties"]] <- NULL
    }
  } else if (supports_cox_ties) {
    args[["cox.ties"]] <- "breslow"
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


.reselection_breslow_basehaz <- function(
    fit,
    centered = FALSE
) {
  if (exists(
        ".selsurv_breslow_basehaz",
        mode = "function",
        inherits = TRUE
      )) {
    return(
      .selsurv_breslow_basehaz(
        fit,
        centered = centered
      )
    )
  }

  # Standalone replication fallback.
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


.reselection_resolve_refit_ties <- function(
    ties,
    rng_mode
) {
  ties <- match.arg(
    ties,
    c("efron", "breslow", "exact")
  )

  # Exact archived NKI70 reproduction always uses Breslow Cox refits.
  if (identical(rng_mode, "legacy_nki70")) {
    return("breslow")
  }

  ties
}


as_numeric_matrix <- function(x) {
  # Normal package execution uses the common input normalizer from utils.R.
  if (exists(".selsurv_as_design_matrix", mode = "function", inherits = TRUE)) {
    return(.selsurv_as_design_matrix(x, arg = "x"))
  }

  # Standalone compatibility for archived replication scripts that source
  # reselection-core.R directly.
  validate_names <- function(nm, what) {
    if (is.null(nm)) return(invisible(NULL))

    if (anyNA(nm) || any(!nzchar(nm))) {
      stop(what, " must be non-missing and non-empty.", call. = FALSE)
    }

    if (anyDuplicated(nm)) {
      stop(what, " must be unique.", call. = FALSE)
    }

    if (any(grepl("`|[[:cntrl:]]", nm))) {
      stop(
        what, " must not contain backticks or control characters.",
        call. = FALSE
      )
    }

    if (any(nm %in% c(".time", ".status"))) {
      stop(
        what, " must not use the reserved names .time or .status.",
        call. = FALSE
      )
    }

    invisible(nm)
  }

  if (is.data.frame(x)) {
    if (nrow(x) < 1L || ncol(x) < 1L) {
      stop("x must contain at least one row and one covariate.", call. = FALSE)
    }

    if (anyNA(x)) {
      stop(
        "x must contain only finite, non-missing covariate values.",
        call. = FALSE
      )
    }

    original_names <- names(x)
    validate_names(original_names, "x column names")

    dat <- lapply(seq_along(x), function(j) {
      z <- x[[j]]
      nm <- original_names[[j]]

      if (is.logical(z) &&
          is.null(dim(z)) &&
          is.null(attr(z, "class"))) {
        z <- factor(z, levels = c(FALSE, TRUE))
      } else if (is.ordered(z)) {
        z <- factor(
          as.character(z),
          levels = levels(z),
          ordered = FALSE
        )
      } else if (is.factor(z)) {
        z <- z
      } else if (is.character(z) &&
                 is.null(dim(z)) &&
                 is.null(attr(z, "class"))) {
        z <- factor(z)
      } else {
        plain_numeric <- is.numeric(z) &&
          is.null(dim(z)) &&
          is.null(attr(z, "class"))

        if (!plain_numeric) {
          stop(
            "x column `", nm, "` has unsupported class `",
            paste(class(z), collapse = "/"), "`.",
            call. = FALSE
          )
        }
      }

      if (is.factor(z)) {
        z <- droplevels(z)

        if (nlevels(z) < 2L) {
          stop(
            "x column `", nm,
            "` must contain at least two observed categorical levels.",
            call. = FALSE
          )
        }
      }

      z
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

    x <- stats::model.matrix(
      ~ .,
      data = dat,
      contrasts.arg = contrasts_arg,
      na.action = stats::na.fail
    )

    x <- x[
      ,
      colnames(x) != "(Intercept)",
      drop = FALSE
    ]
  } else {
    x <- as.matrix(x)
  }

  storage.mode(x) <- "double"

  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(
      "x must contain at least one row and one usable covariate.",
      call. = FALSE
    )
  }

  if (is.null(colnames(x))) {
    colnames(x) <- paste0("x", seq_len(ncol(x)))
  }

  validate_names(colnames(x), "x design-matrix column names")

  if (anyNA(x) || any(!is.finite(x))) {
    stop(
      "x must contain only finite, non-missing numeric values.",
      call. = FALSE
    )
  }

  x
}

check_survival_input <- function(x, time, status) {
  x <- as_numeric_matrix(x)
  time <- as.numeric(time)
  status <- as.integer(status)

  if (nrow(x) != length(time) || length(time) != length(status)) {
    stop("nrow(x), length(time), and length(status) must agree.", call. = FALSE)
  }
  if (anyNA(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("time must be finite, non-missing, and non-negative.", call. = FALSE)
  }
  if (anyNA(status) || any(!status %in% c(0L, 1L))) {
    stop("status must contain only 0/1.", call. = FALSE)
  }
  if (sum(status) < 2L) {
    stop("At least two observed events are required.", call. = FALSE)
  }

  list(x = x, time = time, status = status)
}

balanced_foldid <- function(n, nfolds = 10L) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
      n < 3 || n != floor(n)) {
    stop("n must be an integer >= 3 for cross-validation.", call. = FALSE)
  }
  if (length(nfolds) != 1L || !is.numeric(nfolds) || !is.finite(nfolds) ||
      nfolds < 3 || nfolds != floor(nfolds)) {
    stop("nfolds must be an integer >= 3.", call. = FALSE)
  }

  n <- as.integer(n)
  nfolds <- min(as.integer(nfolds), n)
  sample(rep(seq_len(nfolds), length.out = n))
}

.with_reselection_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(seed)
  force(code)
}

.reselection_resolve_root_seed <- function(
    seed,
    rng_mode
) {
  if (is.null(seed)) {
    if (identical(rng_mode, "legacy_nki70")) {
      stop(
        "`rng_mode = \"legacy_nki70\"` requires an explicit seed for ",
        "reproducibility.",
        call. = FALSE
      )
    }

    # Draw one root seed from the caller's RNG state. All subsequent
    # bootstrap-generation and CV streams are deterministic functions of this
    # root, preserving method-order invariance within the re-selection call.
    return(
      as.integer(
        sample.int(.Machine$integer.max, size = 1L) - 1L
      )
    )
  }

  .reselection_seed(seed)
}


.reselection_seed <- function(seed, replicate = 0L, stream = 0L) {
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
      !is.finite(seed) || seed < 0 || seed > .Machine$integer.max ||
      seed != floor(seed)) {
    stop(
      "seed must be one non-negative integer not exceeding .Machine$integer.max.",
      call. = FALSE
    )
  }

  modulus <- as.double(.Machine$integer.max) + 1
  as.integer((as.double(seed) + 100003 * as.double(replicate) +
                1009 * as.double(stream)) %% modulus)
}

.normalize_reselection_methods <- function(methods) {
  aliases <- c(
    reselect_fixed = "reselect_fixed_lambda",
    reselect_cvmin = "reselect_cv_min",
    reselect_cv1se = "reselect_cv_1se",
    reselect_retuned = "reselect_cv_min"
  )

  methods <- as.character(methods)
  if (!length(methods) || anyNA(methods) || any(!nzchar(methods))) {
    stop("methods must contain at least one non-missing method name.", call. = FALSE)
  }

  hit <- methods %in% names(aliases)
  methods[hit] <- unname(aliases[methods[hit]])

  allowed <- c(
    "fixed_support",
    "reselect_fixed_lambda",
    "reselect_cv_min",
    "reselect_cv_1se",
    "reselect_aic"
  )

  bad <- setdiff(methods, allowed)
  if (length(bad)) {
    stop(
      "Unknown re-selection method(s): ",
      paste(shQuote(bad), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  unique(methods)
}

safe_coef_glmnet <- function(fit, s) {
  b <- as.matrix(stats::coef(fit, s = s))
  drop(b[, 1L])
}

support_from_beta <- function(beta, tol = 1e-10) {
  which(is.finite(beta) & abs(beta) > tol)
}

support_key <- function(support) {
  support <- sort(unique(as.integer(support)))
  if (!length(support)) "<EMPTY>" else paste(support, collapse = ",")
}

jaccard_support <- function(a, b) {
  a <- unique(as.integer(a))
  b <- unique(as.integer(b))
  u <- union(a, b)
  if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}

same_support <- function(a, b) {
  identical(sort(unique(as.integer(a))), sort(unique(as.integer(b))))
}


# -------------------------------------------------------------------------
# 1. Original Cox-Lasso selection and unpenalized Cox refit
# -------------------------------------------------------------------------

fit_original_selection <- function(
    x, time, status,
    nfolds = 10L,
    foldid = NULL,
    standardize = FALSE,
    alpha = 1,
    type_measure = "deviance",
    glmnet_args = list()
) {
  check_reselection_packages()
  z <- check_survival_input(x, time, status)
  x <- z$x; time <- z$time; status <- z$status
  y <- survival::Surv(time, status)

  if (is.null(foldid)) foldid <- balanced_foldid(nrow(x), nfolds)

  args <- c(
    list(
      x = x,
      y = y,
      family = "cox",
      alpha = alpha,
      foldid = foldid,
      type.measure = type_measure,
      standardize = standardize
    ),
    glmnet_args
  )
  args <- .reselection_prepare_cox_glmnet_args(
    args,
    cv = TRUE
  )
  cvfit <- do.call(glmnet::cv.glmnet, args)
  lambda_min <- unname(cvfit$lambda.min)
  beta_lasso <- safe_coef_glmnet(cvfit, s = "lambda.min")
  names(beta_lasso) <- colnames(x)
  support <- support_from_beta(beta_lasso)

  list(
    cvfit = cvfit,
    lambda_min = lambda_min,
    beta_lasso = beta_lasso,
    support = support,
    foldid = foldid,
    standardize = standardize,
    alpha = alpha
  )
}


fit_unpenalized_cox <- function(
    x, time, status, support,
    ties = "efron",
    singular_ok = FALSE
) {
  check_reselection_packages()
  x <- as_numeric_matrix(x)
  time <- as.numeric(time)
  status <- as.integer(status)
  p <- ncol(x)
  support <- sort(unique(as.integer(support)))

  beta_pad <- setNames(rep(0, p), colnames(x))
  se_pad <- setNames(rep(NA_real_, p), colnames(x))

  # Empty selected model is a legitimate selection outcome.
  if (!length(support)) {
    dat <- data.frame(.time = time, .status = status)
    fit <- tryCatch(
      survival::coxph(
        survival::Surv(.time, .status) ~ 1,
        data = dat,
        ties = ties,
        x = TRUE,
        model = TRUE
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      return(list(
        success = FALSE, empty = TRUE, fit = NULL,
        beta = rep(NA_real_, p), se = se_pad,
        support = support, error = conditionMessage(fit)
      ))
    }
    return(list(
      success = TRUE, empty = TRUE, fit = fit,
      beta = beta_pad, se = se_pad,
      support = support, error = NA_character_
    ))
  }

  if (any(support < 1L | support > p)) {
    stop("support contains invalid column indices.", call. = FALSE)
  }

  xs <- as.data.frame(x[, support, drop = FALSE], check.names = FALSE)
  dat <- cbind(data.frame(.time = time, .status = status), xs)
  rhs <- paste(sprintf("`%s`", colnames(xs)), collapse = " + ")
  fml <- stats::as.formula(paste("survival::Surv(.time, .status) ~", rhs))

  warnings_seen <- character()
  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(
        fml, data = dat, ties = ties,
        x = TRUE, model = TRUE,
        singular.ok = singular_ok
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      success = FALSE, empty = FALSE, fit = NULL,
      beta = setNames(rep(NA_real_, p), colnames(x)),
      se = se_pad, support = support,
      error = conditionMessage(fit),
      warnings = warnings_seen
    ))
  }

  b <- stats::coef(fit)
  vv <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  s <- if (is.null(vv)) rep(NA_real_, length(b)) else sqrt(diag(vv))

  bad <- any(!is.finite(b)) || length(b) != length(support)
  if (bad) {
    return(list(
      success = FALSE, empty = FALSE, fit = fit,
      beta = setNames(rep(NA_real_, p), colnames(x)),
      se = se_pad, support = support,
      error = "Non-finite or incomplete unpenalized Cox refit.",
      warnings = warnings_seen
    ))
  }

  beta_pad[support] <- as.numeric(b)
  se_pad[support] <- as.numeric(s)

  list(
    success = TRUE, empty = FALSE, fit = fit,
    beta = beta_pad, se = se_pad,
    support = support, error = NA_character_,
    warnings = warnings_seen
  )
}


# Exact unpenalized Cox-refit convention used by the archived final NKI70
# production script. This is intentionally isolated from the clean/reference
# refit because the historical script used a matrix term (Xs), x=TRUE, y=TRUE,
# coxph.control(iter.max=50), and the default singular.ok behaviour.
fit_unpenalized_cox_legacy_nki70 <- function(x, time, status, support) {
  check_reselection_packages()
  x <- as_numeric_matrix(x)
  time <- as.numeric(time)
  status <- as.integer(status)
  p <- ncol(x)
  support <- sort(unique(as.integer(support)))

  beta_pad <- setNames(rep(0, p), colnames(x))
  se_pad <- setNames(rep(NA_real_, p), colnames(x))

  if (!length(support)) {
    return(list(
      success = TRUE, empty = TRUE, fit = NULL,
      beta = beta_pad, se = se_pad, support = support,
      error = NA_character_, warnings = character()
    ))
  }

  if (any(support < 1L | support > p)) {
    stop("support contains invalid column indices.", call. = FALSE)
  }

  Xs <- x[, support, drop = FALSE]
  warnings_seen <- character()

  fit <- tryCatch(
    withCallingHandlers(
      survival::coxph(
        survival::Surv(time, status) ~ Xs,
        ties = "breslow",
        x = TRUE,
        y = TRUE,
        control = survival::coxph.control(iter.max = 50)
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  if (inherits(fit, "error")) {
    return(list(
      success = FALSE, empty = FALSE, fit = NULL,
      beta = setNames(rep(NA_real_, p), colnames(x)),
      se = se_pad, support = support,
      error = conditionMessage(fit), warnings = warnings_seen
    ))
  }

  cc <- unname(stats::coef(fit))
  vv <- tryCatch(diag(stats::vcov(fit)),
                 error = function(e) rep(NA_real_, length(cc)))
  ss <- sqrt(vv)

  if (length(cc) != length(support) || any(!is.finite(cc))) {
    return(list(
      success = FALSE, empty = FALSE, fit = fit,
      beta = setNames(rep(NA_real_, p), colnames(x)),
      se = se_pad, support = support,
      error = paste(c(warnings_seen,
                      "Non-finite or incomplete Cox-refit coefficients."),
                    collapse = " | "),
      warnings = warnings_seen
    ))
  }

  beta_pad[support] <- cc
  se_pad[support] <- ss

  list(
    success = TRUE, empty = FALSE, fit = fit,
    beta = beta_pad, se = se_pad, support = support,
    error = NA_character_, warnings = warnings_seen
  )
}


fit_original_workflow <- function(
    x, time, status,
    nfolds = 10L,
    foldid = NULL,
    standardize = FALSE,
    alpha = 1,
    type_measure = "deviance",
    ties = "efron",
    glmnet_args = list()
) {
  sel <- fit_original_selection(
    x = x, time = time, status = status,
    nfolds = nfolds, foldid = foldid,
    standardize = standardize, alpha = alpha,
    type_measure = type_measure,
    glmnet_args = glmnet_args
  )

  refit <- fit_unpenalized_cox(
    x = x, time = time, status = status,
    support = sel$support, ties = ties
  )
  if (!isTRUE(refit$success)) {
    stop("Original unpenalized Cox refit failed: ", refit$error, call. = FALSE)
  }

  list(selection = sel, refit = refit)
}


make_original_workflow <- function(
    x, time, status,
    support, lambda,
    beta_lasso = NULL,
    ties = "efron",
    selection_rule = "supplied"
) {
  z <- check_survival_input(x, time, status)
  x <- z$x; time <- z$time; status <- z$status

  if (length(lambda) != 1L || is.na(lambda) || !is.numeric(lambda) ||
      !is.finite(lambda) || lambda <= 0) {
    stop("lambda must be one positive finite number.", call. = FALSE)
  }

  support <- sort(unique(as.integer(support)))
  if (anyNA(support) || any(support < 1L | support > ncol(x))) {
    stop("support contains invalid column indices.", call. = FALSE)
  }

  if (is.null(beta_lasso)) {
    beta_lasso <- setNames(rep(NA_real_, ncol(x)), colnames(x))
  } else {
    beta_lasso <- as.numeric(beta_lasso)
    if (length(beta_lasso) != ncol(x)) {
      stop("beta_lasso must have length ncol(x).", call. = FALSE)
    }
    names(beta_lasso) <- colnames(x)
  }

  refit <- fit_unpenalized_cox(
    x = x, time = time, status = status,
    support = support, ties = ties
  )
  if (!isTRUE(refit$success)) {
    stop("Original unpenalized Cox refit failed: ", refit$error, call. = FALSE)
  }

  list(
    selection = list(
      cvfit = NULL,
      lambda_min = as.numeric(lambda),
      beta_lasso = beta_lasso,
      support = support,
      foldid = NULL,
      standardize = NA,
      alpha = NA,
      rule = selection_rule
    ),
    refit = refit
  )
}


fit_original_fixed_lambda <- function(
    x, time, status, lambda,
    standardize = FALSE,
    alpha = 1,
    ties = "efron",
    tol = 1e-10,
    glmnet_args = list()
) {
  z <- check_survival_input(x, time, status)
  x <- z$x; time <- z$time; status <- z$status

  sel <- select_fixed_lambda(
    x = x, time = time, status = status, lambda = lambda,
    standardize = standardize, alpha = alpha,
    tol = tol, glmnet_args = glmnet_args
  )
  if (!isTRUE(sel$success)) {
    stop("Original fixed-lambda Cox-Lasso fit failed: ", sel$error, call. = FALSE)
  }

  refit <- fit_unpenalized_cox(
    x = x, time = time, status = status,
    support = sel$support, ties = ties
  )
  if (!isTRUE(refit$success)) {
    stop("Original unpenalized Cox refit failed: ", refit$error, call. = FALSE)
  }

  sel$lambda_min <- sel$lambda
  sel$rule <- "fixed_lambda"
  list(selection = sel, refit = refit)
}


.validate_original_workflow <- function(original, x) {
  if (!is.list(original) || !is.list(original$selection) ||
      !is.list(original$refit)) {
    stop("original must contain $selection and $refit lists.", call. = FALSE)
  }

  if (is.null(original$selection$support) ||
      is.null(original$selection$lambda_min)) {
    stop("original$selection must contain support and lambda_min.", call. = FALSE)
  }

  support <- sort(unique(as.integer(original$selection$support)))

  if (anyNA(support) || any(support < 1L | support > ncol(x))) {
    stop("original$selection$support is incompatible with x.", call. = FALSE)
  }

  if (!isTRUE(original$refit$success) ||
      length(original$refit$beta) != ncol(x)) {
    stop(
      "original$refit must be a successful p-dimensional refit.",
      call. = FALSE
    )
  }

  original_names <- names(original$refit$beta)

  if (is.null(original_names) ||
      !identical(original_names, colnames(x))) {
    stop(
      "original$refit design-matrix coordinates must match x exactly, ",
      "including column names and order.",
      call. = FALSE
    )
  }

  original$selection$support <- support
  original
}


# -------------------------------------------------------------------------
# 2. Model-based bootstrap generator
# -------------------------------------------------------------------------

fit_reverse_km_censoring <- function(time, status) {
  # Reverse Kaplan-Meier: censoring is treated as the event.
  survival::survfit(survival::Surv(time, 1L - status) ~ 1)
}

reverse_km_discrete_distribution <- function(time, status) {
  # Exact representation used by the archived NKI70 production script:
  # place mass at reverse-KM censoring-event times and retain residual
  # survival mass at Inf.
  km <- fit_reverse_km_censoring(time, status)
  idx <- which(km$n.event > 0)

  if (!length(idx)) {
    return(list(
      fit = km,
      values = Inf,
      prob = 1
    ))
  }

  c_times <- as.numeric(km$time[idx])
  c_surv_after <- as.numeric(km$surv[idx])
  c_surv_before <- c(1, head(c_surv_after, -1L))
  c_mass <- pmax(c_surv_before - c_surv_after, 0)
  c_tail <- max(utils::tail(c_surv_after, 1L), 0)

  values <- c(c_times, Inf)
  prob <- c(c_mass, c_tail)
  prob <- prob / sum(prob)

  list(
    fit = km,
    values = values,
    prob = prob
  )
}

sample_from_step_survival <- function(km, n) {
  # Clean/reference inverse-transform sampler.
  # If the KM curve does not reach the requested probability, return Inf,
  # corresponding to residual mass beyond the observed support.
  tt <- as.numeric(km$time)
  ss <- as.numeric(km$surv)

  if (!length(tt)) return(rep(Inf, n))

  u <- stats::runif(n)
  out <- rep(Inf, n)

  for (i in seq_len(n)) {
    idx <- which(ss <= u[i])[1L]
    if (!is.na(idx)) out[i] <- tt[idx]
  }
  out
}

sample_from_reverse_km_mass <- function(values, prob, n) {
  # Exact archived NKI70 sampling mechanism. Using sample.int() rather than
  # inverse-transform runif() is important for RNG-level reproducibility.
  values[
    sample.int(length(values), size = n, replace = TRUE, prob = prob)
  ]
}

sample_event_times_breslow <- function(basehaz_df, linear_predictor) {
  # Clean/reference event-time sampler.
  hz <- as.numeric(basehaz_df$hazard)
  tt <- as.numeric(basehaz_df$time)

  keep <- is.finite(hz) & is.finite(tt)
  hz <- hz[keep]; tt <- tt[keep]

  if (!length(hz)) return(rep(Inf, length(linear_predictor)))

  hz <- cummax(hz)
  target <- -log(stats::runif(length(linear_predictor))) /
    exp(pmin(pmax(linear_predictor, -700), 700))

  out <- rep(Inf, length(target))
  for (i in seq_along(target)) {
    idx <- which(hz >= target[i])[1L]
    if (!is.na(idx)) out[i] <- tt[idx]
  }
  out
}

sample_event_times_breslow_legacy_nki70 <- function(
    basehaz_df, linear_predictor
) {
  # Mirrors code/production/nki70/01_nki70_modelbased_reselection.R.
  u <- stats::runif(length(linear_predictor))
  target <- -log(u) / exp(linear_predictor)
  vapply(
    target,
    function(z) {
      ii <- which(basehaz_df$hazard >= z)[1L]
      if (is.na(ii)) Inf else basehaz_df$time[ii]
    },
    numeric(1)
  )
}

fit_modelbased_generator <- function(
    x, time, status, original_refit,
    censoring = c("reverse_km"),
    administrative_tau = max(time),
    centered_basehaz = FALSE,
    generator_mode = c("reference", "legacy_nki70")
) {
  censoring <- match.arg(censoring)
  generator_mode <- match.arg(generator_mode)
  x <- as_numeric_matrix(x)

  if (!isTRUE(original_refit$success) || is.null(original_refit$fit)) {
    stop("A successful original Cox refit is required.", call. = FALSE)
  }

  bh <- if (identical(generator_mode, "legacy_nki70")) {
    # Preserve the archived production implementation byte-for-byte at the
    # algorithmic level.
    survival::basehaz(
      original_refit$fit,
      centered = centered_basehaz
    )
  } else {
    # Clean/reference workflow: Efron coefficient refit, but the model-based
    # bootstrap specification uses a Breslow cumulative baseline hazard.
    .reselection_breslow_basehaz(
      original_refit$fit,
      centered = centered_basehaz
    )
  }

  if (identical(generator_mode, "legacy_nki70")) {
    # Archived NKI70 production kept the final cumulative-hazard value at
    # each tied time, then sorted by time.
    keep_bh <- !duplicated(bh$time, fromLast = TRUE)
    bh <- bh[keep_bh, c("hazard", "time"), drop = FALSE]
    bh <- bh[order(bh$time), , drop = FALSE]

    if (!nrow(bh) || any(!is.finite(bh$hazard)) ||
        max(bh$hazard) <= 0) {
      stop(
        "Could not construct a valid Breslow baseline cumulative hazard.",
        call. = FALSE
      )
    }
  }

  # With centered=FALSE, use the ordinary X beta linear predictor.
  lp <- drop(x %*% original_refit$beta)

  if (identical(generator_mode, "legacy_nki70")) {
    censor_dist <- reverse_km_discrete_distribution(time, status)
    censor_fit <- censor_dist$fit
    censor_values <- censor_dist$values
    censor_prob <- censor_dist$prob
  } else {
    censor_fit <- switch(
      censoring,
      reverse_km = fit_reverse_km_censoring(time, status)
    )
    censor_values <- NULL
    censor_prob <- NULL
  }

  list(
    x = x,
    basehaz = bh,
    beta = original_refit$beta,
    linear_predictor = lp,
    censor_fit = censor_fit,
    censor_values = censor_values,
    censor_prob = censor_prob,
    censoring = censoring,
    administrative_tau = administrative_tau,
    generator_mode = generator_mode
  )
}

generate_modelbased_bootstrap <- function(generator) {
  x <- generator$x
  n <- nrow(x)

  if (identical(generator$generator_mode, "legacy_nki70")) {
    # Exact archived NKI70 ordering:
    #   1) n uniforms for event times;
    #   2) sample.int() for reverse-KM censoring;
    #   3) administrative truncation at tau for every subject.
    t_event <- sample_event_times_breslow_legacy_nki70(
      generator$basehaz,
      generator$linear_predictor
    )

    t_cens <- sample_from_reverse_km_mass(
      generator$censor_values,
      generator$censor_prob,
      n
    )

    w <- pmin(t_event, t_cens, generator$administrative_tau)
    d <- as.integer(
      t_event <= t_cens &
        t_event <= generator$administrative_tau
    )

    if (any(!is.finite(w)) || any(w < 0)) {
      stop("Non-finite generated follow-up time.", call. = FALSE)
    }

    return(list(
      x = x,
      time = w,
      status = d,
      event_time = t_event,
      censor_time = t_cens
    ))
  }

  t_event <- sample_event_times_breslow(
    generator$basehaz,
    generator$linear_predictor
  )

  t_cens <- switch(
    generator$censoring,
    reverse_km = sample_from_step_survival(generator$censor_fit, n)
  )

  # Clean/reference behavior retained unchanged for future analyses.
  both_inf <- is.infinite(t_event) & is.infinite(t_cens)
  if (any(both_inf)) {
    t_cens[both_inf] <- generator$administrative_tau
  }

  w <- pmin(t_event, t_cens)
  d <- as.integer(t_event <= t_cens)

  bad <- !is.finite(w)
  if (any(bad)) {
    w[bad] <- generator$administrative_tau
    d[bad] <- 0L
  }

  list(
    x = x,
    time = w,
    status = d,
    event_time = t_event,
    censor_time = t_cens
  )
}


# -------------------------------------------------------------------------
# 3. Bootstrap selection rules
# -------------------------------------------------------------------------

select_fixed_lambda <- function(
    x, time, status, lambda,
    standardize = FALSE, alpha = 1,
    tol = 1e-10, glmnet_args = list()
) {
  y <- survival::Surv(time, status)

  args <- c(
    list(
      x = x, y = y, family = "cox",
      alpha = alpha,
      lambda = lambda,
      standardize = standardize
    ),
    glmnet_args
  )

  args <- .reselection_prepare_cox_glmnet_args(
    args,
    cv = FALSE
  )
  fit <- tryCatch(
    do.call(glmnet::glmnet, args),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(success = FALSE, support = integer(), beta_lasso = NULL,
                lambda = lambda, fit = NULL, error = conditionMessage(fit)))
  }

  b <- safe_coef_glmnet(fit, s = lambda)
  names(b) <- colnames(x)

  list(
    success = TRUE,
    support = support_from_beta(b, tol = tol),
    beta_lasso = b,
    lambda = lambda,
    fit = fit,
    error = NA_character_
  )
}

select_cv_rule <- function(
    x, time, status,
    rule = c("min", "1se"),
    nfolds = 10L,
    foldid = NULL,
    standardize = FALSE,
    alpha = 1,
    type_measure = "deviance",
    tol = 1e-10,
    glmnet_args = list()
) {
  rule <- match.arg(rule)
  if (is.null(foldid)) foldid <- balanced_foldid(nrow(x), nfolds)
  y <- survival::Surv(time, status)

  args <- c(
    list(
      x = x, y = y, family = "cox",
      alpha = alpha,
      foldid = foldid,
      type.measure = type_measure,
      standardize = standardize
    ),
    glmnet_args
  )

  args <- .reselection_prepare_cox_glmnet_args(
    args,
    cv = TRUE
  )
  fit <- tryCatch(
    do.call(glmnet::cv.glmnet, args),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(success = FALSE, support = integer(), beta_lasso = NULL,
                lambda = NA_real_, fit = NULL, foldid = foldid,
                rule = rule, error = conditionMessage(fit)))
  }

  lambda_name <- if (identical(rule, "min")) "lambda.min" else "lambda.1se"
  lambda_value <- unname(fit[[lambda_name]])
  coef_s <- if (identical(rule, "min")) "lambda.min" else "lambda.1se"
  b <- safe_coef_glmnet(fit, s = coef_s)
  names(b) <- colnames(x)

  list(
    success = TRUE,
    support = support_from_beta(b, tol = tol),
    beta_lasso = b,
    lambda = lambda_value,
    fit = fit,
    foldid = foldid,
    rule = rule,
    error = NA_character_
  )
}

select_cv_min <- function(...) {
  select_cv_rule(..., rule = "min")
}

select_cv_1se <- function(...) {
  select_cv_rule(..., rule = "1se")
}

select_aic_path <- function(
    x, time, status,
    standardize = FALSE,
    alpha = 1,
    tol = 1e-10,
    glmnet_args = list()
) {
  y <- survival::Surv(time, status)

  args <- c(
    list(
      x = x, y = y, family = "cox",
      alpha = alpha,
      standardize = standardize
    ),
    glmnet_args
  )

  args <- .reselection_prepare_cox_glmnet_args(
    args,
    cv = FALSE
  )
  fit <- tryCatch(
    do.call(glmnet::glmnet, args),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(success = FALSE, support = integer(), beta_lasso = NULL,
                lambda = NA_real_, fit = NULL, aic_path = NULL,
                error = conditionMessage(fit)))
  }

  # Manuscript rule:
  #   AIC(lambda) = D(lambda) + 2 * d(lambda),
  # where D is the in-sample Cox deviance from glmnet and d is the number
  # of nonzero Cox-Lasso coefficients.
  dev <- as.numeric(fit$nulldev) * (1 - as.numeric(fit$dev.ratio))
  df <- as.numeric(fit$df)
  aic <- dev + 2 * df

  good <- is.finite(aic) & is.finite(fit$lambda)
  if (!any(good)) {
    return(list(success = FALSE, support = integer(), beta_lasso = NULL,
                lambda = NA_real_, fit = fit,
                aic_path = data.frame(lambda = fit$lambda, deviance = dev, df = df, AIC = aic),
                error = "No finite AIC value on glmnet path."))
  }

  idx_candidates <- which(good)
  idx <- idx_candidates[which.min(aic[good])]
  lambda_aic <- fit$lambda[idx]

  bmat <- as.matrix(stats::coef(fit, s = lambda_aic))
  b <- drop(bmat[, 1L])
  names(b) <- colnames(x)

  list(
    success = TRUE,
    support = support_from_beta(b, tol = tol),
    beta_lasso = b,
    lambda = unname(lambda_aic),
    fit = fit,
    path_index = idx,
    aic_path = data.frame(
      lambda = fit$lambda,
      deviance = dev,
      df = df,
      AIC = aic
    ),
    error = NA_character_
  )
}


# -------------------------------------------------------------------------
# 4. Run one bootstrap method on one generated dataset
# -------------------------------------------------------------------------

analyse_bootstrap_dataset <- function(
    boot,
    method = c(
      "fixed_support", "reselect_fixed_lambda",
      "reselect_cv_min", "reselect_cv_1se", "reselect_aic"
    ),
    original,
    nfolds = 10L,
    cv_foldid = NULL,
    standardize = FALSE,
    alpha = 1,
    type_measure = "deviance",
    ties = "efron",
    tol = 1e-10,
    glmnet_args = list(),
    refit_mode = c("reference", "legacy_nki70")
) {
  method <- .normalize_reselection_methods(method)
  refit_mode <- match.arg(refit_mode)
  if (length(method) != 1L) {
    stop("analyse_bootstrap_dataset() requires exactly one method.", call. = FALSE)
  }
  x <- boot$x
  time <- boot$time
  status <- boot$status

  if (method == "fixed_support") {
    sel <- list(
      success = TRUE,
      support = original$selection$support,
      lambda = original$selection$lambda_min,
      error = NA_character_
    )
  } else if (method == "reselect_fixed_lambda") {
    sel <- select_fixed_lambda(
      x, time, status,
      lambda = original$selection$lambda_min,
      standardize = standardize, alpha = alpha,
      tol = tol, glmnet_args = glmnet_args
    )
  } else if (method == "reselect_cv_min") {
    sel <- select_cv_min(
      x, time, status,
      nfolds = nfolds,
      foldid = cv_foldid,
      standardize = standardize, alpha = alpha,
      type_measure = type_measure,
      tol = tol, glmnet_args = glmnet_args
    )
  } else if (method == "reselect_cv_1se") {
    sel <- select_cv_1se(
      x, time, status,
      nfolds = nfolds,
      foldid = cv_foldid,
      standardize = standardize, alpha = alpha,
      type_measure = type_measure,
      tol = tol, glmnet_args = glmnet_args
    )
  } else {
    sel <- select_aic_path(
      x, time, status,
      standardize = standardize, alpha = alpha,
      tol = tol, glmnet_args = glmnet_args
    )
  }

  if (!isTRUE(sel$success)) {
    p <- ncol(x)
    return(list(
      method = method,
      selection_success = FALSE,
      refit_success = FALSE,
      support = integer(),
      beta = setNames(rep(NA_real_, p), colnames(x)),
      se = setNames(rep(NA_real_, p), colnames(x)),
      lambda = sel$lambda,
      selection_error = sel$error,
      refit_error = NA_character_,
      aic_path = sel$aic_path %||% NULL
    ))
  }

  refit <- if (identical(refit_mode, "legacy_nki70")) {
    fit_unpenalized_cox_legacy_nki70(
      x, time, status, support = sel$support
    )
  } else {
    fit_unpenalized_cox(
      x, time, status, support = sel$support, ties = ties
    )
  }

  list(
    method = method,
    selection_success = TRUE,
    refit_success = isTRUE(refit$success),
    support = sel$support,
    beta = refit$beta,
    se = refit$se,
    lambda = sel$lambda,
    selection_error = NA_character_,
    refit_error = refit$error,
    aic_path = sel$aic_path %||% NULL
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x


# -------------------------------------------------------------------------
# 5. Full bootstrap: same generated dataset across all methods
# -------------------------------------------------------------------------

run_reselection_bootstrap <- function(
    x, time, status,
    B = 500L,
    methods = c(
      "fixed_support",
      "reselect_fixed_lambda",
      "reselect_cv_min",
      "reselect_cv_1se",
      "reselect_aic"
    ),
    conf_level = 0.95,
    seed = NULL,
    nfolds = 10L,
    original = NULL,
    original_support = NULL,
    original_lambda = NULL,
    original_foldid = NULL,
    standardize = TRUE,
    alpha = 1,
    type_measure = "deviance",
    ties = "efron",
    censoring = "reverse_km",
    rng_mode = c("independent_streams", "legacy_nki70"),
    glmnet_args = list(),
    progress = interactive()
) {
  check_reselection_packages()
  z <- check_survival_input(x, time, status)
  x <- z$x; time <- z$time; status <- z$status

  if (length(B) != 1L || is.na(B) || !is.numeric(B) || !is.finite(B) ||
      B < 1 || B > .Machine$integer.max || B != floor(B)) {
    stop("B must be one positive integer within the R integer range.", call. = FALSE)
  }
  B <- as.integer(B)

  if (length(conf_level) != 1L || is.na(conf_level) ||
      !is.numeric(conf_level) || !is.finite(conf_level) ||
      conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be strictly between 0 and 1.", call. = FALSE)
  }

  if (length(alpha) == 1L &&
      is.numeric(alpha) &&
      !is.na(alpha) &&
      is.finite(alpha) &&
      abs(alpha - 1) > 1e-12) {
    warning(
      "coxlasso_boot_reselect() is proceeding with `alpha = ",
      format(alpha, digits = 8),
      "`. selSurv's inferential validation currently concerns Lasso ",
      "screening (`alpha = 1`). Other mixing parameters are available ",
      "for fitting flexibility but are not separately validated ",
      "inferentially.",
      call. = FALSE
    )
  }

  methods <- .normalize_reselection_methods(methods)
  rng_mode <- match.arg(rng_mode)
  seed <- .reselection_resolve_root_seed(
    seed,
    rng_mode
  )
  ties <- .reselection_resolve_refit_ties(
    ties,
    rng_mode
  )

  if (!is.null(original)) {
    if (!is.null(original_support) || !is.null(original_lambda)) {
      stop(
        "Supply either original or original_support/original_lambda, not both.",
        call. = FALSE
      )
    }
    original <- .validate_original_workflow(original, x)
  } else if (!is.null(original_support) || !is.null(original_lambda)) {
    if (is.null(original_support) || is.null(original_lambda)) {
      stop(
        "original_support and original_lambda must be supplied together.",
        call. = FALSE
      )
    }
    original <- make_original_workflow(
      x, time, status,
      support = original_support,
      lambda = original_lambda,
      ties = ties,
      selection_rule = "supplied_support_and_lambda"
    )
  } else {
    # Reserve one deterministic stream for the observed-data CV split.
    if (is.null(original_foldid)) {
      original_foldid <- .with_reselection_seed(
        .reselection_seed(seed, replicate = 0L, stream = 11L),
        balanced_foldid(nrow(x), nfolds)
      )
    }

    original <- fit_original_workflow(
      x, time, status,
      nfolds = nfolds,
      foldid = original_foldid,
      standardize = standardize,
      alpha = alpha,
      type_measure = type_measure,
      ties = ties,
      glmnet_args = glmnet_args
    )
  }

  generator <- fit_modelbased_generator(
    x, time, status,
    original_refit = original$refit,
    censoring = censoring,
    generator_mode = if (identical(rng_mode, "legacy_nki70")) {
      "legacy_nki70"
    } else {
      "reference"
    }
  )

  p <- ncol(x)
  P <- length(methods)

  beta_store <- lapply(methods, function(m) {
    matrix(NA_real_, nrow = B, ncol = p,
           dimnames = list(NULL, colnames(x)))
  })
  names(beta_store) <- methods

  se_store <- lapply(methods, function(m) {
    matrix(NA_real_, nrow = B, ncol = p,
           dimnames = list(NULL, colnames(x)))
  })
  names(se_store) <- methods

  selected_store <- lapply(methods, function(m) {
    matrix(FALSE, nrow = B, ncol = p,
           dimnames = list(NULL, colnames(x)))
  })
  names(selected_store) <- methods

  lambda_store <- lapply(methods, function(m) rep(NA_real_, B))
  names(lambda_store) <- methods

  support_store <- lapply(methods, function(m) vector("list", B))
  names(support_store) <- methods

  selection_success <- lapply(methods, function(m) rep(FALSE, B))
  names(selection_success) <- methods

  refit_success <- lapply(methods, function(m) rep(FALSE, B))
  names(refit_success) <- methods

  errors <- vector("list", B * P)
  err_i <- 0L

  store_answer <- function(ans, m, b) {
    support_store[[m]][[b]] <<- ans$support
    if (length(ans$support)) {
      selected_store[[m]][b, ans$support] <<- TRUE
    }
    lambda_store[[m]][b] <<- ans$lambda
    selection_success[[m]][b] <<- ans$selection_success
    refit_success[[m]][b] <<- ans$refit_success

    if (isTRUE(ans$refit_success)) {
      beta_store[[m]][b, ] <<- ans$beta
      se_store[[m]][b, ] <<- ans$se
    }

    if (!isTRUE(ans$selection_success) || !isTRUE(ans$refit_success)) {
      err_i <<- err_i + 1L
      errors[[err_i]] <<- data.frame(
        b = b,
        method = m,
        selection_success = ans$selection_success,
        refit_success = ans$refit_success,
        selection_error = ans$selection_error %||% NA_character_,
        refit_error = ans$refit_error %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }
    invisible(NULL)
  }

  for (b in seq_len(B)) {
    if (identical(rng_mode, "independent_streams")) {
      # The generated dataset and CV folds have separate replicate-specific
      # streams. Results for a named method are therefore invariant to adding
      # or reordering other methods.
      boot <- .with_reselection_seed(
        .reselection_seed(seed, replicate = b, stream = 101L),
        generate_modelbased_bootstrap(generator)
      )

      cv_foldid <- if (any(c("reselect_cv_min", "reselect_cv_1se") %in% methods)) {
        .with_reselection_seed(
          .reselection_seed(seed, replicate = b, stream = 211L),
          balanced_foldid(nrow(x), nfolds)
        )
      } else {
        NULL
      }

      for (m in methods) {
        ans <- analyse_bootstrap_dataset(
          boot = boot, method = m, original = original,
          nfolds = nfolds, cv_foldid = cv_foldid,
          standardize = standardize, alpha = alpha,
          type_measure = type_measure, ties = ties,
          glmnet_args = glmnet_args, refit_mode = "reference"
        )
        store_answer(ans, m, b)
      }
    } else {
      # Reproduction mode for the archived final NKI70 run. The entire
      # replicate uses set.seed(seed + 100003*b), then generation and methods
      # consume that stream in the supplied order. Keep the historical method
      # order when exact parity with archived outputs is required.
      replicate_answers <- .with_reselection_seed(
        .reselection_seed(seed, replicate = b, stream = 0L),
        {
          boot <- generate_modelbased_bootstrap(generator)
          lapply(methods, function(m) {
            analyse_bootstrap_dataset(
              boot = boot, method = m, original = original,
              nfolds = nfolds, cv_foldid = NULL,
              standardize = standardize, alpha = alpha,
              type_measure = type_measure, ties = ties,
              # Archived NKI70 production defined selected coordinates by
              # abs(beta) > 1e-8 and used glmnet solver controls
              # thresh = 1e-9, maxit = 100000 for fixed-lambda, CV-min,
              # and AIC selection.  Force those settings only in the exact
              # legacy reproduction path; the clean reference defaults are
              # unchanged.
              tol = 1e-8,
              glmnet_args = utils::modifyList(
                glmnet_args,
                list(thresh = 1e-9, maxit = 100000L)
              ),
              refit_mode = "legacy_nki70"
            )
          })
        }
      )
      names(replicate_answers) <- methods
      for (m in methods) store_answer(replicate_answers[[m]], m, b)
    }

    if (isTRUE(progress) && (b %% max(1L, floor(B / 20L)) == 0L || b == B)) {
      message(sprintf("bootstrap %d/%d", b, B))
    }
  }

  errors <- errors[seq_len(err_i)]
  errors <- if (length(errors)) do.call(rbind, errors) else data.frame()

  out <- list(
    call = match.call(),
    x_colnames = colnames(x),
    B = B,
    conf_level = conf_level,
    seed = seed,
    methods = methods,
    original = original,
    generator = generator,
    beta = beta_store,
    se = se_store,
    selected = selected_store,
    lambda = lambda_store,
    support = support_store,
    selection_success = selection_success,
    refit_success = refit_success,
    errors = errors,
    settings = list(
      root_seed = seed,
      nfolds = nfolds,
      original_foldid = original_foldid,
      original_lambda = original$selection$lambda_min,
      original_support = original$selection$support,
      standardize = standardize,
      alpha = alpha,
      type_measure = type_measure,
      ties = ties,
      censoring = censoring,
      rng_mode = rng_mode
    )
  )
  class(out) <- "cox_reselection_bootstrap"
  out
}


# -------------------------------------------------------------------------
# 6. Interval and stability summaries
# -------------------------------------------------------------------------

summarise_reselection <- function(object, conf_level = object$conf_level) {
  stopifnot(inherits(object, "cox_reselection_bootstrap"))

  alpha <- 1 - conf_level
  qlo <- alpha / 2
  qhi <- 1 - alpha / 2
  zcrit <- stats::qnorm(qhi)

  bhat <- object$original$refit$beta
  sehat <- object$original$refit$se
  Mhat <- object$original$selection$support
  p <- length(bhat)
  vars <- names(bhat)

  rows <- list()
  k <- 0L

  # We report coefficient intervals for originally selected coordinates.
  for (m in object$methods) {
    Bmat <- object$beta[[m]]
    Smat <- object$se[[m]]
    Sel <- object$selected[[m]]
    ok <- object$refit_success[[m]]

    sel_success <- object$selection_success[[m]]
    support_jacc <- rep(NA_real_, object$B)
    support_same <- rep(NA, object$B)
    support_size <- rep(NA_integer_, object$B)

    if (any(sel_success)) {
      support_jacc[sel_success] <- vapply(
        object$support[[m]][sel_success],
        jaccard_support,
        numeric(1),
        b = Mhat
      )
      support_same[sel_success] <- vapply(
        object$support[[m]][sel_success],
        same_support,
        logical(1),
        b = Mhat
      )
      support_size[sel_success] <- lengths(object$support[[m]][sel_success])
    }

    for (j in Mhat) {
      k <- k + 1L

      xb <- Bmat[ok, j]
      xb_finite <- xb[is.finite(xb)]
      nboot_eff <- length(xb_finite)

      if (nboot_eff >= 2L) {
        qq <- stats::quantile(
          xb_finite,
          probs = c(qlo, qhi),
          names = FALSE,
          type = 7,
          na.rm = TRUE
        )
        perc_lo <- qq[1L]
        perc_hi <- qq[2L]
        basic_lo <- 2 * bhat[j] - qq[2L]
        basic_hi <- 2 * bhat[j] - qq[1L]
        boot_sd <- stats::sd(xb_finite)
        bootse_lo <- bhat[j] - zcrit * boot_sd
        bootse_hi <- bhat[j] + zcrit * boot_sd
      } else {
        perc_lo <- perc_hi <- basic_lo <- basic_hi <- NA_real_
        boot_sd <- bootse_lo <- bootse_hi <- NA_real_
      }

      # Studentized sensitivity: only bootstrap refits in which j is selected
      # and a finite refit-specific SE exists.
      stud_ok <- ok & Sel[, j] &
        is.finite(Bmat[, j]) &
        is.finite(Smat[, j]) &
        Smat[, j] > 0

      if (sum(stud_ok) >= 2L && is.finite(sehat[j]) && sehat[j] > 0) {
        tstar <- (Bmat[stud_ok, j] - bhat[j]) / Smat[stud_ok, j]
        tq <- stats::quantile(
          tstar, probs = c(qlo, qhi),
          names = FALSE, type = 7, na.rm = TRUE
        )
        stud_lo <- bhat[j] - tq[2L] * sehat[j]
        stud_hi <- bhat[j] - tq[1L] * sehat[j]
      } else {
        stud_lo <- stud_hi <- NA_real_
      }

      naive_lo <- bhat[j] - zcrit * sehat[j]
      naive_hi <- bhat[j] + zcrit * sehat[j]

      k_all <- length(object$refit_success[[m]])
      n_success <- sum(ok)

      rows[[k]] <- data.frame(
        method = m,
        variable = vars[j],
        j = j,
        estimate = bhat[j],
        original_se = sehat[j],

        naive_lower = naive_lo,
        naive_upper = naive_hi,

        percentile_lower = perc_lo,
        percentile_upper = perc_hi,

        basic_lower = basic_lo,
        basic_upper = basic_hi,

        bootstrap_se = boot_sd,
        bootstrap_se_lower = bootse_lo,
        bootstrap_se_upper = bootse_hi,

        studentized_lower = stud_lo,
        studentized_upper = stud_hi,

        B_requested = object$B,
        B_refit_success = n_success,
        refit_success_rate = n_success / k_all,

        # Selection diagnostics are conditional on a successful selection fit.
        # Failed selection fits are not silently counted as empty supports.
        reselection_frequency =
          if (any(sel_success)) mean(Sel[sel_success, j]) else NA_real_,

        # Zero coefficient mass used in the bootstrap coefficient distribution
        # among successful refits, matching interval construction.
        boot_zero_frequency =
          if (n_success) mean(!Sel[ok, j]) else NA_real_,
        # Backward-compatible names retained for downstream manuscript scripts.
        reselection_frequency_all =
          if (any(sel_success)) mean(Sel[sel_success, j]) else NA_real_,
        boot_zero_frequency_successful =
          if (n_success) mean(!Sel[ok, j]) else NA_real_,

        exact_support_reproduction =
          if (any(sel_success)) mean(support_same[sel_success]) else NA_real_,
        mean_jaccard =
          if (any(sel_success)) mean(support_jacc[sel_success]) else NA_real_,
        mean_support_size =
          if (any(sel_success)) mean(support_size[sel_success]) else NA_real_,
        sd_support_size =
          if (sum(sel_success) >= 2L) stats::sd(support_size[sel_success]) else NA_real_,
        mean_lambda = {
          lam <- object$lambda[[m]][sel_success]
          if (any(is.finite(lam))) mean(lam, na.rm = TRUE) else NA_real_
        },
        sd_lambda = {
          lam <- object$lambda[[m]][sel_success]
          if (sum(is.finite(lam)) >= 2L) stats::sd(lam, na.rm = TRUE) else NA_real_
        },

        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, rows)
}


method_level_diagnostics <- function(object) {
  stopifnot(inherits(object, "cox_reselection_bootstrap"))
  Mhat <- object$original$selection$support

  do.call(rbind, lapply(object$methods, function(m) {
    supports <- object$support[[m]]
    sel_ok <- object$selection_success[[m]]
    refit_ok <- object$refit_success[[m]]

    if (any(sel_ok)) {
      successful_supports <- supports[sel_ok]
      exact <- mean(vapply(
        successful_supports, same_support, logical(1), b = Mhat
      ))
      jacc <- mean(vapply(
        successful_supports, jaccard_support, numeric(1), b = Mhat
      ))
      sizes <- lengths(successful_supports)
      mean_size <- mean(sizes)
      sd_size <- stats::sd(sizes)
      lambda_ok <- object$lambda[[m]][sel_ok]
      mean_lambda <- if (any(is.finite(lambda_ok))) mean(lambda_ok, na.rm = TRUE) else NA_real_
      sd_lambda <- if (sum(is.finite(lambda_ok)) >= 2L) stats::sd(lambda_ok, na.rm = TRUE) else NA_real_
    } else {
      exact <- jacc <- mean_size <- sd_size <- mean_lambda <- sd_lambda <- NA_real_
    }

    data.frame(
      method = m,
      exact_support_reproduction = exact,
      mean_jaccard = jacc,
      mean_jaccard_original = jacc,
      mean_n_selected = mean_size,
      mean_n_selected_boot = mean_size,
      sd_n_selected = sd_size,
      sd_n_selected_boot = sd_size,
      selection_success_rate = mean(sel_ok),
      refit_success_rate = mean(refit_ok),
      bootstrap_success_rate = mean(refit_ok),
      refit_success_given_selection =
        if (any(sel_ok)) mean(refit_ok[sel_ok]) else NA_real_,
      mean_lambda = mean_lambda,
      sd_lambda = sd_lambda,
      stringsAsFactors = FALSE
    )
  }))
}


# -------------------------------------------------------------------------
# 7. Simulation target helpers
# -------------------------------------------------------------------------

fit_support_population_target <- function(
    x_ref, time_ref, status_ref, support,
    ties = "efron"
) {
  fit <- fit_unpenalized_cox(
    x_ref, time_ref, status_ref,
    support = support, ties = ties
  )
  if (!isTRUE(fit$success)) {
    stop("Reference-population Cox refit failed for support ",
         support_key(support), ": ", fit$error, call. = FALSE)
  }
  fit$beta
}

make_support_target_cache <- function(
    x_ref, time_ref, status_ref,
    ties = "efron"
) {
  cache <- new.env(parent = emptyenv())

  function(support) {
    key <- support_key(support)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(
        key,
        fit_support_population_target(
          x_ref, time_ref, status_ref,
          support = support, ties = ties
        ),
        envir = cache
      )
    }
    get(key, envir = cache, inherits = FALSE)
  }
}


evaluate_against_target <- function(
    interval_summary,
    beta_target,
    beta_true = NULL,
    interval = c("bootstrap_se", "basic", "percentile", "studentized")
) {
  interval <- match.arg(interval)
  lo_col <- paste0(interval, "_lower")
  hi_col <- paste0(interval, "_upper")

  idx <- interval_summary$j
  target <- beta_target[idx]

  ans <- interval_summary
  ans$target <- target
  ans$covered_target <-
    is.finite(ans[[lo_col]]) &
    is.finite(ans[[hi_col]]) &
    target >= ans[[lo_col]] &
    target <= ans[[hi_col]]
  ans$interval_width <- ans[[hi_col]] - ans[[lo_col]]
  ans$reject0 <-
    is.finite(ans[[lo_col]]) &
    is.finite(ans[[hi_col]]) &
    (ans[[lo_col]] > 0 | ans[[hi_col]] < 0)

  if (!is.null(beta_true)) {
    bt <- beta_true[idx]
    ans$beta_true <- bt
    ans$covered_beta_true <-
      is.finite(ans[[lo_col]]) &
      is.finite(ans[[hi_col]]) &
      bt >= ans[[lo_col]] &
      bt <= ans[[hi_col]]
    ans$abs_selected_support_target_deviation <- abs(target - bt)
  }

  ans
}


legacy_performance_summary <- function(evaluated_df) {
  # One-Monte-Carlo-replicate building block.
  # Across MC replicates, average these columns by method/scenario as needed.
  do.call(rbind, lapply(split(evaluated_df, evaluated_df$method), function(d) {
    data.frame(
      method = d$method[1L],
      target2_coverage_cond_selected = mean(d$covered_target, na.rm = TRUE),
      mean_width_cond_selected = mean(d$interval_width, na.rm = TRUE),
      mean_reselect_freq = mean(d$reselection_frequency, na.rm = TRUE),
      mean_boot_zero_freq = mean(d$boot_zero_frequency, na.rm = TRUE),
      mean_abs_target_drift =
        if ("abs_selected_support_target_deviation" %in% names(d))
          mean(d$abs_selected_support_target_deviation, na.rm = TRUE)
        else NA_real_,
      reject0_rate_cond_selected = mean(d$reject0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}


#' @export
print.cox_reselection_bootstrap <- function(x, ...) {
  cat("Cox-Lasso model-based bootstrap with re-selection\n")
  cat("B:", x$B, "\n")
  cat("Original lambda_min:", signif(x$original$selection$lambda_min, 5), "\n")
  cat("Original support size:", length(x$original$selection$support), "\n")
  cat("Methods:", paste(x$methods, collapse = ", "), "\n")
  if (!is.null(x$settings$rng_mode)) cat("RNG mode:", x$settings$rng_mode, "\n")
  cat("\n")
  print(method_level_diagnostics(x), row.names = FALSE)
  invisible(x)
}
