# Prepare the analysis object ---------------------------------------------

.selsurv_prepare_boot_object <- function(object, inference) {
  if (identical(inference, "debiased")) {
    return(.selsurv_as_debiased(object))
  }

  if (inherits(object, c("selSurv_debias", "selsurv_debiased"))) {
    original_fit <- attr(object, "selsurv_fit", exact = TRUE)

    if (is.null(original_fit) && is.list(object)) {
      original_fit <- object$fit
    }

    if (is.null(original_fit)) {
      stop(
        "Conditional bootstrap requires the original fitted object. ",
        "Pass the result of `coxlasso_fit()` rather than a standalone ",
        "debiased object.",
        call. = FALSE
      )
    }

    return(original_fit)
  }

  object
}


.selsurv_as_debiased <- function(object) {
  if (inherits(object, c("selSurv_debias", "selsurv_debiased"))) {
    return(object)
  }

  debias_fun <- .selsurv_find_backend("coxlasso_debias")
  result <- debias_fun(object)

  # v12 uses class `selSurv_debias`; keep the established class and attach
  # the original fit for round-tripping through the unified interface.
  attr(result, "selsurv_fit") <- object
  result
}


# Dispatch to existing method engines ------------------------------------

.selsurv_run_boot_engine <- function(
    object,
    inference,
    type,
    B,
    conf_level,
    seed,
    dots
) {
  registry <- .selsurv_boot_registry()
  backend_name <- unname(registry[[inference]][[type]])
  backend <- .selsurv_find_backend(backend_name)

  # The unified `conditional` family is deliberately fixed-support.
  # Procedure-level re-selection remains on `coxlasso_boot_parametric()`.
  if (identical(inference, "conditional") && identical(type, "model_based")) {
    if ("workflow" %in% names(dots) && !identical(dots$workflow, "fixed_support")) {
      stop(
        "`coxlasso_boot(inference = \"conditional\")` uses the fixed-support ",
        "model-based bootstrap. Call `coxlasso_boot_parametric()` directly for ",
        "re-selection workflows.",
        call. = FALSE
      )
    }
    dots$workflow <- "fixed_support"
  }

  arguments <- c(
    list(
      object = object,
      B = B,
      conf_level = conf_level,
      seed = seed
    ),
    dots
  )

  result <- .selsurv_call_backend(
    backend = backend,
    arguments = arguments,
    seed = seed
  )

  structure(
    list(
      inference = inference,
      type = type,
      result = result,
      settings = list(
        B = B,
        conf_level = conf_level,
        seed = seed,
        backend = backend_name
      )
    ),
    class = c("selsurv_boot_method", "selsurv_result")
  )
}


.selsurv_find_backend <- function(name) {
  if (!exists(name, mode = "function", inherits = TRUE)) {
    stop(
      "The backend function ",
      shQuote(name),
      " is not available. Ensure that the existing selSurv method ",
      "implementation is loaded before using the unified interface.",
      call. = FALSE
    )
  }

  get(name, mode = "function", inherits = TRUE)
}


.selsurv_call_backend <- function(backend, arguments, seed = NULL) {
  backend_formals <- names(formals(backend))
  supports_dots <- "..." %in% backend_formals

  # The historical functions may use different names for their first object
  # argument. Map the unified `object` argument to that first formal.
  object_position <- which(names(arguments) == "object")
  first_formal <- setdiff(backend_formals, "...")[[1L]]

  if (!first_formal %in% names(arguments)) {
    names(arguments)[object_position] <- first_formal
  }

  # Map the common confidence-level parameter to older conventions.
  if (
    "conf_level" %in% names(arguments) &&
      !"conf_level" %in% backend_formals
  ) {
    if ("level" %in% backend_formals) {
      arguments$level <- arguments$conf_level
    } else if ("alpha" %in% backend_formals) {
      arguments$alpha <- 1 - arguments$conf_level
    }

    arguments$conf_level <- NULL
  }

  # Seed handling is performed outside the backend unless it explicitly
  # provides a seed argument.
  if (!"seed" %in% backend_formals) {
    arguments$seed <- NULL
  }

  if (!supports_dots) {
    arguments <- arguments[
      intersect(names(arguments), backend_formals)
    ]
  }

  .selsurv_with_seed(
    seed,
    do.call(backend, arguments)
  )
}


.selsurv_with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

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
