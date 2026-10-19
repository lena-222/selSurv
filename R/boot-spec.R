# Internal registry --------------------------------------------------------

# Keep the relationship between inferential targets and bootstrap engines in
# one place. New engines should be added here instead of being hard-coded
# throughout the public API.
.selsurv_boot_registry <- function() {
  list(
    debiased = c(
      wild = "coxlasso_wildboot",
      efron = "coxlasso_efron_events"
    ),
    conditional = c(
      model_based = "coxlasso_boot_parametric"
    )
  )
}


.selsurv_boot_defaults <- function() {
  c(
    debiased = "wild",
    conditional = "model_based"
  )
}


# Resolve user specification ---------------------------------------------

.selsurv_resolve_boot_spec <- function(inference = NULL, type = NULL) {
  registry <- .selsurv_boot_registry()
  defaults <- .selsurv_boot_defaults()

  valid_inference <- names(registry)
  all_types <- unlist(lapply(registry, names), use.names = FALSE)

  if (!is.null(inference)) {
    .selsurv_assert_character_scalar(
      inference,
      argument = "inference",
      choices = valid_inference
    )
  }

  if (!is.null(type)) {
    if (!is.character(type) || length(type) < 1L || anyNA(type)) {
      stop(
        "`type` must be a non-empty character vector.",
        call. = FALSE
      )
    }

    if (any(!nzchar(type))) {
      stop("`type` must not contain empty strings.", call. = FALSE)
    }

    if (anyDuplicated(type)) {
      type <- unique(type)
    }
  }

  # Fully implicit default: provide the shortest useful workflow.
  if (is.null(inference) && is.null(type)) {
    return(list(
      inference = "debiased",
      type = unname(defaults[["debiased"]])
    ))
  }

  # An explicit inferential family selects its documented default engine.
  if (!is.null(inference) && is.null(type)) {
    return(list(
      inference = inference,
      type = unname(defaults[[inference]])
    ))
  }

  # `type = "all"` is only meaningful when the inferential family is known.
  if (is.null(inference) && identical(type, "all")) {
    stop(
      "`type = \"all\"` requires an explicit `inference` value.",
      call. = FALSE
    )
  }

  if (!is.null(inference) && identical(type, "all")) {
    return(list(
      inference = inference,
      type = names(registry[[inference]])
    ))
  }

  unknown <- setdiff(type, all_types)
  if (length(unknown) > 0L) {
    stop(
      "Unknown bootstrap type",
      if (length(unknown) > 1L) "s" else "",
      ": ",
      paste(shQuote(unknown), collapse = ", "),
      ". Available types are ",
      paste(shQuote(all_types), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  inferred_targets <- unique(vapply(
    type,
    .selsurv_inference_for_type,
    character(1L),
    registry = registry
  ))

  # A type-only call is concise and unambiguous when all types belong to the
  # same inferential family.
  if (is.null(inference)) {
    if (length(inferred_targets) != 1L) {
      mapping <- paste0(
        "- ", shQuote(type), ": ",
        vapply(
          type,
          .selsurv_inference_for_type,
          character(1L),
          registry = registry
        ),
        collapse = "\n"
      )

      stop(
        "The supplied bootstrap types correspond to different ",
        "inference targets:\n",
        mapping,
        "\nUse separate calls or types from one inference family.",
        call. = FALSE
      )
    }

    return(list(
      inference = inferred_targets,
      type = type
    ))
  }

  # Explicit but incompatible arguments must fail. Silently changing either
  # argument would hide a scientifically meaningful user decision.
  invalid <- setdiff(type, names(registry[[inference]]))
  if (length(invalid) > 0L) {
    stop(
      "Bootstrap type",
      if (length(invalid) > 1L) "s " else " ",
      paste(shQuote(invalid), collapse = ", "),
      " ",
      if (length(invalid) > 1L) "are" else "is",
      " not available for `inference = ",
      shQuote(inference),
      "`.\nAvailable types: ",
      paste(shQuote(names(registry[[inference]])), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  list(
    inference = inference,
    type = type
  )
}


.selsurv_inference_for_type <- function(type, registry = .selsurv_boot_registry()) {
  hits <- names(Filter(
    function(x) type %in% names(x),
    registry
  ))

  if (length(hits) != 1L) {
    stop(
      "Internal bootstrap registry error for type ",
      shQuote(type),
      ".",
      call. = FALSE
    )
  }

  hits
}


# Input validation --------------------------------------------------------

.selsurv_assert_character_scalar <- function(x, argument, choices = NULL) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop(
      "`", argument, "` must be one non-missing character value.",
      call. = FALSE
    )
  }

  if (!is.null(choices) && !x %in% choices) {
    stop(
      "`", argument, "` must be one of ",
      paste(shQuote(choices), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(x)
}


.selsurv_validate_B <- function(B) {
  if (
    length(B) != 1L ||
      is.na(B) ||
      !is.numeric(B) ||
      !is.finite(B) ||
      B < 1 ||
      B != as.integer(B)
  ) {
    stop("`B` must be one positive integer.", call. = FALSE)
  }

  as.integer(B)
}


.selsurv_validate_conf_level <- function(conf_level) {
  if (
    length(conf_level) != 1L ||
      is.na(conf_level) ||
      !is.numeric(conf_level) ||
      !is.finite(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    stop(
      "`conf_level` must be one number strictly between 0 and 1.",
      call. = FALSE
    )
  }

  as.numeric(conf_level)
}


.selsurv_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (
    length(seed) != 1L ||
      is.na(seed) ||
      !is.numeric(seed) ||
      !is.finite(seed) ||
      seed < 0 ||
      seed != as.integer(seed)
  ) {
    stop("`seed` must be NULL or one non-negative integer.", call. = FALSE)
  }

  as.integer(seed)
}

.selsurv_validate_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(dots))
  }

  dot_names <- names(dots)
  if (is.null(dot_names) || any(!nzchar(dot_names))) {
    stop(
      "All arguments passed through `...` must be named.",
      call. = FALSE
    )
  }

  reserved <- c(
    "object",
    "inference",
    "type",
    "B",
    "conf_level",
    "seed"
  )
  conflicts <- intersect(dot_names, reserved)

  if (length(conflicts) > 0L) {
    stop(
      "Arguments reserved by `coxlasso_boot()` must not be supplied ",
      "through `...`: ",
      paste(shQuote(conflicts), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(dots)
}

