#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}

B <- as.integer(get_arg("B", "20"))
outdir <- get_arg("outdir", "nki70_parity")
core <- get_arg("core", "R/reselection-core.R")
api <- get_arg("api", "R/reselection-api.R")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source(core)
if (file.exists(api)) source(api)

for (pkg in c("survival", "glmnet", "penalized")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg)
  }
}

data("nki70", package = "penalized", envir = environment())
stopifnot(is.data.frame(nki70), nrow(nki70) == 144L)

X <- scale(as.matrix(nki70[, 8:77]))
storage.mode(X) <- "double"
time <- as.numeric(nki70$time)
status <- as.integer(nki70$event)

expected_support <- c("QSCN6L1", "ZNF533", "IGFBP5.1", "PRC1")
archived_lambda <- 0.0983

original <- fit_original_fixed_lambda(
  X, time, status,
  lambda = archived_lambda,
  standardize = FALSE,
  ties = "breslow"
)
observed_support <- colnames(X)[original$selection$support]

cat("Observed support:\n")
print(observed_support)
if (!identical(observed_support, expected_support)) {
  stop(
    "Observed support does not match archived NKI70 support. Got: ",
    paste(observed_support, collapse = ", ")
  )
}

expected_refit <- c(
  QSCN6L1 = 0.2000196,
  ZNF533 = -0.2884900,
  IGFBP5.1 = 0.3190741,
  PRC1 = 0.5780386
)
refit_now <- original$refit$beta[expected_support]
refit_diff <- refit_now - expected_refit
cat("\nOriginal refit audit:\n")
print(data.frame(variable = expected_support, current = refit_now,
                 expected = expected_refit, difference = refit_diff))
if (max(abs(refit_diff)) > 5e-4) {
  stop("Original NKI70 refit is not numerically aligned with the archive.")
}

methods <- c(
  "fixed_support",
  "reselect_fixed_lambda",
  "reselect_cv_min",
  "reselect_aic"
)

fit <- run_reselection_bootstrap(
  X, time, status,
  B = B,
  methods = methods,
  conf_level = 0.90,
  seed = 20260827,
  nfolds = 10,
  original = original,
  standardize = FALSE,
  ties = "breslow",
  rng_mode = "legacy_nki70",
  progress = TRUE
)

diag <- method_level_diagnostics(fit)
ci <- summarise_reselection(fit, conf_level = 0.90)

cat("\nGenerator audit:\n")
cat("mode =", fit$generator$generator_mode %||% "<unknown>", "\n")
cat("administrative tau =", fit$generator$administrative_tau, "\n")
if (!is.null(fit$generator$censor_values)) {
  cat("reverse-KM support points =", length(fit$generator$censor_values), "\n")
}

write.csv(diag, file.path(outdir, "method_diagnostics.csv"), row.names = FALSE)
write.csv(ci, file.path(outdir, "coefficient_summary.csv"), row.names = FALSE)
saveRDS(fit, file.path(outdir, sprintf("nki70_reselection_B%d.rds", B)))

cat("\nMethod diagnostics:\n")
print(diag)
cat("\nBootstrap-SE Wald intervals:\n")
print(ci[, c("method", "variable", "estimate",
             "bootstrap_se_lower", "bootstrap_se_upper",
             "reselection_frequency", "boot_zero_frequency")])

if (B == 20L) {
  expected_diag20 <- data.frame(
    method = methods,
    exact_support_reproduction = c(1.000, 0.000, 0.000, 0.000),
    mean_jaccard_original = c(1.0000000, 0.3809463, 0.2768506, 0.1825172),
    mean_n_selected_boot = c(4.00, 5.80, 8.05, 42.95),
    sd_n_selected_boot = c(0.000000, 1.735087, 3.734265, 28.757562),
    refit_success_rate = c(1.00, 1.00, 1.00, 0.95),
    stringsAsFactors = FALSE
  )

  got20 <- merge(
    expected_diag20[, "method", drop = FALSE],
    diag,
    by = "method",
    sort = FALSE
  )
  cols20 <- setdiff(names(expected_diag20), "method")
  diff20 <- sapply(cols20, function(v) {
    max(abs(got20[[v]] - expected_diag20[[v]]), na.rm = TRUE)
  })

  parity20 <- list(
    diagnostic_max_abs_difference = diff20,
    diagnostic_tolerance = 1e-5,
    passed = all(diff20 <= 1e-5)
  )
  dput(parity20, file = file.path(outdir, "PARITY_B20_RESULT.R"))
  cat("\nArchived B=20 diagnostic parity:\n")
  print(parity20)
  if (!isTRUE(parity20$passed)) quit(status = 2L)
}


if (B == 500L) {
  expected_diag <- data.frame(
    method = methods,
    exact_support_reproduction = c(1.000, 0.006, 0.000, 0.000),
    mean_jaccard_original = c(1.0000000, 0.3354380, 0.2563048, 0.1043163),
    mean_n_selected_boot = c(4.000, 5.880, 9.238, 51.522),
    sd_n_selected_boot = c(0.000000, 2.063520, 4.364798, 23.476452),
    refit_success_rate = c(1.00, 1.00, 1.00, 0.93),
    stringsAsFactors = FALSE
  )

  expected_ci <- data.frame(
    method = rep(c("fixed_support", "reselect_fixed_lambda", "reselect_cv_min"), each = 4),
    variable = rep(expected_support, times = 3),
    lower = c(
      -0.071, -0.584,  0.080, 0.257,
      -0.034, -0.643, -0.086, 0.129,
      -0.068, -0.684, -0.123, 0.079
    ),
    upper = c(
       0.471, 0.007, 0.558, 0.899,
       0.434, 0.066, 0.724, 1.027,
       0.468, 0.107, 0.761, 1.077
    ),
    stringsAsFactors = FALSE
  )

  got_diag <- merge(
    expected_diag[, "method", drop = FALSE], diag,
    by = "method", sort = FALSE
  )
  diag_cols <- setdiff(names(expected_diag), "method")
  diag_diff <- sapply(diag_cols, function(v) {
    max(abs(got_diag[[v]] - expected_diag[[v]]), na.rm = TRUE)
  })

  got_ci <- ci[ci$method %in% unique(expected_ci$method) &
                 ci$variable %in% expected_support,
               c("method", "variable", "bootstrap_se_lower", "bootstrap_se_upper")]
  got_ci <- merge(expected_ci[, c("method", "variable")], got_ci,
                  by = c("method", "variable"), sort = FALSE)
  ci_diff <- c(
    lower = max(abs(got_ci$bootstrap_se_lower - expected_ci$lower), na.rm = TRUE),
    upper = max(abs(got_ci$bootstrap_se_upper - expected_ci$upper), na.rm = TRUE)
  )

  parity <- list(
    diagnostic_max_abs_difference = diag_diff,
    interval_max_abs_difference = ci_diff,
    diagnostic_tolerance = 0.03,
    interval_tolerance = 0.03,
    passed = all(diag_diff <= 0.03) && all(ci_diff <= 0.03)
  )
  dput(parity, file = file.path(outdir, "PARITY_RESULT.R"))
  cat("\nArchived B=500 parity differences:\n")
  print(parity)
  if (!isTRUE(parity$passed)) quit(status = 2L)
}

cat("\nNKI70 parity run completed.\n")
