library(selSurv)
library(survival)

set.seed(1)
n <- 120
p <- 10

X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", seq_len(p))

beta <- c(0.8, -0.6, 0.5, rep(0, p - 3))
linpred <- drop(X %*% beta)

time <- rexp(n, rate = exp(linpred) * 0.05)
cens <- rexp(n, rate = 0.03)

y <- pmin(time, cens)
status <- as.integer(time <= cens)

fit <- coxlasso_fit(X = X, time = y, status = status, nfolds = 5)
ref <- coxlasso_refit(fit)

boot_fixed <- coxlasso_boot_parametric(
  fit = fit,
  B = 100,
  workflow = "fixed_support",
  nfolds = 5,
  seed = 123
)

boot_t2 <- coxlasso_boot_parametric(
  fit = fit,
  B = 100,
  workflow = "reselect_fixed_lambda",
  nfolds = 5,
  seed = 123
)

boot_retuned <- coxlasso_boot_parametric(
  fit = fit,
  B = 100,
  workflow = "reselect_retuned",
  nfolds = 5,
  seed = 123
)

compare_diag <- compare_selection_diagnostics(
  fixed = boot_fixed,
  reselect_fixed_lambda = boot_t2,
  reselect_retuned = boot_retuned,
  digits = 2
)

print(compare_diag)
print(boot_t2$ci)
print(boot_t2$success_rate)


# Optional diagnostic plots in an interactive session:
# plot_selection_diagnostics(compare_diag, metric = "reselect", sort_by = "t2_reselect")
# plot_selection_diagnostics(compare_diag, metric = "zero", sort_by = "t2_zero")


# Prettier diagnostic dot plots:
# plot_selection_diagnostics(compare_diag, metric = "reselect", plot_type = "dot", sort_by = "t2_reselect")
# plot_selection_diagnostics(compare_diag, metric = "zero", plot_type = "dot", sort_by = "t2_zero")


# Robust legend placement:
# plot_selection_diagnostics(compare_diag, metric = "reselect", plot_type = "bar", legend_position = "top_panel")
