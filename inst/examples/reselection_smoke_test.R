## Development smoke test for the procedure-level re-selection API.
## Intended for local/CI execution, not as an expensive help-page example.

set.seed(20261007)
n <- 100L
p <- 6L
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", seq_len(p))
beta <- c(0.8, -0.6, 0.4, 0, 0, 0)
T <- rexp(n, rate = 0.06 * exp(drop(X %*% beta)))
C <- rexp(n, rate = 0.04)
time <- pmin(T, C)
status <- as.integer(T <= C)

fit <- coxlasso_boot_reselect(
  X, time, status,
  B = 5L,
  methods = c("fixed_support", "reselect_fixed_lambda", "reselect_cv_min"),
  conf_level = 0.90,
  seed = 17L,
  nfolds = 5L,
  standardize = FALSE,
  ties = "breslow",
  progress = FALSE
)

s <- summary(fit)
print(s)

stopifnot(inherits(fit, "selsurv_reselection"))
stopifnot(inherits(s, "summary.selsurv_reselection"))
stopifnot(is.data.frame(s$diagnostics))
stopifnot(nrow(s$diagnostics) == 3L)
