# selSurv

`selSurv` provides methods for inference after Cox-Lasso variable selection. The methods answer different questions; choosing an interval type alone does not determine its target.

## Implemented methods

Choose a method according to whether the original selected support is held fixed or selection is repeated.

| Function | Role | Inferential target |
|---|---|---|
| `coxlasso_fit()` | Cox-Lasso fitting and selection | common starting point |
| `coxlasso_debias()` | ridge-stabilized one-step correction | approximate coefficient-level target |
| `coxlasso_boot()` | score resampling or fixed-support model-based resampling | approximate debiased or fixed-support target |
| `coxlasso_boot_reselect()` | model-based resampling with repeated selection and refitting | selection-and-refit procedure |
| `coxlasso_split()` | sample splitting for selection and inference | split-sample target |

Lower-level method-specific functions remain available for compatibility and
development, but these functions define the recommended user-facing interface.

## Target distinction

`selSurv` separates three targets.

1. **Approximate coefficient-level target:** a one-step correction is applied to the Cox-Lasso estimate. This does not, by itself, establish valid high-dimensional confidence intervals.
2. **Fixed-support target:** the model selected in the observed data is held fixed during bootstrap refitting. This is not exact conditioning on the selection event.
3. **Procedure-level target:** selection and refitting are repeated in bootstrap samples. Unselected coordinates are recorded as zero, while failed fits remain missing.

The procedure-level target is exposed separately through `coxlasso_boot_reselect()`, because this procedure repeats Cox-Lasso selection and refitting rather than keeping the observed selected support fixed.

For `coxlasso_boot(..., inference = "conditional")`, the term `conditional` means that the support selected in the observed data is held fixed during bootstrap analysis. It does **not** mean formal conditioning on the complete selection event.

## Minimal example

```r
library(selSurv)

# Reproducible continuous-time example (without tied event times).
set.seed(2026)
n <- 120L
X <- matrix(rnorm(n * 4L), nrow = n)
colnames(X) <- paste0("x", seq_len(ncol(X)))
event_time <- rexp(n, rate = 0.12 * exp(X[, 1] - 0.7 * X[, 2]))
censor_time <- rexp(n, rate = 0.06)
time <- pmin(event_time, censor_time)
status <- as.integer(event_time <= censor_time)

fit <- coxlasso_fit(X, time, status, lambda_choice = "min", nfolds = 3)

fixed <- coxlasso_boot(
  fit,
  inference = "conditional",
  B = 200,
  seed = 2026
)

reselected <- coxlasso_boot_reselect(
  X, time, status,
  B = 200,
  methods = c(
    "fixed_support",
    "reselect_fixed_lambda",
    "reselect_cv_min"
  ),
  seed = 2026
)

db <- coxlasso_debias(fit)

wb <- coxlasso_boot(
  fit,
  inference = "debiased",
  type = "wild",
  B = 500,
  seed = 2026
)
sp <- coxlasso_split(X, time, status, seed = 1)
```

The bootstrap counts above illustrate the interface; they do not establish
adequate Monte Carlo precision. Assess simulation error and convergence for
the intended analysis.

## Research-software status

The bootstrap modules are research implementations. Computable intervals should not be read as general coverage guarantees. Report the target, selection rule, handling of tied event times, and relevant diagnostics, such as re-selection frequencies, support overlap, interval width, and fit success rates.


## Re-selection diagnostics

Procedure-level re-selection returns coefficient summaries together with
diagnostics describing how the selected model changes across bootstrap samples.

```r
rs_summary <- summary(reselected)

rs_summary$diagnostics
rs_summary$intervals
```

The diagnostics summarize support reproduction and overlap, selected-model
size, selection and refit success, and tuning-parameter variability where
applicable.

For originally selected coordinates, non-selection in a bootstrap sample is
represented by a zero coefficient. This is distinct from an unavailable
coefficient caused by a failed selection or refit. Re-selection frequency,
zero frequency, and numerical failure can therefore be interpreted separately.

`rng_mode = "legacy_nki70"` is retained only for reproducing the archived NKI70 production analysis. It is not the recommended default for new analyses.

## Vignettes

The package documentation gives each public inferential route a dedicated entry point:

- `getting-started.Rmd`: common Cox-Lasso fit and target selection;
- `debiased-bootstrap.Rmd`: debiased Wald inference, wild bootstrap, and Efron's bootstrap;
- `fixed-support-bootstrap.Rmd`: model-based bootstrap with the observed support held fixed;
- `sample-splitting.Rmd`: selection and inference on disjoint samples;
- `reselection.Rmd`: procedure-level repeated selection and refitting.

The vignettes are organised by inferential target rather than by a single preferred method.

<!-- BEGIN selSurv workflow guide -->
## Which inference method should I use?

`selSurv` supports several Cox-Lasso post-selection methods with **different inferential targets**. Start with `coxlasso_fit()` and choose the next step based on the scientific question.

| Goal | Recommended function | Interpretation |
| --- | --- | --- |
| Approximate full-model-coordinate inference | `coxlasso_debias(fit)` and optionally `coxlasso_boot(fit, inference = "debiased", type = c("wild", "efron"))` | Ridge-stabilized one-step correction and score-level bootstrap; exploratory research implementation, not a general high-dimensional validity guarantee |
| Fixed-support model-based bootstrap | `coxlasso_boot(fit, inference = "conditional")` | Resample and refit on the **original selected support**; does not condition formally on the complete selection event |
| Procedure-level re-selection | `coxlasso_boot_reselect(x, time, status, methods = c("fixed_support", "reselect_fixed_lambda", "reselect_cv_min"))` | Repeat the selection-and-refit procedure and inspect zero-padded coefficients, selection frequencies, support overlap and refit success |

**Additional interface:** `coxlasso_boot_parametric()` remains available for its own fixed-support and retuned re-selection procedures. It is not interchangeable with the three recommended entry points and has its own argument and result conventions.

**Argument naming and penalty scope:** `alpha` is an elastic-net mixing parameter in `coxlasso_fit()` and `coxlasso_boot_reselect()` (Lasso: `alpha = 1`), but a significance level in `coxlasso_boot_parametric()` and the lower-level score-bootstrap functions. `coxlasso_split()` instead uses `penalty_alpha` and `conf_level`. The documented inferential procedures have been examined primarily with Lasso screening (`alpha = 1`); this is not a general validity guarantee. Other penalties remain available for fitting. Debiased inference requires `alpha = 1`; split-sample, model-based bootstrap, and re-selection procedures warn when other mixing parameters are used because those settings have not been assessed separately.

For the complete decision guide in the R help system, run `?selSurv_workflows`. Cox-Lasso selection and unpenalized Cox refitting may use different conventions for tied event times. The *Efron event-score bootstrap* resamples event-level score contributions; it is unrelated to the *Efron approximation* used in Cox partial likelihood. The debiasing score calculation still requires particular caution when event times are tied.
<!-- END selSurv workflow guide -->
