# selSurv implementation notes

This package was created from the IWSM paper draft `Schemet_IWSM2026_clean_submission.tex`, which described the intended package interface but did not itself contain R package code.

## Implemented modules

### `coxlasso_fit()`
Fits a Cox-Lasso model with `glmnet::cv.glmnet()`, stores the selected support, and refits an unpenalized Cox model on the selected variables.

### `coxlasso_refit()`
Refits a Cox model on a supplied support and returns a zero-padded coefficient vector.

### `coxlasso_debias()`
Implements a ridge-stabilized one-step correction using a Cox score and information approximation. This is a research implementation intended to provide a common workflow interface; it is not a fully validated high-dimensional debiasing theorem implementation.

### `coxlasso_wildboot()` and `coxlasso_efron_events()`
Use score-contribution perturbations of the debiased estimator. These target the debiased estimator, not the selection-refit procedure-level target.

### `coxlasso_boot_parametric()`
Implements model-based bootstrap workflows:

- `workflow = "fixed_support"`: freezes the original selected support and refits the Cox model.
- `workflow = "reselect_fixed_lambda"`: reuses the original Cox-Lasso penalty and re-selects inside each bootstrap sample.
- `workflow = "reselect_retuned"`: retunes the penalty inside each bootstrap sample.

The re-selection workflows return zero-padded bootstrap coefficient vectors and diagnostics such as re-selection and zero frequencies.

### `coxlasso_split()`
Implements split-sample selection/refit inference.

## Target interpretation

The package distinguishes:

1. full-model/debiased targets,
2. fixed-support selected-model targets,
3. zero-padded procedure-level targets.

The procedure-level target is the relevant target for model-based bootstrap re-selection.

## Validation status

The current environment did not provide `Rscript`, so I could not run `R CMD check` here. The package structure, function files, namespace, README, vignette, and tests are included. Before using the package for final results, run:

```r
install.packages(c("survival", "glmnet", "MASS", "testthat", "knitr", "rmarkdown"))
devtools::document()
devtools::test()
devtools::check()
```

or from shell:

```bash
R CMD build selSurv
R CMD check selSurv_*.tar.gz
```

## Important methodological caveat

The bootstrap and debiasing functions are deliberately labelled as research implementations. Final paper claims should still be phrased in terms of the target and the assumptions: fixed-support, full-model/debiased, or procedure-level with support-score replicability.
