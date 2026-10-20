# selSurv 0.1.1

Initial CRAN release.

- Updated the simulation-study reference to the published journal DOI.
- Clarified public method descriptions and documentation.
- Simplified the README while retaining development smoke tests.

# selSurv 0.1.0

Initial public GitHub release.

- Provides Cox-Lasso screening, support-restricted Cox refitting,
  split-sample inference, approximately debiased coefficient estimates,
  and score-level and model-based bootstrap workflows.
- Supports fixed-support inference and procedure-level re-selection,
  with explicit distinctions between the corresponding inferential targets.
- Corrects the Cox score and observed information calculations for
  tied event times under the Breslow and Efron conventions.
- Includes numerical regression tests against an independent
  Cox partial-likelihood reference implementation.
- Documents the supported Lasso penalty scope and the conventions
  used for fitting, refitting, and bootstrap inference.

# selSurv 0.1.0.9016

- Clarified the supported penalty scope of the inferential workflows. `coxlasso_fit()` continues to expose the `glmnet` elastic-net mixing parameter for fitting flexibility, while the inferential validation in selSurv currently concerns Lasso screening (`alpha = 1`).
- Restricted `coxlasso_debias()` and the associated wild and Efron score-bootstrap backends to Cox-Lasso pilots with `alpha = 1`; non-Lasso pilots now fail explicitly rather than returning apparently validated debiased inference.
- Kept split-sample, fixed-support model-based bootstrap, and procedure-level re-selection workflows computationally available for non-Lasso mixing parameters, but added explicit warnings that these settings have not been separately validated inferentially in selSurv.
- Added a package workflow guide documenting inferential targets, bootstrap interfaces, penalty-parameter naming, penalty scope, and tie-handling distinctions.
- Added focused regression tests for Lasso-only debiasing safeguards and non-Lasso warning behavior.

# selSurv 0.1.0.9015

- Made Cox tie conventions explicit across the model-based bootstrap workflows: penalized Cox-Lasso selection uses the Breslow convention, ordinary post-selection Cox refits use Efron's method, and model-based event generation uses an explicit Breslow cumulative baseline hazard.
- Added compatibility handling for installed `glmnet` versions that do not formally support `cox.ties`, avoiding silently ignored tie arguments while preserving the historical Breslow selection convention.
- Made the Cox-Lasso cross-validation convention explicit as deviance loss with grouped cross-validation, while retaining CV-min, 10 folds, and Lasso mixing `alpha = 1` as the general fitting defaults.
- Updated the general re-selection defaults to `conf_level = 0.95`, `seed = NULL`, and `standardize = TRUE`. A randomly generated root seed is stored when `seed = NULL`, after which bootstrap-generation and CV streams remain deterministic and method-order invariant.
- Retained `rng_mode = "legacy_nki70"` and its historical Breslow, random-stream, solver-tolerance, and standardization conventions solely for exact reproduction of the archived NKI70 analysis.
- Revised `coxlasso_split()` to separate interval confidence through `conf_level = 0.95` from Cox-Lasso mixing through `penalty_alpha = 1`. Historical `alpha` calls remain supported as a deprecated confidence-level bridge.
- Kept the historical direct bootstrap backend defaults unchanged for compatibility; the unified `coxlasso_boot()` interface remains the recommended bootstrap entry point with `B = 500`, `conf_level = 0.95`, and `seed = NULL`.
- Added regression tests for tie handling, Breslow baseline generation, cross-validation conventions, re-selection seed semantics, general-purpose defaults, and split-sample compatibility.

# selSurv 0.1.0.9014

- Standardized categorical data-frame inputs to stable treatment coding: factors, ordered factors, character variables, and logical variables now contribute \(k-1\) design-matrix coordinates using the first level as the reference, independent of data-frame column order.
- Ordered factors are treated as ordinary categorical predictors rather than using polynomial contrasts; logical predictors use `FALSE` as the reference level. Numeric and integer covariates remain numeric.
- Character predictors are converted to unordered factors. Users who require a particular reference category should supply a factor with an explicit level order or use `relevel()` before fitting.
- Selection and re-selection continue to operate on individual design-matrix coordinates, not on whole categorical variables as groups. Consequently, Lasso selection can depend on the chosen reference level.
- Added explicit checks for constant categorical predictors, unsupported data-frame column classes, empty inputs, missing or non-finite values, invalid or duplicated design-column names, and generated dummy-name collisions.
- Added design-coordinate checks so that post-selection refits and externally supplied original workflows cannot silently use reordered or differently coded covariates.
- Expanded regression tests for factor, ordered-factor, character, logical, numeric, sparse-matrix, reference-level, contrast-option, input-boundary, and procedure-level re-selection behavior.

# selSurv 0.1.0.9013

- Centralized covariate design-matrix normalization across the main Cox-Lasso workflow and procedure-level re-selection.
- `coxlasso_boot_reselect()` now accepts data-frame covariates using the same `model.matrix(~ . - 1)` expansion already used by `coxlasso_fit()`, including factor and logical columns.
- Added explicit rejection of missing and non-finite covariate values before model fitting.
- Added regression tests for numeric, factor/logical data-frame, and sparse-matrix inputs.
- Preserved the existing coding convention, numerical behavior, and archived NKI70 re-selection parity.

# selSurv 0.1.0.9012

- Added an experimental procedure-level model-based bootstrap with Cox-Lasso re-selection through `coxlasso_boot_reselect()`.
- Added fixed-support, fixed-penalty re-selection, bootstrap CV-min, CV-1SE, and AIC tuning variants.
- Non-selected coordinates are recorded as zero; failed non-empty unpenalized Cox refits remain missing rather than being converted to zero.
- Added method-level diagnostics for support overlap, model size, selection/refit success, tuning parameters, and coordinate-level re-selection frequencies.
- Added method-order-invariant random streams for normal use.
- Added `rng_mode = "legacy_nki70"` solely for reproducibility of the archived NKI70 analysis. The B=20 regression gate reproduces the archived support diagnostics, including the historical AIC refit failure pattern.
- The existing `coxlasso_boot(inference = "conditional")` semantics remain fixed-support and are not broadened to include procedure-level re-selection.

# selSurv 0.1.0.9011

- Integrated the August 2026 unified `coxlasso_boot()` API with the recovered
  v12 package source.
- Added dispatch between debiased wild bootstrap and Efron's bootstrap and fixed-support
  model-based bootstrap.
- Added compatibility with the established v12 `selSurv_debias` class.
- Kept procedure-level re-selection outside the first unified bootstrap API;
  re-selection remains available through `coxlasso_boot_parametric()`.
- Added validation/dispatch tests and development use-case documentation.
- Archived the historical selSurv IWSM text as an unpublished draft and noted
  that the retained IWSM 2026 proceedings contribution was the PoSI work.
- Removed local `.RData`/`.Rhistory` workspace artifacts from the source tree.

# selSurv 0.1.0.9000

Initial research implementation.

- Added Cox-Lasso fitting workflow.
- Added support-restricted Cox refitting.
- Added approximate debiasing workflow.
- Added multiplier wild and event-level score bootstraps.
- Added model-based bootstrap with fixed-support, fixed-lambda re-selection, and retuned re-selection modes.
- Added split-sample workflow.
- Added target and method summary helpers.
- Added vignette explaining target definitions and workflow interpretation.

# selSurv 0.1.0.9001

- Added user-facing `coef` and `se` fields to `coxlasso_fit()` and `coxlasso_refit()` outputs.
- Added support-specific `coef_support` and `se_support` fields to `coxlasso_refit()`.


## selSurv 0.1.0.9002

- Added documentation for S3 print methods.
- Rewrote workflow Rd usage/argument documentation.
- Moved local smoke-test examples and implementation notes under `inst/`.
- Reduced R CMD check warnings/notes caused by package documentation layout.


## selSurv 0.1.0.9003

- Fixed Rd code/documentation mismatches for `coxlasso_wildboot()` and
  `coxlasso_efron_events()`.
- Added `.Rbuildignore` entries for stale top-level example/development files.
- Shortened the `coxlasso_fit()` Rd usage block.


## selSurv 0.1.0.9004

- Added explicit roxygen usage block for `coxlasso_fit()` to remove Rd line-width NOTE.
- Broke the `coxlasso_fit()` function signature across lines.
- Removed empty `inst/extdata` from the source tree.


## selSurv 0.1.0.9005

- Removed the obsolete combined `coxlasso_workflows.Rd` file to avoid duplicate Rd aliases after roxygen generates individual workflow documentation.
- Removed the manual `coxlasso_fit.Rd` file so roxygen can regenerate it from the explicit multi-line `@usage` block.


## selSurv 0.1.0.9006

- Added `plot_selection_diagnostics()` for base-R visualization of re-selection and zero-frequency diagnostics.
- Added tests for diagnostic comparison and plotting helpers.


## selSurv 0.1.0.9007

- Improved `plot_selection_diagnostics()` layout; default legend is now a compact top legend and neutral greys are used by default.


## selSurv 0.1.0.9008

- Added a prettier default Cleveland-style dot plot to `plot_selection_diagnostics()`.
- Kept grouped bars available via `plot_type = "bar"`.
- Improved plot readability for many variables and small RStudio plot panes.


## selSurv 0.1.0.9009

- Added a colorblind-friendly color palette to `plot_selection_diagnostics()`.
- Added a robust `legend_position = "side"` default using a separate legend panel.
- Kept grouped bars available and improved the color bar plot option.


## selSurv 0.1.0.9010

- Replaced the default side legend by a top legend panel in `plot_selection_diagnostics()`.
- This avoids legend clipping in the RStudio zoom window and keeps the main plot area readable.
