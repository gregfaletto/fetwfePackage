# Run twfeCovs on Simulated Data

This function runs the bridge-penalized extended two-way fixed effects
estimator
([`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md))
on simulated data. It is simply a wrapper for
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md):
it accepts an object of class `"FETWFE_simulated"` (produced by
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md))
and unpacks the necessary components to pass to
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md).
So the outputs match
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md),
and the needed inputs match their counterparts in
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md).

## Usage

``` r
twfeCovsWithSimulatedData(
  simulated_obj,
  verbose = FALSE,
  alpha = 0.05,
  add_ridge = FALSE,
  allow_no_never_treated = TRUE,
  se_type = "default",
  ci_type = c("simultaneous", "pointwise")
)
```

## Arguments

- simulated_obj:

  An object of class `"FETWFE_simulated"` containing the simulated panel
  data and design matrix.

- verbose:

  Logical; if TRUE, more details on the progress of the function will be
  printed as the function executes. Default is FALSE.

- alpha:

  Numeric; function will calculate (1 - `alpha`) confidence intervals
  for the cohort average treatment effects that will be returned in
  `catt_df`.

- add_ridge:

  (Optional.) Logical; if TRUE, adds a small amount of ridge
  regularization to the (untransformed) coefficients to stabilize
  estimation. Default is FALSE.

- allow_no_never_treated:

  (Optional.) Logical; if `TRUE` (default) and the input panel contains
  no never-treated units, the panel is auto-truncated by dropping time
  periods at and after the latest cohort's start time — the units in
  that latest cohort then serve as the never-treated comparison group in
  the retained sub-panel — with a warning naming the dropped periods. If
  `FALSE`, the estimator stops with an error in this case (the package's
  behavior prior to version 1.5.6). The argument has no effect when the
  input already contains never-treated units. Default is `TRUE`.

- se_type:

  Character; one of `"default"`, `"conservative"`, or `"cluster"`.
  `"default"` returns the tight Gaussian variance
  `sqrt(att_var_1 + att_var_2)` from Theorem (c\$'\$) under Assumption
  (Psi-IF) (asymptotically exact for the package's default cohort
  sample-proportions estimator); `"conservative"` returns the
  Cauchy-Schwarz upper bound from Theorem (c) (use only when the
  propensity-score estimator violates (Psi-IF)); `"cluster"` is an
  *experimental* unit-clustered Liang-Zeger sandwich SE on the
  OLS-selected support (see the companion vignette `inference_vignette`
  for the formula, the assumptions, and the theory-pending caveat).
  Default is `"default"`. v1.12.0 introduced the tight Gaussian default;
  versions \<= 1.11.7 used the conservative Cauchy-Schwarz formula as
  the default.

- ci_type:

  Character; one of `"simultaneous"` (default) or `"pointwise"`.
  Controls the confidence-interval bounds reported for the
  cohort-specific ATTs (in `catt_df`). `"simultaneous"` reports
  parametric simultaneous (family-wise, uniform) bands computed via
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md):
  the band covers all cohort effects jointly with probability
  `1 - alpha`, matching the default presentation of
  `did::aggte(cband = TRUE)`. `"pointwise"` reports per-effect Wald
  intervals (each covers its own effect with probability `1 - alpha`, no
  joint guarantee — the behavior of versions \<= 1.15.1). Both the
  interval bounds and the per-cohort p-values (`p_value`) follow
  `ci_type` (single-step max-T multiplicity-adjusted under
  `"simultaneous"`, per-cohort Wald under `"pointwise"`; \#200); the
  standard errors (`se`) are identical under both settings. `twfeCovs`
  estimates a single pooled effect per cohort, so only the cohort family
  is affected (it has no event-study surface). When standard errors are
  unavailable (e.g., a rank-deficient design) the bounds are `NA` under
  both settings. Default is `"simultaneous"`.

## Value

An object of class `twfeCovs` containing the following elements:

- att_hat:

  The estimated overall average treatment effect for a randomly selected
  treated unit.

- att_se:

  A standard error for the ATT. If `indep_counts` was provided, this
  standard error is asymptotically exact; otherwise, it is
  asymptotically conservative. If the Gram matrix is not invertible,
  this will be NA.

- att_p_value:

  A two-sided p-value for the overall ATT against the null
  `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
  `att_se` is zero or `NA`. Standard post-OLS interpretation; `twfeCovs`
  does not perform selection.

- catt_hats:

  A named vector containing the estimated average treatment effects for
  each cohort.

- catt_ses:

  A named vector containing the (asymptotically exact, non-conservative)
  standard errors for the estimated average treatment effects within
  each cohort. If the Gram matrix is not invertible, the entries are NA.

- cohort_probs:

  A vector of the estimated probabilities of being in each cohort
  conditional on being treated, which was used in calculating `att_hat`.
  If `indep_counts` was provided, `cohort_probs` was calculated from
  that; otherwise, it was calculated from the counts of units in each
  treated cohort in `pdata`.

- catt_df:

  A data frame (with S3 class `c("catt_df", "data.frame")`) displaying
  the cohort names (`cohort`), average treatment effects (`estimate`),
  standard errors (`se`), `1 - alpha` confidence interval bounds
  (`ci_low`, `ci_high`), and per-cohort p-values (`p_value`). No
  `selected` column; `twfeCovs` does not perform selection. The
  `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0
  Title-Case column names (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`,
  `ConfIntHigh`, `P_value`) [`stop()`](https://rdrr.io/r/base/stop.html)
  with a migration message pointing to the new name. See `NEWS.md` for
  the rename table.

- beta_hat:

  The full vector of estimated coefficients.

- treat_inds:

  The indices of `beta_hat` corresponding to the treatment effects for
  each cohort.

- treat_int_inds:

  The indices of `beta_hat` corresponding to the interactions between
  the treatment effects for each cohort and the covariates.

- sig_eps_sq:

  Either the provided `sig_eps_sq` or the estimated one, if a value
  wasn't provided.

- sig_eps_c_sq:

  Either the provided `sig_eps_c_sq` or the estimated one, if a value
  wasn't provided.

- X_ints:

  The design matrix created containing all interactions, time and cohort
  dummies, etc.

- y:

  The vector of responses, containing `nrow(X_ints)` entries.

- X_final:

  The design matrix after applying the change in coordinates to fit the
  model and also multiplying on the left by the square root inverse of
  the estimated covariance matrix for each unit.

- y_final:

  The final response after multiplying on the left by the square root
  inverse of the estimated covariance matrix for each unit.

- N:

  The final number of units that were in the data set used for
  estimation (after any units may have been removed because they were
  treated in the first time period).

- T:

  The number of time periods in the final data set.

- G:

  The final number of treated cohorts that appear in the final data set.

- R:

  Deprecated alias for `G`, retained for backward compatibility;
  populated with the same value. Use `G`. Will be removed in a future
  release.

- d:

  The final number of covariates that appear in the final data set
  (after any covariates may have been removed because they contained
  missing values or all contained the same value for every unit).

- p:

  The final number of columns in the full set of covariates used to
  estimate the model.

- calc_ses:

  Logical indicating whether standard errors were calculated.

- cohort_probs_overall:

  A vector of the estimated cohort probabilities on the overall sample
  (treated and untreated), used in computing the variance of the overall
  ATT.

- indep_counts_used:

  Logical scalar; `TRUE` if a valid `indep_counts` argument was provided
  and used for asymptotically-exact ATT inference, `FALSE` otherwise.

- se_type:

  Character scalar; the `se_type` argument the user passed (`"default"`,
  `"conservative"`, or `"cluster"`).

- alpha:

  The alpha level used for confidence intervals.

- ci_type:

  Character scalar; the `ci_type` argument the user passed
  (`"simultaneous"` or `"pointwise"`), controlling whether the reported
  `catt_df` confidence-interval bounds are simultaneous (family-wise) or
  pointwise.

- y_mean:

  Numeric scalar; mean of the original (pre-centering) response. Stored
  so downstream methods (`augment()`,
  [`predict()`](https://rdrr.io/r/stats/predict.html)) can return fitted
  values on the original-response scale.

- response_col_name:

  Character scalar; the response column name in the original `pdata`.

- time_var, unit_var, treatment:

  Character scalars; the corresponding arguments the user passed.

- covs:

  Character vector; the original `covs` argument (pre-factor-
  expansion).

- internal:

  A list containing internal outputs that are typically not needed for
  interpretation, packaged here for parity with
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  so downstream consumers can use a single canonical access path across
  all four estimator classes (#144). The first five sub-slots (`X_ints`,
  `y`, `X_final`, `y_final`, `calc_ses`) are also duplicated at top
  level for backward compat; `variance_components` and `first_year` live
  only under `$internal`:

  X_ints

  :   The design matrix containing all interactions, time and cohort
      dummies, etc. Same value as top-level `X_ints`.

  y

  :   The vector of responses. Same as top-level `y`.

  X_final

  :   The design matrix after the change-of-coordinates step. Same as
      top-level `X_final`.

  y_final

  :   The transformed response vector. Same as top-level `y_final`.

  calc_ses

  :   Logical indicating whether standard errors were calculated. Same
      as top-level `calc_ses`.

  variance_components

  :   A list exposing the two variance pieces (`att_var_1`, `att_var_2`)
      plus paper-notation counterparts (`V_1`, `V_2`) and unit-scaled
      variance estimators (`tilde_v_N`, `hat_v_N`, `tilde_v_N_C`,
      `tilde_v_N_C_pi_hat`, `tilde_v_N_C_pi_hat_cons`,
      `tilde_v_N_cons`). The Wald CI is
      `[hat_T_N +- qnorm(1-alpha/2) * sqrt(tilde_v_N / N)]` (paper Eq.
      `conf.int.form`). New in v1.12.0 (issue \#141 + \#146).

  first_year

  :   Integer or numeric scalar; the first (earliest) `time_var` value
      in the panel after `idCohorts()` processing. Consumed by
      [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
      to map `cohort_probs`' cohort labels (treatment-start years) to
      1-based panel-time-index offsets when the labels are
      integer-coercible. New in v1.13.3 (issue \#174).

The returned object is an S3-classed `"twfeCovs"` list with
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html), and
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
methods, matching the three sibling estimators.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) is
intentionally not defined —
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
estimates one pooled effect per cohort, so there is no per-(cohort,
time) / event-study structure to plot. `augment()` is intentionally not
defined — the coefficient vector lives in a reduced cohort-level basis
that `augment()`'s fitted-value path does not match. Both raise an
informative error (#58).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate coefficients
  coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)

  # Simulate data using the coefficients
  sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)

  result <- twfeCovsWithSimulatedData(sim_data)
} # }
```
