# Run FETWFE on Simulated Data

This function runs the fused extended two-way fixed effects estimator
([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md))
on simulated data. It is simply a wrapper for
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md):
it accepts an object of class `"FETWFE_simulated"` (produced by
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md))
and unpacks the necessary components to pass to
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).
So the outputs match
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
and the needed inputs match their counterparts in
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).

## Usage

``` r
fetwfeWithSimulatedData(
  simulated_obj,
  lambda.max = NA,
  lambda.min = NA,
  nlambda = 100,
  q = 0.5,
  verbose = FALSE,
  alpha = 0.05,
  add_ridge = FALSE,
  allow_no_never_treated = TRUE,
  se_type = "default",
  lambda_selection = "cv",
  cv_folds = 10L,
  cv_seed = NULL,
  ci_type = c("simultaneous", "pointwise"),
  fusion_structure = c("cohort", "event_study"),
  fusion_matrix = NULL,
  gls = TRUE
)
```

## Arguments

- simulated_obj:

  An object of class `"FETWFE_simulated"` containing the simulated panel
  data and design matrix.

- lambda.max:

  (Optional.) Numeric. A penalty parameter `lambda` will be selected
  over a grid search by BIC in order to select a single model. The
  largest `lambda` in the grid will be `lambda.max`. If no `lambda.max`
  is provided, one will be selected automatically. For `lambda <= 1`,
  the model will be sparse, and ideally all of the following are true at
  once: the smallest model (the one corresponding to `lambda.max`)
  selects close to 0 features, the largest model (the one corresponding
  to `lambda.min`) selects close to `p` features, `nlambda` is large
  enough so that models are considered at every feasible model size, and
  `nlambda` is small enough so that the computation doesn't become
  infeasible. You may want to manually tweak `lambda.max`, `lambda.min`,
  and `nlambda` to try to achieve these goals, particularly if the
  selected model size is very close to the model corresponding to
  `lambda.max` or `lambda.min`, which could indicate that the range of
  `lambda` values was too narrow. You can use the function outputs
  `lambda.max_model_size`, `lambda.min_model_size`, and
  `lambda_star_model_size` to try to assess this. Default is NA.

- lambda.min:

  (Optional.) Numeric. The smallest `lambda` penalty parameter that will
  be considered. See the description of `lambda.max` for details.
  Default is NA.

- nlambda:

  (Optional.) Integer. The total number of `lambda` penalty parameters
  that will be considered. See the description of `lambda.max` for
  details. Default is 100.

- q:

  (Optional.) Numeric; determines what `L_q` penalty is used for the
  fusion regularization. `q` = 1 is the lasso, and for 0 \< `q` \< 1, it
  is possible to get standard errors and confidence intervals. `q` = 2
  is ridge regression. See Faletto (2025) for details. Default is 0.5.

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
  (Psi-IF); this is asymptotically exact for the package's default
  cohort sample-proportions estimator and for every standard
  propensity-score estimator that satisfies (Psi-IF) (multinomial logit,
  any GLM on `W | X`, kernel/series regression of `1{W = g}` on `X`).
  `"conservative"` returns the Cauchy-Schwarz upper bound
  `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` from
  Theorem (c); use only if the propensity-score estimator violates
  (Psi-IF) (e.g., a Robins-Rotnitzky-augmented doubly-robust estimator,
  which the package does not currently implement). `"cluster"` is an
  *experimental* unit-clustered Liang-Zeger sandwich SE on the
  bridge-selected support (see the companion vignette
  `inference_vignette` for the formula, the assumptions, and the
  theory-pending caveat); only meaningful when `q < 1` (the bridge
  oracle property is required), and for `q >= 1` the SE will be `NA`
  regardless of `se_type`. The default value of `"default"` corresponds
  to the new tight Gaussian default introduced in version 1.12.0;
  previous versions used the conservative Cauchy-Schwarz formula as the
  default. To recover the prior conservative default behavior, pass
  `se_type = "conservative"`.

- lambda_selection:

  Character; method for selecting the bridge penalty parameter `lambda`.
  Either `"cv"` (10-fold cross-validation on `cv.grpreg`; the v1.13.0+
  default) or `"bic"` (BIC over the `grpreg` lambda grid; the prior
  default for v1.12.0 and earlier). The default changed in v1.13.0 to
  address a finite-sample bias issue documented in simulation studies
  (see issue \#164): under the prior BIC default, the overall-ATT
  estimator was biased toward zero at moderate sample sizes, producing
  95% confidence intervals whose empirical coverage was as low as 0.00
  in some regimes. Cross-validation restores near-nominal coverage in
  every regime tested. To recover the prior behavior — for example, when
  reproducing analyses run against v1.12.0 or earlier — pass
  `lambda_selection = "bic"`. See the inference vignette section
  "Choosing the bridge penalty parameter" for details.

- cv_folds:

  Integer; number of folds for the CV path. Ignored when
  `lambda_selection = "bic"`. Default is 10.

- cv_seed:

  Integer or `NULL`; the seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  the `cv.grpreg()` call, controlling fold assignment. If `NULL` (the
  default), the seed defaults internally to `as.integer(N * T)` so
  consecutive calls on the same dataset are reproducible without the
  user having to specify a seed. The seed actually used is stored on the
  returned object as `cv_seed`. Ignored when `lambda_selection = "bic"`.

- ci_type:

  Character; one of `"simultaneous"` (default) or `"pointwise"`.
  Controls the confidence-interval bounds reported for the
  cohort-specific ATTs (in `catt_df`) and the event-study effects (from
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md),
  shown by `print` / `summary` / `plot`, and surfaced by
  [`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) on
  the fitted object and on the
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
  /
  [`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
  outputs). `"simultaneous"` reports parametric simultaneous
  (family-wise, uniform) bands computed via
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md):
  each family's band covers all of its effects jointly with probability
  `1 - alpha`, matching the default presentation of
  `did::aggte(cband = TRUE)`. `"pointwise"` reports per-effect Wald
  intervals (each covers its own effect with probability `1 - alpha`, no
  joint guarantee — the behavior of versions \<= 1.15.1). Both the
  interval bounds and the per-cohort p-values (`p_value`) follow
  `ci_type`: under `"simultaneous"` the `p_value` is the single-step
  max-T multiplicity- adjusted p-value matching the band, under
  `"pointwise"` the per-cohort Wald p-value (#200). The standard errors
  (`se`) and selection flags (`selected`) are identical under both
  settings, and the overall-ATT confidence interval (a single scalar) is
  unaffected. When standard errors are unavailable (`q >= 1`, or a
  rank-deficient design) the bounds are `NA` under both settings.
  Default is `"simultaneous"`.

- fusion_structure:

  Character; one of `"cohort"` or `"event_study"`. `"cohort"` (the
  default) uses the within-cohort / between-cohort two-way fusion
  penalty. `"event_study"` instead fuses treatment effects at the same
  time since treatment (event time `e = t - g`) across cohorts. The
  event-study penalty carries the same theoretical guarantees as the
  default (Faletto 2025); only the treatment-effect fusion structure
  changes.

- fusion_matrix:

  (Optional.) Numeric matrix or `NULL` (the default). An advanced-use
  override: a user-supplied `num_treats x num_treats` forward
  differences matrix `D_N` for the treatment-effect block, encoding an
  arbitrary fusion structure beyond the two built-ins. When non-`NULL`
  it overrides `fusion_structure` for the treatment-effect block only
  (the fixed-effect blocks are unchanged); the estimator uses
  `solve(fusion_matrix)` internally. The rows/columns are interpreted in
  the cohort-major `(g, t)` order used internally for the treatment
  effects (the order `getFirstInds()` / `getTreatInds()` encode):
  row/column `i` corresponds to base treatment effect `i`, with cohort
  `g` occupying rows `first_inds[g]:(first_inds[g + 1] - 1)` ordered by
  event time. `num_treats` equals `T * G - G * (G + 1) / 2`.
  `fusion_matrix` must be a finite, invertible numeric matrix of that
  exact dimension (otherwise
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  stops). Under the paper's fixed-dimension scoping, *any* finite
  invertible `D_N` of that dimension inherits the paper's inferential
  guarantees: the theory depends on `D_N` only through its invertibility
  and singular-value bounds (Assumption (D) of Faletto 2025), which a
  fixed invertible matrix automatically satisfies, and swapping in a
  different `D_N` from this class changes only constant factors. A
  numerically near-singular (ill-conditioned) `D_N` still yields a valid
  point estimator but emits a
  [`warning()`](https://rdrr.io/r/base/warning.html) that its inverse
  may be unreliable. Default is `NULL` (use the built-in
  `fusion_structure`).

- gls:

  (Optional.) Logical; default `TRUE`. When `TRUE`, the design is
  GLS-whitened using REML-estimated (or supplied) variance components —
  the standard, efficient path. When `FALSE`,
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  **skips GLS whitening and variance-component estimation entirely**,
  fitting on the un-whitened (fusion-transformed) design. This is the
  high-dimensional (`p >= NT`) path: REML cannot estimate the variance
  components there (the `p < N(T - 1)` REML guard would otherwise stop
  the fit), and whitening buys efficiency, not validity — the
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  cluster-robust sandwich standard error needs no `Omega` (paper
  Decision D1). A `gls = FALSE` fit has `calc_ses = FALSE` (no
  within-selection oracle standard errors) and un-whitened
  `internal$X_final` / `internal$y_final`; pass it to
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  for a valid cluster-robust SE. `add_ridge = TRUE` is not supported,
  and supplied `sig_eps_sq` / `sig_eps_c_sq` are ignored, under
  `gls = FALSE`.

## Value

An object of class `fetwfe` containing the following elements:

- att_hat:

  The estimated overall average treatment effect for a randomly selected
  treated unit.

- att_se:

  If `q < 1`, a standard error for the ATT. Under the default
  `se_type = "default"`, the SE is the tight Gaussian variance
  `sqrt(att_var_1 + att_var_2)` (Theorem (c\$'\$) under Assumption
  (Psi-IF); paper line 1233 onwards). Assumption (Psi-IF) is satisfied
  by the package's default cohort sample-proportions estimator
  `hat_pi_g = N_g / N` (and by multinomial logit, any GLM on `W | X`,
  and kernel/series regression of `1{W = g}` on `X`), so the default SE
  is asymptotically exact for the package's default estimator. Under
  `se_type = "conservative"` (or in version \<= 1.11.7 by default), the
  SE is the Cauchy-Schwarz upper bound
  `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` from
  Theorem (c). When `indep_counts` is provided, the two-sample exact
  formula `sqrt(att_var_1 + att_var_2)` is used regardless of `se_type`.
  If `q >= 1`, this will be NA.

- att_p_value:

  A two-sided p-value for the overall ATT against the null
  `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
  `att_se` is zero or `NA` (e.g., under the bridge solver's selected-out
  fallback). See the package vignette section "Testing the zero-effect
  null" for interpretation guidance under selection consistency.

- att_selected:

  Logical scalar; `TRUE` if `att_hat` is not exactly zero (i.e., at
  least one cohort's bridge-penalized coefficient survived selection),
  `FALSE` otherwise. Under FETWFE Theorem 6.2 (restriction selection
  consistency), `att_selected = FALSE` is the asymptotic statement that
  the truth is zero. For ridge (`q = 2`) the bridge solver does not zero
  coefficients, so this will typically be `TRUE`.

- catt_hats:

  A named vector containing the estimated average treatment effects for
  each cohort.

- catt_ses:

  If `q < 1`, a named vector containing the (asymptotically exact,
  non-conservative) standard errors for the estimated average treatment
  effects within each cohort.

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
  (`ci_low`, `ci_high`), per-cohort p-values (`p_value`), and a
  `selected` logical flag (`TRUE` when the bridge penalty left the
  cohort's CATT nonzero). For selected-out cohorts (`selected = FALSE`),
  `p_value` is `NA` — the inferential content lives in `selected`. The
  `catt_df` S3 class makes `[[` / `$` / `[` access on the pre-1.11.0
  Title-Case column names (`Cohort`, `Estimated TE`, `SE`, `ConfIntLow`,
  `ConfIntHigh`, `P_value`) [`stop()`](https://rdrr.io/r/base/stop.html)
  with a migration message pointing to the new name. See `NEWS.md` for
  the rename table.

- beta_hat:

  The full vector of estimated coefficients.

- treat_inds:

  The indices of `beta_hat` corresponding to the treatment effects for
  each cohort at each time.

- treat_int_inds:

  The indices of `beta_hat` corresponding to the interactions between
  the treatment effects for each cohort at each time and the covariates.

- sig_eps_sq:

  Either the provided `sig_eps_sq` or the estimated one, if a value
  wasn't provided.

- sig_eps_c_sq:

  Either the provided `sig_eps_c_sq` or the estimated one, if a value
  wasn't provided.

- lambda.max:

  Either the provided `lambda.max` or the one that was used, if a value
  wasn't provided. (This is returned to help with getting a reasonable
  range of `lambda` values for grid search.)

- lambda.max_model_size:

  The number of selected features (excluding the always-present
  intercept) at `lambda.max` (for `q <= 1`, this will be the smallest
  model size). As mentioned above, for `q <= 1` ideally this value is
  close to 0.

- lambda.min:

  Either the provided `lambda.min` or the one that was used, if a value
  wasn't provided.

- lambda.min_model_size:

  The number of selected features (excluding the always-present
  intercept) at `lambda.min` (for `q <= 1`, this will be the largest
  model size). As mentioned above, for `q <= 1` ideally this value is
  close to `p`.

- lambda_star:

  The value of `lambda` chosen by the selection method recorded in
  `lambda_selection`. If this value is close to `lambda.min` or
  `lambda.max`, that could suggest that the range of `lambda` values
  should be expanded.

- lambda_star_model_size:

  The number of selected features (excluding the always-present
  intercept) in the chosen model. If this value is close to
  `lambda.max_model_size` or `lambda.min_model_size`, that could suggest
  that the range of `lambda` values should be expanded.

- lambda_selection:

  Character scalar; either `"cv"` (10-fold cross-validation on
  `cv.grpreg`; v1.13.0+ default) or `"bic"` (BIC over the `grpreg`
  lambda grid; the prior default). Mirrors the `lambda_selection`
  argument the user passed.

- cv_folds:

  Integer scalar; the `cv_folds` value used when
  `lambda_selection = "cv"`, `NA_integer_` when
  `lambda_selection = "bic"`.

- cv_seed:

  Integer scalar; the seed actually fed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  `cv.grpreg()` was called. Defaults to `as.integer(N * T)` when the
  user did not pass a seed. `NA_integer_` when
  `lambda_selection = "bic"`.

- fusion_structure:

  Character scalar; the `fusion_structure` argument the user passed
  (`"cohort"` or `"event_study"`), recording which fusion-penalty
  differences matrix was used for the treatment effects.

- fusion_matrix:

  The user-supplied custom forward differences matrix `D_N` (a
  `num_treats x num_treats` numeric matrix), or `NULL` if none was
  supplied. When non-`NULL` it overrode `fusion_structure` for the
  treatment-effect block; the estimator used `solve(fusion_matrix)`
  internally. See the `fusion_matrix` argument.

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

- alpha:

  The alpha level used for confidence intervals.

- calc_ses:

  Logical indicating whether standard errors were calculated. Same as
  `$internal$calc_ses`; duplicated at top level for parity with
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md),
  and
  [`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
  (#180).

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

- ci_type:

  Character scalar; the `ci_type` argument the user passed
  (`"simultaneous"` or `"pointwise"`), controlling whether the reported
  `catt_df` /
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
  confidence-interval bounds are simultaneous (family-wise) or
  pointwise.

- y_mean:

  Numeric scalar; the mean of the original (pre-centering) response.
  Stored so downstream methods (`augment()`,
  [`predict()`](https://rdrr.io/r/stats/predict.html)) can return fitted
  values on the original-response scale.

- response_col_name:

  Character scalar; the name of the response column in the original
  `pdata`. Consumed by `augment.<class>()`.

- time_var, unit_var, treatment:

  Character scalars; the `time_var` / `unit_var` / `treatment` arguments
  the user passed. Consumed by `augment.<class>()` when auto-aligning a
  user-supplied panel to the fitted design (e.g., dropping
  first-period-treated units the estimator removed internally, and
  sorting rows to match the design matrix's internal `(unit, time)`
  order).

- covs:

  Character vector; the original `covs` argument the user passed (before
  any factor expansion the estimator performed internally). Consumed by
  `augment.<class>()`.

- internal:

  A list containing internal outputs that are typically not needed for
  interpretation:

  X_ints

  :   The design matrix created containing all interactions, time and
      cohort dummies, etc.

  y

  :   The vector of responses, containing `nrow(X_ints)` entries.

  X_final

  :   The design matrix after applying the change in coordinates to fit
      the model and also multiplying on the left by the square root
      inverse of the estimated covariance matrix for each unit.

  y_final

  :   The final response after multiplying on the left by the square
      root inverse of the estimated covariance matrix for each unit.

  theta_hat

  :   The vector of estimated coefficients in the transformed (fused)
      space, including the intercept as the first element.

  calc_ses

  :   Logical indicating whether standard errors were calculated.

  variance_components

  :   A list exposing the two variance pieces (`att_var_1`, `att_var_2`)
      plus their paper-notation counterparts (`V_1`, `V_2`) and the
      unit-scaled variance estimators (`tilde_v_N`, `hat_v_N`,
      `tilde_v_N_C`, `tilde_v_N_C_pi_hat`, `tilde_v_N_C_pi_hat_cons`,
      `tilde_v_N_cons`) catalogued at paper line 2006. The Wald CI is
      `[hat_T_N +- qnorm(1-alpha/2) * sqrt(tilde_v_N / N)]` (paper Eq.
      `conf.int.form`). New in v1.12.0 (issue \#141 + \#146).

  first_year

  :   Integer or numeric scalar; the first (earliest) `time_var` value
      in the panel after `idCohorts()` processing. Consumed by
      [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
      to map `cohort_probs`' cohort labels (treatment-start years) to
      1-based panel-time-index offsets when the labels are
      integer-coercible. New in v1.13.3 (issue \#174).

  d_inv_treat

  :   The inverted custom treatment-effect fusion block
      `solve(fusion_matrix)` (a `num_treats x num_treats` numeric
      matrix), or `NULL` if no `fusion_matrix` was supplied. Consumed by
      [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
      and
      [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
      so the access-time bands reuse the same fusion block the fit used
      (#236).

The object has methods for
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`coef()`](https://rdrr.io/r/stats/coef.html). By default,
[`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) only show the
essential outputs. To see internal details, use
`print(x, show_internal = TRUE)` or `summary(x, show_internal = TRUE)`.
The [`coef()`](https://rdrr.io/r/stats/coef.html) method returns the
vector of estimated coefficients (`beta_hat`).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate coefficients
  coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)

  # Simulate data using the coefficients
  sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)

  result <- fetwfeWithSimulatedData(sim_data)
} # }
```
