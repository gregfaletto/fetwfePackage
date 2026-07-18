# Bridge-penalized extended two-way fixed effects

Implementation of extended two-way fixed effects with a bridge penalty.
Estimates overall ATT as well as CATT (cohort average treatment effects
on the treated units).

## Usage

``` r
betwfe(
  pdata,
  time_var,
  unit_var,
  treatment,
  response,
  covs = c(),
  indep_counts = NA,
  sig_eps_sq = NA,
  sig_eps_c_sq = NA,
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
  ci_type = c("simultaneous", "pointwise")
)
```

## Arguments

- pdata:

  Dataframe; the panel data set. Each row should represent an
  observation of a unit at a time. Should contain columns as described
  below.

- time_var:

  Character; the name of a single column containing a variable for the
  time period. This column is expected to contain integer values (for
  example, years). Recommended encodings for dates include format YYYY,
  YYYYMM, or YYYYMMDD, whichever is appropriate for your data.

- unit_var:

  Character; the name of a single column containing a variable for each
  unit. This column is expected to contain character values (i.e. the
  "name" of each unit).

- treatment:

  Character; the name of a single column containing a variable for the
  treatment dummy indicator. This column is expected to contain integer
  values, and in particular, should equal 0 if the unit was untreated at
  that time and 1 otherwise. Treatment should be an absorbing state;
  that is, if unit `i` is treated at time `t`, then it must also be
  treated at all times `t` + 1, ..., `T`. Any units treated in the first
  time period will be removed automatically. Please make sure yourself
  that at least some units remain untreated at the final time period
  ("never-treated units").

- response:

  Character; the name of a single column containing the response for
  each unit at each time. The response must be an integer or numeric
  value.

- covs:

  (Optional.) Either a character vector containing the names of the
  columns for covariates (e.g., `covs = c("x1", "x2")`), or a one-sided
  formula (e.g., `covs = ~ x1 + x2`) – the formula form mirrors the
  convention used by `did::att_gt(xformla = ...)`. Only additive bare
  variable names are supported in the formula form; for derived
  variables, compute them in the data frame first and pass via the
  character-vector form. All of these columns are expected to contain
  integer, numeric, or factor values, and any categorical values will be
  automatically encoded as binary indicators. If no covariates are
  provided, the treatment effect estimation will proceed, but it will
  only be valid under unconditional versions of the parallel trends and
  no anticipation assumptions. Default is c().

- indep_counts:

  (Optional.) Integer; a vector. If you have a sufficiently large number
  of units, you can optionally randomly split your data set in half
  (with `N` units in each data set). The data for half of the units
  should go in the `pdata` argument provided above. For the other `N`
  units, simply provide the counts for how many units appear in the
  untreated cohort plus each of the other `G` cohorts in this argument
  `indep_counts`. The benefit of doing this is that the standard error
  for the average treatment effect will be (asymptotically) exact
  instead of conservative. The length of `indep_counts` must equal 1
  plus the number of treated cohorts in `pdata`. All entries of
  `indep_counts` must be strictly positive (if you are concerned that
  this might not work out, maybe your data set is on the small side and
  it's best to just leave your full data set in `pdata`). The sum of all
  the counts in `indep_counts` must match the total number of units in
  `pdata`. Default is NA (in which case conservative standard errors
  will be calculated if `q < 1`.)

- sig_eps_sq:

  (Optional.) Numeric; the variance of the row-level IID noise assumed
  to apply to each observation. See Section 2 of Faletto (2025) for
  details. It is best to provide this variance if it is known (for
  example, if you are using simulated data). If this variance is
  unknown, this argument can be omitted, and the variance will be
  estimated by REML on the linear mixed-effects model
  `y ~ X + (1 | unit)` via
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) (Bates et al.
  2015; Patterson & Thompson 1971). Default is NA.

- sig_eps_c_sq:

  (Optional.) Numeric; the variance of the unit-level IID noise (random
  effects) assumed to apply to each observation. See Section 2 of
  Faletto (2025) for details. It is best to provide this variance if it
  is known (for example, if you are using simulated data). If this
  variance is unknown, this argument can be omitted, and the variance
  will be estimated by REML via
  [`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) on the linear
  mixed-effects model `y ~ X + (1 | unit)` (Bates et al. 2015; Patterson
  & Thompson 1971). Default is NA.

- lambda.max:

  (Optional.) Numeric. A penalty parameter `lambda` will be selected
  over a grid search by BIC in order to select a single model. The
  largest `lambda` in the grid will be `lambda.max`. If no `lambda.max`
  is provided, one will be selected automatically. When `q <= 1`, the
  model will be sparse, and ideally all of the following are true at
  once: the smallest model (the one corresponding to `lambda.max`)
  selects close to 0 features, the largest model (the one corresponding
  to `lambda.min`) selects close to `p` features, `nlambda` is large
  enough so that models are considered at every feasible model size, and
  `nlambda` is small enough so that the computation doesn't become
  infeasible. You may want to manually tweak `lambda.max`, `lambda.min`,
  and `nlambda` to try to achieve these goals, particularly if the
  selected model size is very close to the model corresponding to
  `lambda.max` or `lambda.min`, which could indicate that the range of
  `lambda` values was too narrow or coarse. You can use the function
  outputs `lambda.max_model_size`, `lambda.min_model_size`, and
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
  regularization. `q` = 1 is the lasso, and for 0 \< `q` \< 1, it is
  possible to get standard errors and confidence intervals. `q` = 2 is
  ridge regression. See Faletto (2025) for details. Default is 0.5.

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
  `"conservative"` returns the Cauchy-Schwarz upper bound from Theorem
  (c); use only if the propensity-score estimator violates (Psi-IF).
  `"cluster"` is an *experimental* unit-clustered Liang-Zeger sandwich
  SE on the bridge-selected support (see the companion vignette
  `inference_vignette` for details). `"cluster"` is only meaningful when
  `q < 1` (the bridge oracle property is required); for `q >= 1` the SE
  will be `NA` regardless of `se_type`. The default `"default"` was the
  conservative Cauchy-Schwarz formula in versions \<= 1.11.7; v1.12.0
  switched the default to the tight Gaussian variance. Default is
  `"default"`.

- lambda_selection:

  Character; method for selecting the bridge penalty parameter `lambda`.
  Either `"cv"` (10-fold cross-validation on `cv.grpreg`; the v1.13.0+
  default) or `"bic"` (BIC over the `grpreg` lambda grid; the prior
  default for v1.12.0 and earlier). The default changed in v1.13.0 to
  address a finite-sample bias issue documented in simulation studies
  (see issue \#164). Pass `lambda_selection = "bic"` to recover the
  prior behavior. See the inference vignette section "Choosing the
  bridge penalty parameter" for details.

- cv_folds:

  Integer; number of folds for the CV path. Ignored when
  `lambda_selection = "bic"`. Default is 10.

- cv_seed:

  Integer or `NULL`; the seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  the `cv.grpreg()` call. If `NULL` (the default), the seed defaults
  internally to `as.integer(N * T)`. Ignored when
  `lambda_selection = "bic"`.

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

## Value

An object of class `betwfe` containing the following elements:

- att_hat:

  The estimated overall average treatment effect for a randomly selected
  treated unit.

- att_se:

  If `q < 1`, a standard error for the ATT. If `indep_counts` was
  provided, this standard error is asymptotically exact; if not, it is
  asymptotically conservative. If `q >= 1`, this will be NA.

- att_p_value:

  A two-sided p-value for the overall ATT against the null
  `H_0: tau = 0`, computed as `2 * pnorm(-|att_hat / att_se|)`. `NA` if
  `att_se` is zero or `NA` (e.g., under the bridge solver's selected-out
  fallback).

- att_selected:

  Logical scalar; `TRUE` if `att_hat` is not exactly zero, `FALSE`
  otherwise. BETWFE uses bridge regression directly on the coefficients
  (rather than on the fused restrictions used by FETWFE); under the
  bridge oracle property of Kock (2013), `att_selected = FALSE` is an
  analogous asymptotic statement that the truth is zero under a sparsity
  assumption different from the one Theorem 6.2 establishes for FETWFE.
  For ridge (`q = 2`) the bridge solver does not zero coefficients, so
  this will typically be `TRUE`.

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
  `p_value` is `NA`. The `catt_df` S3 class makes `[[` / `$` / `[`
  access on the pre-1.11.0 Title-Case column names (`Cohort`,
  `Estimated TE`, `SE`, `ConfIntLow`, `ConfIntHigh`, `P_value`)
  [`stop()`](https://rdrr.io/r/base/stop.html) with a migration message
  pointing to the new name. See `NEWS.md` for the rename table.

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
  intercept) at `lambda.max` (for `q <= 1`, the smallest model). As
  mentioned above, for `q <= 1` ideally this value is close to 0.

- lambda.min:

  Either the provided `lambda.min` or the one that was used, if a value
  wasn't provided.

- lambda.min_model_size:

  The number of selected features (excluding the always-present
  intercept) at `lambda.min` (for `q <= 1`, the largest model). As
  mentioned above, for `q <= 1` ideally this value is close to `p`.

- lambda_star:

  The value of `lambda` chosen by the method recorded in
  `lambda_selection`. If this value is close to `lambda.min` or
  `lambda.max`, that could suggest that the range of `lambda` values
  should be expanded.

- lambda_star_model_size:

  The number of selected features (excluding the always-present
  intercept) in the chosen model. If this value is close to
  `lambda.max_model_size` or `lambda.min_model_size`, that could suggest
  that the range of `lambda` values should be expanded.

- lambda_selection:

  Character scalar; either `"cv"` or `"bic"`. Mirrors the
  `lambda_selection` argument the user passed.

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

- ci_type:

  Character scalar; the `ci_type` argument the user passed
  (`"simultaneous"` or `"pointwise"`), controlling whether the reported
  `catt_df` /
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
  confidence-interval bounds are simultaneous (family-wise) or
  pointwise.

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

- y_mean:

  Numeric scalar; mean of the original (pre-centering) response. Stored
  so downstream methods (`augment()`,
  [`predict()`](https://rdrr.io/r/stats/predict.html)) can return fitted
  values on the original-response scale.

- response_col_name:

  Character scalar; the response column name in the original `pdata`.
  Consumed by
  [`augment.betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/augment.betwfe.md).

- time_var, unit_var, treatment:

  Character scalars; the corresponding arguments the user passed.
  Consumed by
  [`augment.betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/augment.betwfe.md)
  when auto-aligning a user-supplied panel to the fitted design.

- covs:

  Character vector; the original `covs` argument (pre-factor-
  expansion). Consumed by
  [`augment.betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/augment.betwfe.md).

- alpha:

  The alpha level used for confidence intervals.

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

## References

Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions. *arXiv preprint
arXiv:2312.05985*. <https://arxiv.org/abs/2312.05985>.

Bates, D., Maechler, M., Bolker, B., & Walker, S. (2015). Fitting Linear
Mixed-Effects Models Using lme4. *Journal of Statistical Software*,
67(1), 1-48.
[doi:10.18637/jss.v067.i01](https://doi.org/10.18637/jss.v067.i01) .

Patterson, H. D., & Thompson, R. (1971). Recovery of inter-block
information when block sizes are unequal. *Biometrika*, 58(3), 545-554.

Pinheiro, J. C., & Bates, D. M. (2000). *Mixed-Effects Models in S and
S-PLUS*. Springer.

## Author

Gregory Faletto

## Examples

``` r
# `bacondecomp` (which supplies the `castle` data) is a Suggests-only
# dependency, so guard the example on its availability.
if (requireNamespace("bacondecomp", quietly = TRUE)) {
  library(bacondecomp)

  data(castle)

  # Response: the log homicide rate. Treatment: `cdl` records the share of
  # the year the castle-doctrine law was in effect, so `cdl > 0` gives the
  # absorbing 0/1 treatment indicator. No `covs`: castle's smallest
  # adoption cohorts contain a single state, so the design is
  # rank-deficient once any covariate is added.
  castle$l_homicide <- log(castle$homicide)
  castle$treated <- as.integer(castle$cdl > 0)

  # On this panel betwfe's bridge penalty selects every cohort out, so the
  # estimated ATT and cohort effects below are all zero.
  res <- betwfe(
      pdata = castle,
      time_var = "year",
      unit_var = "state",
      treatment = "treated",
      response = "l_homicide",
      verbose = TRUE)

  # Average treatment effect on the treated units (in percentage point
  # units)
  100 * res$att_hat

  # Conservative 95% confidence interval for ATT (in percentage point units)

  low_att <- 100 * (res$att_hat - qnorm(1 - 0.05 / 2) * res$att_se)
  high_att <- 100 * (res$att_hat + qnorm(1 - 0.05 / 2) * res$att_se)

  c(low_att, high_att)

  # Cohort average treatment effects and confidence intervals (in percentage
  # point units)

  catt_df_pct <- res$catt_df
  catt_df_pct[["estimate"]] <- 100 * catt_df_pct[["estimate"]]
  catt_df_pct[["se"]] <- 100 * catt_df_pct[["se"]]
  catt_df_pct[["ci_low"]] <- 100 * catt_df_pct[["ci_low"]]
  catt_df_pct[["ci_high"]] <- 100 * catt_df_pct[["ci_high"]]

  catt_df_pct
}
#> No covariates provided; skipping covariate processing.
#> Getting omega sqrt inverse estimate...
#> Done! Time to estimate noise variances:
#> 0.249144077301025
#> Time to get sqrt inverse matrix:
#> 0.000271320343017578
#> Estimating bridge regression with 10-fold CV...
#> Done! Time for estimation:
#> 0.106780767440796
#> No treatment features selected; all treatment effects estimated to be 0.
#>   cohort estimate se ci_low ci_high p_value selected
#> 1   2005        0  0      0       0      NA    FALSE
#> 2   2006        0  0      0       0      NA    FALSE
#> 3   2007        0  0      0       0      NA    FALSE
#> 4   2008        0  0      0       0      NA    FALSE
#> 5   2009        0  0      0       0      NA    FALSE
```
