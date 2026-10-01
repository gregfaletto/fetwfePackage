# Fused extended two-way fixed effects

Implementation of fused extended two-way fixed effects. Estimates
overall ATT as well as CATT (cohort average treatment effects on the
treated units).

The treatment-effect fusion penalty defaults to a within-/between-cohort
geometry (`fusion_structure = "cohort"`) and also supports an
event-study geometry (`fusion_structure = "event_study"`, fusing effects
at the same time since treatment across cohorts) or a fully custom
`fusion_matrix`. See the `fusion_structure` / `fusion_matrix` arguments
below and
[`vignette("fusion_structure_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/fusion_structure_vignette.md)
for guidance on choosing.

## Usage

``` r
fetwfe(
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
  ci_type = c("simultaneous", "pointwise"),
  fusion_structure = c("cohort", "event_study"),
  fusion_matrix = NULL,
  gls = TRUE
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
  2015; Patterson & Thompson 1971). When supplied, the value also sets
  the scale of the bridge penalty's grid (the response is divided by
  `sqrt(sig_eps_sq)` before the fit), so at `q != 1` a value supplied in
  the wrong units moves the estimates, not only their standard errors.
  Default is NA.

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

  (Optional.) Numeric. Used only on the BIC route
  (`lambda_selection = "bic"`), which selects `lambda` by BIC over a
  grid. The largest `lambda` in the grid will be `lambda.max`. If no
  `lambda.max` is provided, one will be selected automatically. When
  `q <= 1`, the model will be sparse, and ideally all of the following
  are true at once: the smallest model (the one corresponding to
  `lambda.max`) selects close to 0 features, the largest model (the one
  corresponding to `lambda.min`) selects close to `p` features,
  `nlambda` is large enough so that models are considered at every
  feasible model size, and `nlambda` is small enough so that the
  computation doesn't become infeasible. You may want to manually tweak
  `lambda.max`, `lambda.min`, and `nlambda` to try to achieve these
  goals, particularly if the selected model size is very close to the
  model corresponding to `lambda.max` or `lambda.min`, which could
  indicate that the range of `lambda` values was too narrow or coarse.
  You can use the function outputs `lambda.max_model_size`,
  `lambda.min_model_size`, and `lambda_star_model_size` to try to assess
  this. Default is NA.

- lambda.min:

  (Optional.) Numeric. Used only on the BIC route: the smallest `lambda`
  penalty parameter considered. See the description of `lambda.max` for
  details. Default is NA.

- nlambda:

  (Optional.) Integer. Used only on the BIC route: the total number of
  `lambda` penalty parameters considered. See the description of
  `lambda.max` for details. Default is 100.

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
  every regime tested. See the inference vignette section "Choosing the
  bridge penalty parameter" for details.

- cv_folds:

  Integer; number of folds for the CV path. Ignored when
  `lambda_selection = "bic"`. Default is 10.

- cv_seed:

  Integer or `NULL`; the seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) immediately before
  the `cv.grpreg()` call, controlling fold assignment. If supplied, must
  be within `+/- .Machine$integer.max`. If `NULL` (the default), the
  seed defaults internally to `as.integer(N * T)` so consecutive calls
  on the same dataset are reproducible without the user having to
  specify a seed. The seed actually used is stored on the returned
  object as `cv_seed`. Ignored when `lambda_selection = "bic"`.

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
  exact dimension (otherwise `fetwfe()` stops). Under the paper's
  fixed-dimension scoping, *any* finite invertible `D_N` of that
  dimension inherits the paper's inferential guarantees: the theory
  depends on `D_N` only through its invertibility and singular-value
  bounds (Assumption (D) of Faletto 2025), which a fixed invertible
  matrix automatically satisfies, and swapping in a different `D_N` from
  this class changes only constant factors. A numerically near-singular
  (ill-conditioned) `D_N` still yields a valid point estimator but emits
  a [`warning()`](https://rdrr.io/r/base/warning.html) that its inverse
  may be unreliable. Default is `NULL` (use the built-in
  `fusion_structure`).

- gls:

  (Optional.) Logical; default `TRUE`. When `TRUE`, the design is
  GLS-whitened using REML-estimated (or supplied) variance components —
  the standard, efficient path. When `FALSE`, `fetwfe()` **skips GLS
  whitening and variance-component estimation entirely**, fitting on the
  un-whitened (fusion-transformed) design. This is the high-dimensional
  (`p >= NT`) path: REML cannot estimate the variance components there
  (the `p < N(T - 1)` REML guard would otherwise stop the fit), and
  whitening buys efficiency, not validity — the
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

  The largest `lambda` of the grid used, which is a supplied
  `lambda.max` only on the BIC route. (This is returned to help with
  getting a reasonable range of `lambda` values for grid search.)

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

- catt_band_applied:

  Logical scalar; `TRUE` when the fit-time cohort-family simultaneous
  band was computed and written into `catt_df`, and `FALSE` otherwise —
  including when `ci_type = "pointwise"`, when standard errors were
  unavailable, and when the band construction degraded. It records only
  whether the band was applied: it does not imply the band is wider than
  the pointwise interval (they coincide when fewer than two effects have
  positive variance) and does not imply the band is informative (a
  degenerate fit can carry an all-zero band with `TRUE`).

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

## See also

[`vignette("fusion_structure_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/fusion_structure_vignette.md)
for guidance on choosing between the cohort (default) and event-study
fusion penalties and on supplying a custom `fusion_matrix`.

## Author

Gregory Faletto

## Examples

``` r
# `bacondecomp` (which supplies the `divorce` data) is a Suggests-only
# dependency, so guard the example on its availability. The fit is wrapped in
# \donttest{} because it is slower than a toy example.
# \donttest{
if (requireNamespace("bacondecomp", quietly = TRUE)) {
  library(bacondecomp)

  data(divorce)

  # Stevenson & Wolfers (2006): the effect of unilateral ("no-fault") divorce
  # reforms on female suicide rates. Restrict to the female subset
  # (`sex == 2`); `changed` is already the absorbing 0/1 reform indicator, and
  # the elasticity-scaled female suicide rate is the response.
  divorce_f <- divorce[divorce$sex == 2, ]

  # Reproduces the empirical application in Faletto (2025, Sec. 8.2). The 9
  # states already treated by 1964 are auto-dropped as first-period-treated,
  # and `murderrate` is auto-dropped (missing in 1964 for one state); both are
  # reported as (expected) warnings. The noise variances are supplied
  # (precomputed by REML) to keep the example fast and reproducible; the
  # default lambda_selection is "cv" (10-fold cross-validation).
  res <- fetwfe(
      pdata = divorce_f,
      time_var = "year",
      unit_var = "st",
      treatment = "changed",
      covs = c("murderrate", "lnpersinc", "afdcrolls"),
      response = "suiciderate_elast_jag",
      sig_eps_sq = 0.0344,
      sig_eps_c_sq = 0.1507,
      add_ridge = TRUE,
      q = 0.5)

  # FETWFE estimates an overall ATT of roughly -6% on the elasticity-scaled
  # female suicide rate, with a 95% confidence interval that excludes zero.
  # The selection step retains heterogeneous cohort effects (several cohorts
  # are pruned to exactly zero), rather than fusing to a single common effect.
  print(res, max_cohorts = Inf)
}
#> Warning: 9 units were removed because they were treated in the first time period: AK, LA, MD, NC, OK, UT, VA, VT, WV
#> Warning: 1 covariate(s) were removed because they contained missing values in the first time period for at least one unit:  murderrate
#> Fused Extended Two-Way Fixed Effects Results
#> ===========================================
#> 
#> Overall Average Treatment Effect (ATT):
#>   Estimate:   -0.0596
#>   Std. Error: 0.0186
#>   P-value:    0.001327
#>   Selected:   TRUE
#>   95% CI:    [-0.0960, -0.0232]
#> 
#> Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
#>  cohort    estimate          se       ci_low      ci_high      p_value selected
#>    1969  0.00000000 0.000000000  0.000000000  0.000000000           NA    FALSE
#>    1970 -0.44275152 0.046410168 -0.569265971 -0.316237066 0.000000e+00     TRUE
#>    1971 -0.03269001 0.019825552 -0.086734613  0.021354583 5.628621e-01     TRUE
#>    1972 -0.01624466 0.009356597 -0.041750806  0.009261495 4.947428e-01     TRUE
#>    1973 -0.06295968 0.013023899 -0.098462923 -0.027456439 1.069370e-05     TRUE
#>    1974 -0.03051791 0.012975096 -0.065888117  0.004852291 1.392168e-01     TRUE
#>    1975  0.00000000 0.000000000  0.000000000  0.000000000           NA    FALSE
#>    1976 -0.03584431 0.063581218 -0.209167164  0.137478549 9.988446e-01     TRUE
#>    1977 -0.12340880 0.024169900 -0.189296123 -0.057521483 2.634343e-06     TRUE
#>    1980  0.00000000 0.000000000  0.000000000  0.000000000           NA    FALSE
#>    1984  0.00000000 0.000000000  0.000000000  0.000000000           NA    FALSE
#>    1985  0.14777033 0.050861041  0.009122763  0.286417890 2.890215e-02     TRUE
#> 
#> Event-Study Average Treatment Effects (per event time) [simultaneous 95% CI]:
#>  event_time n_cohorts      estimate          se       ci_low    ci_high
#>           0        12  0.0000000000 0.000000000  0.000000000 0.00000000
#>           1        12  0.0116460157 0.012611495 -0.022434196 0.04572623
#>           2        12  0.0038910962 0.004665326 -0.008716077 0.01649827
#>           3        12 -0.0068098902 0.009707715 -0.033043177 0.01942340
#>           4        12 -0.0008192013 0.011643412 -0.032283349 0.03064495
#>           5        12 -0.0008192013 0.011643412 -0.032283349 0.03064495
#>           6        12 -0.0008192013 0.011643412 -0.032283349 0.03064495
#>           7        12 -0.0156592321 0.012529765 -0.049518582 0.01820012
#>           8        12 -0.0484842687 0.022567145 -0.109467766 0.01249923
#>           9        12 -0.0484842687 0.022567145 -0.109467766 0.01249923
#>    p_value
#>         NA
#>  0.9183921
#>  0.9508661
#>  0.9817454
#>  1.0000000
#>  1.0000000
#>  1.0000000
#>  0.7320385
#>  0.1852246
#>  0.1859000
#>   ... and 22 more event times.
#> 
#> Model Details:
#>   Units (N)           : 42
#>   Time periods (T)    : 33
#>   Treated cohorts (G) : 12
#>   Covariates (d)      : 2
#>   Features (p)        : 908
#>   Selected size       : 34
#>   Lambda*             : 0.0004
# }
```
