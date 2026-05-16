# NEWS

## Version 1.9.4 (2026-05-16)

- Fixed a latent bug in `augment.fetwfe()` / `augment.etwfe()` /
  `augment.betwfe()` (`R/broom_methods.R::.augment_estimator_output()`):
  when a user passed an un-trimmed panel that contained first-period-
  treated units (so the augment auto-trim path fired), the returned
  data frame silently dropped the `treatment` column. This was caused
  by augment assigning `data <- ret$df` from `idCohorts()`, which
  strips `treat_var` from its output (`R/utility.R:149`). The no-trim
  path (pre-trimmed input where `nrow(data) == nrow(X_ints)`) was
  unaffected — it preserved every user column. The fix keeps the
  unit-level filtering that `idCohorts()` performs but identifies the
  surviving units from `idCohorts()`'s output and filters the user's
  original `data` directly, so the auto-trim path now preserves the
  `treatment` column and any other extra columns the user supplied
  (e.g., state names, panel IDs) — matching the no-trim path. Point
  estimates, standard errors, and the no-trim path are unchanged.
  Resolves GitHub #51.

## Version 1.9.3 (2026-05-16)

- Fixed a latent bug in two internal helpers — `.truncate_catt()`
  (`R/class_helpers.R`) and `.tidy_estimator_output()`
  (`R/broom_methods.R`) — both of which sorted `catt_df` rows via
  `order(catt_df$Cohort)`. The `Cohort` column is character, so R's
  default `order()` sorted lexicographically: `"10"` and `"11"` came
  before `"2"` and `"3"`. Two visible consequences on panels with
  ≥ 10 unique cohort adoption times: (a) `tidy.<class>()` row order was
  alphabetical rather than chronological; (b) `print.<class>()` and
  `summary.<class>()` (which call `.truncate_catt()` with
  `max_cohorts`) silently dropped the wrong cohorts — the *earliest*
  ones, treated as "later" alphabetically. Both helpers now use a
  composite key `order(suppressWarnings(as.numeric(.)), .)` that
  sorts numerically with a character tiebreak. Behavior on panels with
  ≤ 9 cohorts (every existing test fixture) is bit-identical to
  pre-fix. Resolves GitHub #53.

## Version 1.9.2 (2026-05-16)

- Fixed a bug in `betwfe()`'s standard-error calculation under partial
  selection. The internal `getPsiRUnfused()` helper was normalizing the
  weight vector `psi_r` by `k_sel` (the count of *selected* coefficients
  in a cohort's treatment block) rather than `k_full` (the full block
  size). The cohort point estimate
  `cohort_tes[r] = mean(tes[first_ind_r:last_ind_r])` averages over the
  full block (unselected entries are exact zeros post-bridge), so the
  variance formula's `psi_r` weighting targeted the wrong estimand
  whenever the bridge solver zeroed out some-but-not-all of a cohort's
  coefficients. The bug violated BETWFE's documented "asymptotically
  exact, non-conservative" SE promise. Effect on user-visible output: on
  panels with partial selection, `betwfe()`'s reported `catt_df$SE`,
  `catt_df$ConfIntLow`, `catt_df$ConfIntHigh`, and `catt_df$P_value`
  shift — per-cohort SEs shrink by factor `k_full / k_sel`. The overall
  `att_se` direction is counterintuitive: it can *grow*, because the
  corrected per-cohort CATT vector is more dispersed across cohorts than
  the pre-fix mean-over-selected-only version, and the `att_var_2`
  component is a quadratic form in cohort dispersion. Empirical
  magnitudes are ~2× shifts in either direction on representative
  partial-selection panels. Point estimates (`att_hat`,
  `catt_df$Estimated TE`) are unaffected. FETWFE, ETWFE, and `twfeCovs()`
  are unaffected (FETWFE uses a different `psi_r` constructor; ETWFE and
  twfeCovs never partial-select). Users with `se_type = "cluster"` also
  see SEs shift on partial-selection BETWFE fits (the cluster-robust
  path uses the same `psi_r` weighting).

## Version 1.9.1 (2026-05-15)

- Fixed two coupled bugs in the internal variance-component estimator
  `estOmegaSqrtInv()` (`R/core_funcs.R`). The previous implementation used
  the within-estimator procedure of Pesaran (2015, §26.5.1) with: (a) a
  wrong-LHS in the alpha-hat residual formula that collapsed it to a
  constant across units, and (b) a copy-on-modify trap in the in-place
  X-demeaning loop that left `glmnet` fitting on un-demeaned X.
  Combined, these caused `sig_eps_c_sq_hat` to return as essentially zero
  (~1e-32) regardless of truth, which degenerated the GLS weighting step
  to the identity. Model-based standard errors for `fetwfe()` / `etwfe()` /
  `betwfe()` therefore treated clustered errors as i.i.d. whenever the
  user did not pre-supply `sig_eps_c_sq`. The estimator has been replaced
  by REML on the linear mixed-effects model `y ~ X + (1 | unit)` via
  `lme4::lmer` (Bates, Maechler, Bolker, & Walker 2015; Patterson &
  Thompson 1971; Pinheiro & Bates 2000). REML is the standard procedure
  for variance components in random-intercept models, handles FETWFE's
  design (with its time-invariant columns) natively, and aligns the
  implementation with the random-effects GLS framework the paper assumes.

  `lme4` has been added to `Suggests:` (not `Imports:`) and is gated by
  `requireNamespace("lme4", quietly = TRUE)`. When the package is not
  installed and the user has not supplied `sig_eps_sq` / `sig_eps_c_sq`
  directly, `fetwfe()` (and siblings) error with a clear message pointing
  at `install.packages("lme4")` or supplying the variance components
  manually.

  **User-visible impact: both standard errors and point estimates will
  shift after upgrade.** Reported SE values for `fetwfe()` / `etwfe()` /
  `betwfe()` will be larger on panels with substantial unit-level
  variance, reflecting clustered noise rather than i.i.d. noise. Point
  estimates will also shift, because the GLS weighting step --- which
  was effectively a no-op under the bug --- now actually re-weights the
  design matrix before the bridge regression. On the package's flagship
  divorce-data example (`vignettes/vignette.Rmd`), the overall ATT
  changes from ~-1.84% to 0.00% with every cohort selected out under
  the corrected estimator (the previous value was an artifact of the bug
  preserving cohort 1970's between-unit jump that should have been
  absorbed by the random intercept). The vignette has been updated
  accordingly. The paper's empirical section will be updated in a
  forthcoming revision to the arXiv preprint.

  Users who supplied `sig_eps_sq` and `sig_eps_c_sq` manually are
  unaffected. Users using `se_type = "cluster"` get unchanged
  cluster-robust SEs (the cluster path doesn't depend on
  `estOmegaSqrtInv`), but their point estimates still shift since the
  GLS-transformed regression is shared across SE methods.

  Resolves the bug surfaced in the May 2026 internal code review.

## Version 1.9.0 (2026-05-15)

- Added broom-package S3 methods `tidy()` / `glance()` / `augment()` for
  `fetwfe()` / `etwfe()` / `betwfe()` outputs, plus `tidy()` for the
  outputs of `event_study()` and `getTes()`, so users can pipe estimator
  output directly into `ggplot2` / `modelsummary` / `gt` workflows
  without reading the package class documentation. `tidy()` returns a
  long data frame with broom-standard columns (`term`, `estimate`,
  `std.error`, `statistic`, `p.value`, `conf.low`, `conf.high`; plus
  `selected` for `fetwfe` / `betwfe`), one row for the overall ATT and
  one per cohort. `glance()` returns a one-row model-level summary (13
  columns for `fetwfe` / `betwfe`, 11 for `etwfe` without the
  `lambda_star*` columns). `augment(x, data)` appends `.fitted` and
  `.resid` to the user-supplied panel and auto-trims first-period-
  treated units (the same drop the estimator applies internally during
  fitting), so the call works with either the original `pdata` passed
  to the estimator or a pre-trimmed panel. `tidy.fetwfe_event_study()`
  accepts `conf.int` / `conf.level` for CI control. The estimator
  outputs gain six additive metadata slots — `y_mean`,
  `response_col_name`, `time_var`, `unit_var`, `treatment`, `covs` —
  populated at fit time and consumed by `augment()` (and a future
  `predict()`) so users don't have to re-pass any of them. `getTes()`
  output gains a `cohort_times` slot (the simulator's adoption-time
  labels, `2..(R+1)`) so `tidy.FETWFE_tes()` labels cohorts the same
  way `tidy.<estimator>()` labels them on a fitted panel. Methods are
  registered via the [generics](https://cran.r-project.org/package=generics)
  package (added to `Imports:`); `broom` joins `Suggests:` for the
  vignette demo. Resolves #27.

## Version 1.8.0 (2026-05-13)

- Fixed an off-by-one index in the Jacobian construction inside the
  internal helper `getSecondVarTermDataApp()`, which backs FETWFE's
  `att_var_2` — the cohort-probability variance component of the
  overall ATT standard error. The corrected formula uses the inner-
  loop / column-index marginal cohort probability per paper
  Theorem 6.3 (`paper_arxiv.tex:2577-2592`). The same fix is applied
  to the parallel `.event_study_var2_fetwfe()` helper introduced in
  1.7.0. **The reported `att_se` for `fetwfe()` shifts numerically**;
  the bug could either over- or under-estimate `att_var_2` depending
  on the joint structure of cohort probabilities and per-cohort ATT
  magnitudes. Observed magnitudes: ~30% underestimate on small
  simulated panels with unequal cohorts (Test 2 of the new test
  file); ~13% overestimate on the divorce-data example in the main
  vignette (which has `R = 12` cohorts with very unequal sizes); a
  few percent on near-uniform panels. The shift also propagates
  through `att_p_value` (computed from `att_se`) and through the
  event-study output (`event_study()` and `plot.fetwfe()` columns
  `se`, `ci_low`, `ci_high`, `p_value`). Empirical Monte Carlo
  validation against multinomial resampling of cohort probabilities
  confirms the corrected formula matches the true asymptotic variance
  to within MC noise (under 1% at 10k draws). The bug was masked in
  prior simulation studies because `simulateData()` uses a uniform
  cohort-probability prior, and the buggy and correct formulas are
  algebraically identical when all marginal cohort probabilities are
  equal. `etwfe()`, `betwfe()`, and `twfeCovs()` are unaffected
  (their parallel `getSecondVarTermOLS()` did not have this bug).
  FETWFE cohort-level CATT SEs in `catt_df` are unaffected (they use
  only the regression-coefficient variance, not the cohort-
  probability variance). Closes #46.

## Version 1.7.0 (2026-05-13)

- Added an exported `event_study(x, alpha)` function and S3 methods
  `plot.fetwfe()`, `plot.etwfe()`, `plot.betwfe()` that compute and
  visualize pooled event-time treatment-effect estimates with
  confidence intervals. Pooling weights are sample-cohort-size weights
  (matching `did::aggte(type = "dynamic")` convention); the variance
  combines a regression-coefficient term and a cohort-probability term,
  dispatching on `se_type` for the regression term (model-based under
  default, cluster-robust sandwich under `se_type = "cluster"`) and on
  `indep_counts` for the conservative-vs-asymptotically-exact
  combination. Only post-treatment event times (`e >= 0`) are surfaced
  in this initial release.
- Added `ggplot2` to `Suggests:` (used by the new `plot.*` methods via
  `requireNamespace()` guards; the package still installs and the
  estimators still run without `ggplot2`).
- Added the slot `cohort_probs_overall` (length-`R` numeric vector of
  marginal cohort probabilities, using the indep-sample variant when
  `indep_counts` was provided, else the in-sample variant) to the S3
  outputs of `fetwfe()`, `etwfe()`, `betwfe()` and the list output of
  `twfeCovs()`. Also added `indep_counts_used` (logical scalar) to the
  same four outputs. The `fetwfe()` S3 output additionally gains
  `theta_hat` in `$internal` (the full fused-coordinate coefficient
  vector with intercept). These slots are required by the new
  `event_study()` machinery and are useful for downstream tooling.

## Version 1.6.1 (2026-05-12)

- Added a no-covariate side-by-side demonstration of `etwfe()` and
  `fetwfe()` to `vignettes/etwfe_betwfe_vignette.Rmd`. The main
  vignette `vignettes/vignette.Rmd` now points readers at that
  section near the divorce-data application. An explicit
  `covs = c()` end-to-end test of the public `fetwfe()` entry point
  was added to `tests/testthat/test-fetwfe.R`. The package already
  supported the no-covariate case (the simulator's `d = 0` path is
  exercised by `tests/testthat/test-genCoefs.R` and
  `tests/testthat/test-simulateData.R`); this release closes the
  showcase gap.

## Version 1.6.0 (2026-05-12)

- Added an experimental `se_type = c("default", "cluster")` argument to
  `fetwfe()`, `betwfe()`, `etwfe()`, `twfeCovs()` and their
  `*WithSimulatedData()` wrappers. Setting `se_type = "cluster"`
  replaces the package's Assumption-F1-based standard errors with a
  unit-clustered Liang-Zeger sandwich on the bridge-selected support
  (for FETWFE/BETWFE) or the OLS-selected support (for ETWFE/twfeCovs);
  both `att_se` and `catt_ses` (and the derived confidence intervals
  and p-values) are recomputed from this sandwich. For FETWFE/BETWFE,
  `se_type = "cluster"` is meaningful only when `q < 1` (the bridge
  oracle property is required); for `q >= 1` the SE remains `NA` as
  before. The new option is experimental and the underlying theory is
  pending; see the companion vignette `inference_vignette` for the
  formula, the assumptions it relaxes from F1, and the theory-pending
  caveat. Default is `"default"`, so existing callers see no behavior
  change. The `print()` and `summary()` methods now label the SE line
  "Std. Error (cluster-robust)" / "SE (cluster-robust)" when
  `se_type = "cluster"` was used.

## Version 1.5.6 (2026-05-11)

- The four estimators (`fetwfe()`, `betwfe()`, `etwfe()`, `twfeCovs()`)
  and their `*WithSimulatedData()` wrappers now accept a new
  `allow_no_never_treated = TRUE` argument. When the input panel has
  zero never-treated units (previously a hard error), the package
  auto-truncates the panel by dropping time periods at and after the
  latest cohort's start time --- the units in that latest cohort then
  serve as the never-treated comparison group in the retained
  sub-panel --- and issues a clear warning naming the dropped periods.
  Pass `allow_no_never_treated = FALSE` to restore the previous
  hard-error behavior. The argument has no effect when the input
  already contains never-treated units.

## Version 1.5.5 (2026-05-11)

- `fetwfe()`, `betwfe()`, `etwfe()`, and `twfeCovs()` now surface a
  per-cohort `P_value` column in `catt_df` and a top-level
  `att_p_value` field on the returned list. For `fetwfe()` and
  `betwfe()` (which perform bridge regression with selection), the
  `catt_df` also gains a `selected` logical column and the return
  list gains a top-level `att_selected` field; selected-out cohorts
  have `P_value = NA` and `selected = FALSE`. The previously
  documented but non-functional `order_by = "pvalue"` option in the
  S3 `print()` methods now works as advertised. See the new
  "Testing the zero-effect null" section in the package vignette
  for interpretation guidance under FETWFE's restriction selection
  consistency.

## Version 1.5.4 (2026-05-10)

- Added a new vignette `vignettes/etwfe_betwfe_vignette.Rmd`
  ("Comparison Estimators: ETWFE and BETWFE") demonstrating
  `etwfe()` and `betwfe()` on simulated panel data (with
  comparison to true treatment effects via `getTes()`), plus
  practical guidance on when to use each estimator vs. the
  recommended `fetwfe()`. Resolves #17.
- Output of `betwfe()` is now an S3 object of class `betwfe`
  with `print()`, `summary()`, and `coef()` methods, matching
  the pattern already in place for `fetwfe()` and `etwfe()`.
  Existing field accessors (`$att_hat`, `$att_se`, `$catt_df`,
  etc.) are unchanged.

## Version 1.5.3 (2026-05-10)

- Resolved 11 "Could not resolve link to topic" warnings from
  `devtools::document()` against `R/fetwfe_core.R` and `R/etwfe_core.R`
  by replacing internal-helper Markdown link forms `[func()]` with
  `\code{func()}`, plus one LaTeX-in-`@return` fix to use `\eqn{}`
  instead of `\(...\)`. Internal cleanup; no behavior change. Resolves
  #21.

## Version 1.5.2 (2026-05-10)

- Fixed two redirected URLs in the main vignette (`vignettes/vignette.Rmd`)
  that pointed at a bookdown.org domain now 301-redirecting to Posit
  Connect cloud. Resolves #20.

## Version 1.5.1 (2026-05-10)

- Output of `getTes()` is now an S3 object of class `FETWFE_tes` with
  `print()` and `summary()` methods. Existing field accessors
  (`$att_true`, `$actual_cohort_tes`) are unchanged.

## Version 1.5.0 (2025-07-01)

- Updates to documentation, various minor tweaks to prepare for CRAN submission.

## Version 1.4.1 (2025-06-06)

- Add argument `guarantee_rank_condition` to `simulateData()` to provide option to ensure that final design matrix has full column rank.

## Version 1.4.0 (2025-06-06)

- Add S3 class for `fetwfe()`.

## Version 1.3.2 (2025-06-01)

- Correct some math details, update documentation.

## Version 1.3.1 (2025-06-01)

- Fix some bugs, update documentation.

## Version 1.3.0 (2025-05-25)

- Add functions:
  * `twfeCovs()` (implementation of two-way fixed effects with covariates)
  * `twfeCovsWithSimulatedData()` (analogous to `fetwfeWithSimulatedData()`)

## Version 1.2.0 (2025-05-24)

- Add functions:
  * `btwfe()` (implementation of bridge-penalized extended two-way fixed effects)
  * `btwfeWithSimulatedData()` (analogous to `fetwfeWithSimulatedData()`)

## Version 1.1.0 (2025-05-23)

- Add functions:
  * `etwfe()` (implementation of extended two-way fixed effects)
  * `etwfeWithSimulatedData()` (analogous to `fetwfeWithSimulatedData()`)
  * `attgtToFetwfeDf()` (converts data.frame formatted for `did::att_gt()` to format required for `fetwfe()` and `fewtfe::etwfe()`
  * `ewtfeToFetwfeDf()` (converts data.frame formatted for `etwfe::etwfe()` to format required for `fetwfe()` and `fetwfe::etwfe()`

## Version 1.0.0 (2025-05-14)

- Update vignettes.

## Version 0.11.2 (2025-05-14)

- Update documentation.

## Version 0.11.1 (2025-05-13)

- Update documentation.

## Version 0.11.0 (2025-05-11)

- Various minor bug fixes and refactoring.

## Version 0.10.4 (2025-05-08)

- Fix issue #1 [https://github.com/gregfaletto/fetwfePackage/issues/1](https://github.com/gregfaletto/fetwfePackage/issues/1).

## Version 0.10.3 (2025-05-08)

- Fix bug from previous centering and scaling implementation.

## Version 0.10.2 (2025-05-07)

- Add centering and scaling of covariates before ridge penalty added.

## Version 0.10.1 (2025-03-08)

- Add support for factor covariates in `fetwfe()`.

## Version 0.10.0 (2025-03-04)

- Add back in functions `genCoefsCore()` and `simulateDataCore()` for export, and add tests for them. (Though these more flexible functions are now available, emphasis for end users will still be on using wrapper functions `genCoefs()` and `simulateData()`.)

## Version 0.9.0 (2025-03-04)

- Changed implementation of `genCoefs()` to a wrapper function that creates an object that can be piped into `simulateData()`, itself a new wrapper function for the function that was previously called `genRandomData()`. The output of `simulateData()` can then be piped into another new wrapper function, `fetwfeWithSimulatedData()`.

## Version 0.8.0 (2025-03-02)

- Added argument `add_ridge` to `fetwfe()`, which adds a ridge regularization term to the (untransformed) coefficients in estimation, similarly to the elastic net.

## Version 0.7.0 (2025-03-01)

- Added function `getTes()` to get treatment effects from a vector of coefficients generated by `genCoefs()` and/or estimated on data generated from `genRandomData()`.

## Version 0.6.0 (2025-02-28)

- Added functions `genRandomData()`, for generating random panel data suitable for `fetwfe()`, and `genCoefs()`, which generates a random vector of coefficients that is needed as an input to `genRandomData()`.

## Version 0.5.1 (2025-02-26)

- Slight updates to handle cases when only one column or row of matrix is selected more smoothly.

## Version 0.5.0 (2025-02-25)

- Add support for data with no covariates.

## Version 0.4.5 (2025-02-22)

- Fix typo in definition of `lambda.max_model_size`

## Version 0.4.4 (2025-02-20)

- Update DESCRIPTION to align with CRAN requirements.

## Version 0.4.3 (2025-02-19)

- Change URL to DOI in DESCRIPTION.

## Version 0.4.2 (2025-02-16)

- Added more detail to DESCRIPTION.
- Changed how messages are printed to use `message()` instead of `print()`, and modified code so that all messages are suppressed when the argument `verbose` is FALSE.

## Version 0.4.1 (2025-02-07)

- Modified `.Rbuildignore` to conform with CRAN requirements.

## Version 0.4.0 (2025-02-07)

- Initial release of the `fetwfe` package implementing the Fused Extended Two-Way Fixed Effects estimator.
