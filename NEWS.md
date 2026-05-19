# NEWS

## Version 1.9.16 (2026-05-18)

- Extracted two micro-helpers to `R/utility.R` and standardized one
  threshold convention (GitHub #83):
  (1) `.cohort_block_inds(r, R, first_inds, num_treats)` replaces 6
  inline `first_ind_r <- first_inds[r]; last_ind_r <- if (r < R) ...
  else num_treats` constructions across `R/core_funcs.R`,
  `R/fetwfe_core.R`, `R/gen_funcs.R`, and `R/ols_calcs.R`.
  (2) `.multinomial_cov(probs)` replaces 3 inline `Sigma_pi_hat <-
  -outer(probs, probs); diag(...) <- probs * (1 - probs)` blocks
  across `R/fetwfe_core.R`, `R/ols_calcs.R`, and `R/event_study.R`.
  (3) `R/core_funcs.R`'s 4 occurrences of `10^(-6)` replaced with
  `1e-6` to match the rest of the package (bit-identical value;
  pure cosmetic). No user-visible behavior change.

## Version 1.9.15 (2026-05-18)

- Consolidated the cluster-robust sandwich assembly across 4 call
  sites (`getCohortATTsFinalOLS`, `getCohortATTsFinal`,
  `.event_study_etwfe_betwfe`, `.event_study_fetwfe`) into a single
  internal helper `.assemble_cluster_robust_sandwich()` in
  `R/ols_calcs.R`. No user-visible behavior change; the sandwich
  matrix and `treat_block_mask` outputs are bit-identical to v1.9.14.
  The consolidation reduces the drift surface that produced the
  v1.9.5 `10e-6` threshold typo (#56) from 4 sites to 1.
  Resolves #78.

## Version 1.9.14 (2026-05-18)

- Refactored the four "no features selected" early-exit blocks
  in `fetwfe()` and `betwfe()` (intercept-only and no-treatment
  cases for each estimator) to share a single internal helper
  `.build_selected_out_result()` in `R/core_funcs.R`. No
  user-visible behavior change: the return-list shape, field
  ordering, and values are byte-identical to v1.9.13. The
  consolidation eliminates a 4-way drift surface that produced
  9 stale `q < 1` references in v1.9.5 (#56). Resolves #80.

## Version 1.9.13 (2026-05-18)

- Tightened `idCohorts()`'s panel-balance check to catch units
  with the right row count but duplicate or missing time periods
  (GitHub #75). A unit with rows at times `(1, 2, 2)` and `T = 3`
  previously passed silently, causing a confusing downstream
  failure in `processCovs()`. The new check requires both
  `nrow == T` and exact time coverage; the violation message
  includes both counts. The parallel balance check in
  `processCovs()` (defense-in-depth) is tightened to match.

## Version 1.9.12 (2026-05-18)

- Fixed a silent bug in `betwfe(..., add_ridge = TRUE)`
  (GitHub #74): the L2 ridge augmentation rows were applied
  in the FETWFE fusion-transform basis (`D_inverse`) instead
  of the identity basis BETWFE requires. Default-off, but
  affected users who set `add_ridge = TRUE` to escape
  rank-deficient designs. Also made `is_fetwfe` a required
  argument in the internal `prep_for_etwfe_regression()`
  helper so the same bug class cannot recur silently.

## Version 1.9.11 (2026-05-18)

- Refactored the three `print.<class>` / `summary.<class>` /
  `print.summary.<class>` method bodies (across `fetwfe()`,
  `etwfe()`, `betwfe()`) to share three new internal helpers in
  `R/class_helpers.R`. Output is byte-identical to v1.9.10 (locked
  by the snapshot guardrail added in PR #90). Resolves #77.

## Version 1.9.10 (2026-05-17)

- Fixed a cross-method consistency bug in `event_study()` (GitHub #73):
  the function previously reported finite SEs and p-values for `q >= 1`
  fits even though the corresponding `fetwfe()` / `betwfe()` objects
  correctly return `att_se = NA` in that regime (the bridge oracle
  property is required for the model-based SEs). After this fix,
  `event_study()` propagates the fit's `calc_ses` status: when the fit
  has no valid SEs, `event_study()` returns NA SEs and NA p-values,
  matching the parent object's contract. The point estimates remain
  valid (and finite) in both regimes. Same fix applies under
  `se_type = "cluster"`.
- Added internal method-entry preconditions for cross-class consumers
  (`event_study`, `augment`, `tidy`, `glance`, `plot`, `coef`). Each
  such method now re-validates the input fitted object before computing
  on it, catching hand-modified or programmatically-malformed inputs
  early with a clear error message. Foundation PR beta of three (with
  Foundation alpha #85 already shipped and Foundation gamma #87
  drift-sentinel subagent queued). Resolves GitHub #86.

## Version 1.9.9 (2026-05-17)

- Added internal constructor validators for each estimator class
  (`fetwfe`, `etwfe`, `betwfe`, `twfeCovs`). Each entry-point function
  now runs `.validate_<class>(out)` immediately before class assignment.
  The validators encode the documented cross-slot contracts: slot
  inventory (matches `@return` docs), SE consistency (if `calc_ses`
  is FALSE then `att_se` and `catt_ses` must all be NA), selection
  consistency (`att_selected` iff `att_hat != 0` for the classes that
  have it), p-value NA-derivation, catt_df shape, cohort-probability
  structural sanity, beta_hat / y / X_ints dimensions, lambda
  monotonicity (for the bridge estimators), and type sanity. Resolves
  GitHub #85. This is the first of three foundation PRs (with #86
  method-entry preconditions and #87 drift-sentinel subagent) designed
  to prevent the drift class that produced multiple recent cleanup PRs.
- No public API change. No behavior change on well-formed objects: the
  validators only fire if an object would be malformed (which currently
  cannot happen via the public estimator path; verified empirically
  against q in {0.5, 1, 2} and the all-zero-theta fallback). Shared
  contract-checking helpers live in `R/class_helpers.R`; the per-class
  validators live in `R/<class>_class.R` (or `R/twfeCovs.R` for
  twfeCovs pending a class file in #76).

## Version 1.9.8 (2026-05-16)

- Documentation sweep: refreshed `@return` blocks on `etwfe()`,
  `betwfe()`, `twfeCovs()`, and their `*WithSimulatedData()` wrappers
  to match current slot inventories (slots added since v1.5.5 /
  v1.7.0 / v1.9.0 were missing from the docs). Also added the
  `theta_hat` sub-slot to the `internal` listing for `fetwfe()` and
  `fetwfeWithSimulatedData()`. Added `@examples` blocks to `etwfe()`
  and `twfeCovs()` (the other public estimator entry points already
  had them). Removed copy-paste references to a bridge-penalty `q`
  parameter from `twfeCovs.R` and `etwfe_core.R` docs (neither
  function takes a `q` argument; both are pure OLS). Removed a
  self-referential `@inheritParams getTeResults2` line from
  `getTeResults2`'s own documentation block. Cleaned up `etwfe_core()` and
  `twfeCovs_core()` `@return` blocks, which previously listed
  bridge-only slots (`theta_hat`, `lambda.*`) those functions never
  return. Updated a stale mock value in `test-class-helpers.R`
  (`se_type = "standard"` → `"default"` to mirror the live
  `match.arg()` choice). Resolves GitHub #55.
- `attgtToFetwfeDf()` and `etwfeToFetwfeDf()` gain a `verbose = FALSE`
  argument; the previously-unconditional "Dropped N unit-period(s)
  treated in the first period" message is now gated on `verbose`.
  Default-silenced to match the project convention that
  non-interactive code paths should not emit console output without
  opt-in. Users who want the message back can call with
  `verbose = TRUE`.

## Version 1.9.7 (2026-05-16)

- Internal cleanup pass: dropped four unused parameters from internal
  helpers, removed one block of redundant self-validating assertions,
  and corrected a typo in an internal helper name. No public API
  changes; no behavior changes. Resolves GitHub #59.
- `idCohorts()` (`R/utility.R`) no longer accepts the unused `covs`
  argument. Two internal call sites updated.
- `genTreatVarsSim()` (`R/gen_funcs.R`) no longer accepts the unused
  `d` argument. One internal call site updated.
- `getPsiRUnfused()` (`R/ols_calcs.R`) no longer accepts the unused
  `gram_inv` argument. Three internal assertions on `gram_inv`
  dimensions were redundant with conditions already enforced
  upstream and have been removed. Two source call sites and two
  test call sites updated.
- `checkFetwfeInputs()` (`R/fetwfe_core.R`) no longer accepts the
  unused `nlambda` argument. Three call sites updated. `nlambda`
  itself remains unchanged on the public entry points (`fetwfe()`,
  `betwfe()`, `twfeCovs()`) where it is forwarded to `gBridge`.
- `getFirstInds()` (`R/utility.R`) no longer runs a block of
  assertions that re-derived its own closed-form formula and
  compared the result to itself. The first three assertions on the
  cardinality of `f_inds` are preserved.
- Renamed the internal helper `prep_for_etwfe_regresion` to
  `prep_for_etwfe_regression` (typo correction). The function is
  `@noRd`, not exported, so the rename is internal.

## Version 1.9.6 (2026-05-16)

- `R/utility.R::idCohorts()` now reports every malformed unit at once
  rather than stopping on the first. Previously a user with multiple
  malformed units had to fix one, re-run, hit the next stop, fix
  again, re-run — one error round-trip per bad unit. After this
  change, both the balance check (every unit must appear in exactly
  T periods) and the absorbing-state check (treatment, once 1, must
  stay 1) collect violations across the entire unit loop and produce
  a single grouped error message with all offending units listed
  alphabetically (lexicographically; the unit name column is required
  to be character). Listings over 20 entries per bucket are truncated with
  a trailing summary count. The error message wording starts with
  the canonical prefixes "Panel does not appear to be balanced" and
  "Treatment does not appear to be an absorbing state" so downstream
  grep-based consumers still match. The cohort-level stops at
  `R/utility.R:191-195` ("all units treated in the first period") and
  `:211-214` are end-of-loop validation, not per-unit, and remain
  unchanged. Resolves GitHub #64. This was item 5 of #56, deferred
  from v1.9.5 because it changes user-visible error wording.

## Version 1.9.5 (2026-05-16)

- Small consistency cleanups surfaced by an internal code review;
  bundled as one PR. None of these fire under current package usage,
  but each is wrong on inspection and would bite a future refactor of
  the surrounding machinery. Resolves GitHub #56.
- (Item 1) `R/event_study.R::.assemble_event_study_df()` now calls
  the canonical `.compute_p_values()` helper (`R/utility.R:508`)
  rather than re-implementing the p-value computation inline. The
  previous inline form used `is.na(ses) | ses == 0` to mask NA / zero
  SEs; the canonical helper additionally masks negative SEs. The
  package's own SE machinery never produces negative SEs (everything
  comes out of `sqrt(...)`), so on every code path users currently
  exercise the output is bit-identical. The change is a defensive
  alignment against future refactors of the SE machinery.
- (Item 2) `R/fetwfe_core.R::getSecondVarTermDataApp()` line 2022
  relaxes `stopifnot(length(sel_inds[[r]]) < length(sel_inds[[r-1]]))`
  to `<=`, matching the OLS analogue at `R/ols_calcs.R:449`. Both
  assertions guard the same invariant; under the package's standard
  cohort layout (`first_inds = getFirstInds(R, T)`) block sizes are
  strictly decreasing and both forms pass. Behavior unchanged on
  current usage; the relaxation is a no-op except for any future
  non-standard cohort layout.
- (Item 3) `R/ols_calcs.R::getSecondVarTermOLS()` line 465 changes
  the `Sigma_pi_hat`-degeneracy threshold literal from `10e-6`
  (= 1e-5) to `1e-6`. The previous literal reads naturally as `1e-6`
  to an R reader who hasn't done the arithmetic; the change aligns
  the typed literal with that reading. The assertion is
  `stopifnot(sum(cohort_probs_overall) < 1 - eps)` and the change
  *loosens* the rejection band from width `1e-5` to width `1e-6`.
  The relaxed band is still well above floating-point precision, and
  the guard is unreachable in normal operation: under
  `allow_no_never_treated = TRUE` (default) the panel is auto-truncated
  so `sum(cohort_probs_overall) <= (N-1)/N`, with margin `>= 1/N`
  for `N >= 50`. The symmetric guard
  (`stopifnot(length(cohort_probs_overall) == R)` and
  `stopifnot(sum(cohort_probs_overall) < 1 - 1e-6)`) is added to the
  FETWFE analogue `R/fetwfe_core.R::getSecondVarTermDataApp()`,
  which previously had no such guard.
- (Item 4) `R/utility.R::idCohorts()` panel-balance `stop()` message
  now has balanced parentheses. The previous message read
  `"Panel does not appear to be balanced (unit X does not have
  exactly T observations for T = 3"` — the open parenthesis after
  "balanced" was never closed. Cosmetic.
- One additional item from the same review — `idCohorts()` collects
  all malformed units rather than stopping on the first — is deferred
  to a follow-up because it changes user-visible error wording in a
  way that deserves its own framing.

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
