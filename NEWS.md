# NEWS

## Version 1.34.0

### New features

- `simultaneousCIs(method = "bootstrap")` now extends the high-dimensional
  (`p >= NT`) **desparsified** path to `family = "event_study"` and **debiases
  the high-dimensional band centers** (#299). A high-dim `event_study` fit now
  uses the full-design `riesz_lasso` regression channel together with the
  per-unit propensity channel `F_pi` (the two-channel `Sigma = Sigma_1 + Sigma_2`
  bootstrap), rather than falling back to the fixed-`p` selected-support
  construction. In the high-dimensional regime all bands are now centered on the
  **debiased** estimate (`estimates + colSums(F_reg)/(NT)`, the Theorem 6.6
  correction, equal to `debiasedATT()`'s point estimate for the matching
  contrast) instead of the post-selection bridge estimate; this also re-centers
  the high-dimensional `cohort` / `all_post_treatment` / `custom` bands added in
  1.32.0. The `p >= NT` path remains **experimental** (coverage not yet
  simulation-validated). Fixed-`p` behavior is unchanged.
- In the high-dimensional (`p >= NT`) regime, a `simultaneousCIs()` fit whose
  bridge penalty zeroes *every* treatment effect now emits a `warning()` rather
  than a silent `message()`: an all-zero band there is degenerate and does not
  imply the effects are zero (selection is not consistent at `p >= NT`, and the
  debiased band center is bypassed; see #304). Fixed-`p` behavior (a `message()`)
  is unchanged.

## Version 1.33.0

### New features

- `simultaneousCIs(method = "bootstrap")` now supports `family = "event_study"`,
  completing the multiplier bootstrap's family coverage (#142). The event-study
  family is the one whose simultaneous-band variance carries a non-zero
  cohort-probability (propensity) term `Sigma_2`, so the bootstrap perturbs a
  per-unit propensity influence function (the cohort-sample-proportion summand
  `xi_i = e_{W_i} - pi_hat`, `pi_hat_g = N_g / N`) with its **own independent**
  multiplier stream alongside the regression channel -- a two-channel bootstrap
  matching the analytic `Sigma = Sigma_1 + Sigma_2` (paper Theorem C.1 Step 5
  shows the two channels are asymptotically uncorrelated, so their variances add
  with no cross-term). On a `se_type = "cluster"` fit the event-study bootstrap
  per-effect standard errors equal the analytic ones to machine precision. The
  fixed-`p` selected-support construction is used (including for high-dimensional
  `p >= NT` fits whose selected support is low-dimensional); the high-dimensional
  desparsified event-study combination remains a planned follow-up.

## Version 1.32.0

### New features

- The `simultaneousCIs()` multiplier bootstrap (`method = "bootstrap"`, added in
  1.31.0) now extends to the **high-dimensional (`p >= NT`) regime**, where the
  analytic Gram inverse need not exist so there is no analytic sup-t band. There
  the bootstrap builds each effect's debiasing direction by the full-design
  **desparsified (nodewise / `riesz_lasso`) construction** of `debiasedATT()`
  (#31, paper Theorem 6.6), generalized to the effect family -- uniformly valid,
  **not** post-selection. It returns per-effect `feasibility` / `converged` /
  `lambda_node` diagnostics and a `regime` field. The path is **experimental**
  (`fetwfe()` fits only; its finite-sample coverage is not yet simulation-
  validated). The penalty constant `lambda_c` defaults to the theory scale `1.0`,
  which keeps the nodewise directions feasible; smaller values can leave them
  infeasible, in which case a warning fires and the affected bands are
  unreliable. New arguments `lambda_c` / `riesz_max_iter` / `riesz_tol` control
  the solver (ignored when `p < NT`). The simulator can now generate the
  `p > NT` panels needed to exercise this path (see 1.30.0 / #293). `event_study`
  and the fixed-`p` path are unchanged (#142).

## Version 1.31.0

### New features

- `simultaneousCIs()` gains a **multiplier-bootstrap** method
  (`method = "bootstrap"`), the Chernozhukov-Chetverikov-Kato (2013) sup-t
  bootstrap, as an alternative to the default analytic `mvtnorm::qmvnorm()`
  critical value (`method = "analytic"`, unchanged and still the default). It
  perturbs the per-unit influence functions with random multipliers
  (`multiplier = "rademacher"` or `"mammen"`, `B` replicates, optional `seed`)
  and reads the simultaneous critical value off the bootstrap distribution. The
  two methods are asymptotically equivalent; the bootstrap scales better to
  large effect families (where `qmvnorm` strains) and is
  heteroskedasticity/cluster-robust by construction. The per-unit influence
  matrix `F` reproduces the cluster-robust joint covariance (`crossprod(F) /
  (NT)^2` equals it up to the `N/(N-1)` finite-sample factor), so on a
  `se_type = "cluster"` fit the bootstrap per-effect standard errors match the
  analytic ones exactly. This
  release covers `family = "cohort"`, `"all_post_treatment"`, and `"custom"`;
  `"event_study"` (whose cohort-probability variance term needs a per-unit
  propensity influence function) and the high-dimensional `p >= NT` regime are
  planned follow-ups (#142).

## Version 1.30.0

### Improvements

- `simulateData()` / `simulateDataCore()` can now generate **high-dimensional
  (`p > NT`) panels**. Previously `genCohortTimeFE()` enforced
  `N >= (G + 1) * (d + 1)` **unconditionally** (every cohort needed at least
  `d + 1` units), which structurally confined the generator to the `p < NT`
  regime. That stricter bound is now gated behind `guarantee_rank_condition`
  (consistent with every other rank guard in the generator); with the default
  `guarantee_rank_condition = FALSE` only `>= 1` unit per cohort is required, so
  small cohorts -- and therefore `p > NT` designs -- are reachable. This lets the
  high-dimensional `debiasedATT()` path (added in 1.29.0) be exercised on
  synthetic `p > NT` data with a known data-generating coefficient vector (for
  coverage studies of the method; recovering the truth requires an adequate
  number of units). Generation with
  `N >= (G + 1) * (d + 1)` is unchanged, and `guarantee_rank_condition = TRUE`
  still enforces the full-rank bound (#293).

## Version 1.29.0

### New features

- `debiasedATT()` now also supports **high-dimensional (`p >= NT`)** fits, where
  the unrestricted ETWFE is infeasible and the fixed-`p` exact-inverse debiasing
  direction is invalid. In that regime the debiasing direction is instead the
  **nodewise (desparsified-lasso) relaxed inverse** of the design Gram (paper
  Theorem 6.6, `debiased.highdim.thm`). It is the *same* estimator and standard
  error as the fixed-`p` case — only the construction of the debiasing direction
  `v` changes — so existing `p < NT` results are unchanged (byte-identical). The
  `p >= NT` error introduced in 1.28.0 is lifted. New arguments `lambda_c`,
  `riesz_max_iter`, and `riesz_tol` control the nodewise solver (ignored when
  `p < NT`), and for `p >= NT` fits the returned list gains `feasibility`,
  `converged`, and `lambda_node` diagnostics. The high-dimensional path is
  **experimental**: its finite-sample coverage is not yet validated by
  simulation and the penalty constant `lambda_c` is not yet tuned — inspect the
  returned diagnostics and treat `lambda_c` as a tunable knob (#31).

## Version 1.28.0

### New features

- New accessor `debiasedATT()` (with a `debiasedATTWithSimulatedData()` wrapper)
  returns a *debiased* overall-ATT point estimate and a uniformly-valid standard
  error from a fitted `fetwfe()` object (`q < 1`). It is a complement to the
  fused plug-in interval (`fit$att_hat` / `fit$att_se`), which under-covers the
  aggregated overall ATT; the debiased interval restores asymptotically nominal
  coverage at roughly the efficiency of unrestricted ETWFE. The debiased point
  estimate differs from the fused `att_hat` (by the OLS identity it equals the
  ETWFE/OLS estimate in the ATT direction), so it is exposed as its own accessor
  rather than a `se_type` option. The standard error has a per-unit-clustered
  regression channel plus a cohort-weight channel that add under marginal cohort
  assignment / (Psi-IF); the assumptions are documented on `?debiasedATT`. Errors
  cleanly for `q >= 1` (no selection). (At this version it also errored for the
  high-dimensional `p >= NT` regime; that regime is supported as of 1.29.0.)
  (#291).

## Version 1.27.3

### Bug fixes

- `eventStudy()`, `cohortTimeATTs()`, and `simultaneousCIs()` no longer error on a
  fit produced with `indep_counts` when the treated cohorts adopt *later* than the
  earliest panel periods (e.g. cohorts adopting in 2005–2009 within a 2000–2010
  panel). The accessors previously failed with `first_inds[G] <= n_vars is not TRUE`
  (`eventStudy()` / `cohortTimeATTs()`) or `subscript out of bounds`
  (`simultaneousCIs()`). The cause: the `indep_counts` path stored `cohort_probs`
  unnamed, so the cohort adoption timing could not be recovered and the accessors
  fell back to the earliest-timing layout, inconsistent with the actual number of
  treatment effects. `att_hat` / `att_se` / `catt_df` are unchanged by this fix
  (it only attaches the previously-missing cohort labels) (#288).
  Older saved `indep_counts` fits also now recover their cohort timing (from
  `catt_df`) without re-fitting. Note: with `indep_counts` the `eventStudy()`
  event-time estimates and the overall `att_hat` are (as before) re-weighted by the
  independent-sample cohort probabilities; the per-cell `cohortTimeATTs()` estimates
  are unchanged by `indep_counts`.

## Version 1.27.2

### Improvements

- `genCoefs()` and `genCoefsCore()` now accept `density = 1`, which yields a
  fully dense (non-sparse) coefficient vector — every entry of the initial
  `theta` is nonzero. Previously `density` had to be *strictly* in `(0, 1)`; the
  accepted range is now `(0, 1]`. This supports simulating data in which the
  FETWFE sparsity assumption is maximally violated. `density <= 0` (which would
  zero every coefficient) and `density > 1` still error (#286).

## Version 1.27.1

### Documentation

- Surfaced the fusion-penalty structure options on high-traffic documentation
  surfaces (#282). The choice of fusion geometry --- `fusion_structure = "cohort"`
  (default) vs. `"event_study"`, and the custom `fusion_matrix` override --- is a
  core differentiator of FETWFE but was previously reachable only through a
  buried `@param`. It is now named in a new "Penalty and fusion structures"
  section of the README, cross-linked from the introductory and simulation
  vignettes to the existing "Choosing a fusion structure" vignette, pulled into
  `fetwfe()`'s `@description`, and given a `@seealso` on `fetwfe()` and
  `genCoefs()` (the two functions that expose `fusion_structure`).

## Version 1.27.0

### New features

- New accessor `cohortTimeATTs()` returns the fully disaggregated treatment
  effects --- one row per `(cohort, time)` cell, with no averaging over either
  axis --- complementing `cohortStudy()` (per-cohort, time-averaged) and
  `eventStudy()` (per-event-time, cohort-averaged). Standard errors are the
  per-cell regression standard errors, recomputed from the fit exactly as
  `eventStudy()` does (matching `simultaneousCIs(family = "all_post_treatment")`
  to floating-point precision); confidence intervals and p-values are pointwise.
  A cell zeroed out by the fusion/bridge penalty reports `estimate = 0`,
  `se = 0`, and `selected = FALSE`. Supports `fetwfe()`, `etwfe()`, and
  `betwfe()` (like `eventStudy()`); `broom::tidy()` is supported via
  `tidy.cohortTimeATTs()` (`term = "cohort_<g>_time_<t>"`) (#280).

## Version 1.26.5

### Bug fixes

- `print.twfeCovs()` now accepts a `max_event_times` argument, for parity with the
  other estimators' print methods, instead of erroring with "argument matched by
  multiple actual arguments" when one was supplied. The argument is a no-op for
  `twfeCovs` (its output has no event-study section) (#267).
- `attgtToFetwfeDf()` / `etwfeToFetwfeDf()` now validate the `out_names` argument
  (it must be a list containing the four keys `time`, `unit`, `treatment`, and
  `response`) and raise a clear error naming that contract, instead of failing later
  with a cryptic internal message — including the case of a named character vector,
  which carries the right names but is indexed with `$` downstream (#268).

## Version 1.26.4

### Improvements

- The reported model-size diagnostics `lambda.max_model_size`,
  `lambda.min_model_size`, and `lambda_star_model_size` (from `fetwfe()` /
  `betwfe()` and their `*WithSimulatedData` wrappers) now count the number of
  selected *features* and exclude the always-present intercept, matching the
  documentation ("ideally close to 0 / `p`") and `glmnet`'s `df` convention.
  Previously each also counted the intercept, so the smallest model reported 1
  instead of 0. **This shifts each of the three reported sizes down by 1.** It does
  *not* change `lambda_star`, the coefficients, the treatment-effect estimates, or
  their standard errors: the intercept is a constant offset on the BIC, so model
  selection is identical. Resolves the doc/behavior inconsistency in #269.

## Version 1.26.3

### Bug fixes

- The estimators (`fetwfe()` / `etwfe()` / `betwfe()` / `twfeCovs()`) now emit a
  warning when exactly one of `sig_eps_sq` / `sig_eps_c_sq` is supplied and the
  other is left `NA`. The two noise variances are estimated jointly by REML, so
  supplying only one silently discarded it (re-estimating both); the warning
  names the ignored value and notes that they must be supplied together or both
  omitted (#266).

## Version 1.26.2

### Bug fixes

- `print()` and `summary()` for the `FETWFE_tes` object returned by `getTes()`
  now label cohorts by calendar adoption time (cohort `g` adopts at time `g + 1`),
  matching `tidy.FETWFE_tes()` and the fitted estimators' `catt_df$cohort`.
  Previously they used the 1-based loop index, so the same cohort showed as e.g.
  "Cohort 1" under `print()` / `summary()` but "Cohort 2" under `tidy()`. The
  convention now lives in one shared helper (#261).

## Version 1.26.1

### Documentation

- `fetwfe()` and `betwfe()` runnable `@examples` now guard their
  `library(bacondecomp)` / `data(castle)` usage behind
  `if (requireNamespace("bacondecomp", quietly = TRUE)) { ... }`. `bacondecomp`
  is a Suggests-only dependency, so previously `example(fetwfe)` / `example(betwfe)`
  errored for users without it installed (#260).
- Clarified `twfeCovs()`'s title/description: it estimates a single pooled
  treatment effect per cohort (not per-(cohort, time)), matching the body text;
  added the "Optional" qualifier to `genCoefs()`'s `@param G` (#257).

## Version 1.26.0

### New features

- `twfeCovs()` now returns a fully-featured S3 object: a styled `print()` /
  `summary()` plus `tidy()` / `glance()` (broom), matching the three sibling
  estimators (`fetwfe()`, `etwfe()`, `betwfe()`). `plot()` and `augment()` are
  intentionally not defined --- `twfeCovs()` estimates one pooled effect per
  cohort (so there is no per-(cohort, time) / event-study structure to plot)
  and its coefficient vector lives in a reduced cohort-level basis (so
  `augment()`'s fitted-value path does not match); both raise an informative
  error (#58).

## Version 1.25.0

### Bug fixes

- `getTes()` no longer perturbs the caller's random-number generator. Under a
  covariate-dependent (`"multinomial"` / `"ordered"`) assignment DGP it
  integrates the expected cohort propensities by Monte Carlo, and that draw
  previously advanced the caller's `.Random.seed` as a side effect. `getTes()`
  is an analysis/truth function, so the integration now runs under a restored
  RNG (the per-cohort truth is unchanged): the caller's `.Random.seed` is saved
  and restored, leaving the user's stream untouched (#254).

### New features

- `genCoefs()`, `genCoefsCore()`, and `simulateDataCore()` now accept
  `seed = NA` to draw from the ambient random-number generator silently (no
  `set.seed()` call), matching the behavior already available on
  `simulateData()`. As before, an explicit numeric `seed` calls `set.seed()`
  (reproducible draw, global RNG left advanced) and `seed = NULL` draws from the
  ambient stream; an invalid `seed` now reports a single canonical message
  (`"seed must be NULL, NA, or a single numeric value"`) (#254).

## Version 1.24.0

### Breaking changes

- `simulateData()` now draws from the ambient random-number generator by
  default and **no longer reuses `coefs_obj$seed`** (the seed `genCoefs()` used
  to build the coefficients). Previously `simulateData()` silently called
  `set.seed(coefs_obj$seed)` internally, so repeated calls on one coefficients
  object returned byte-identical panels and a preceding `set.seed()` was
  ignored. The new `seed` argument controls the panel RNG: the default
  (`seed = NULL`) draws from the ambient generator (respecting any preceding
  `set.seed()`) and emits a warning that this default changed; pass
  `seed = <integer>` to restore reproducibility (the same integer always yields
  the same panel), or `seed = NA` to use the ambient generator silently. To
  vary the panel across Monte Carlo replications, pass a different `seed` each
  replication instead of mutating `coefs_obj$seed` (#250).

## Version 1.23.0

### Features

- `fetwfe()`, `etwfe()`, `betwfe()`, and `twfeCovs()` (and the
  `*WithSimulatedData()` wrappers), along with `genCoefs()` and `simulateData()`,
  now support a single treatment effect (one post-treatment period for the
  treated cohort), down to the true 2x2 case (`T = 2`, `G = 1`) and a
  late-adopting single cohort (`G = 1`, `T >= 3`, adopting in the final period).
  This removes the `num_treats >= 2` requirement added in #112, giving continuity
  for simulation sweeps over `T` (#251). At `T = 2` -- the unique configuration
  with nothing to fuse -- `fetwfe()` now warns that no fusion penalty applies and
  the bridge penalty reduces to individual shrinkage; the other estimators and
  the late-adopting case (which still fuses the time fixed effects) do not warn.

## Version 1.22.0

### Features

- `fetwfe()`, `etwfe()`, `betwfe()`, and `twfeCovs()` (and the
  `*WithSimulatedData()` wrappers), along with `genCoefs()` and
  `simulateData()`, now support a single treated cohort (`G = 1`): one cohort of
  units that adopt treatment plus never-treated units, rather than requiring at
  least two treated cohorts (#112). At least two treatment effects (two
  post-treatment periods for the treated cohort) are still required so the
  dynamic-effect fusion has something to fuse; a panel with a single
  post-treatment period (e.g. `T = 2`, or a cohort adopting in the final period)
  is rejected with a clear message.

## Version 1.21.1

### Minor improvements

- The estimator input validation now rejects `sig_eps_sq = 0` with a clear
  "must be positive" error. Previously only negative values were rejected, and
  a zero noise variance produced an opaque downstream failure (the GLS
  transform `1 / sqrt(sig_eps_sq)` is non-finite). Part of a small internal
  cleanup batch (#185) that also hardens several internal invariants and the
  integer return types of internal index helpers.

## Version 1.21.0

### Features

- `fetwfe()` (and `fetwfeWithSimulatedData()`) gain a `fusion_matrix` argument:
  an optional user-supplied `num_treats x num_treats` forward-differences matrix
  `D_N` for the treatment-effect block. When supplied it overrides
  `fusion_structure` for that block (the estimator uses `solve(fusion_matrix)`
  internally), letting advanced users encode arbitrary identifying-restriction
  structures beyond the two built-ins. The rows/columns follow the cohort-major
  `(g, t)` order of `getFirstInds()` / `getTreatInds()`. Under the paper's
  fixed-dimension scoping any finite, invertible `D_N` inherits the inferential
  guarantees (Assumption (D) of Faletto 2025), so `fetwfe()` validates only that
  `fusion_matrix` is finite and invertible (else it stops); a numerically
  near-singular `D_N` still fits but triggers a `warning()`. `fusion_matrix =
  NULL` (the default) is byte-identical to prior behavior (#236).

## Version 1.20.0

### Features

- `genCoefs()` and `genCoefsCore()` gain a `fusion_structure` argument
  (`"cohort"` or `"event_study"`), the simulation-side companion to the
  `fetwfe()` argument added in 1.19.0. Under `"event_study"` the generated true
  coefficients are sparse in the event-study basis (treatment effects that share
  the same time since treatment are fused across cohorts), letting simulation
  studies generate truth matched to the event-study penalty. `"cohort"` (the
  default) is byte-identical to prior behavior. The returned `FETWFE_coefs`
  object carries (and prints) a `fusion_structure` slot, and the
  "Choosing a fusion structure" vignette now generates its event-study example
  natively via this argument (#237).

## Version 1.19.1

### Documentation

- Added a vignette ("Choosing a fusion structure: cohort vs. event-study
  penalties") demonstrating the two `fetwfe()` `fusion_structure` options
  side-by-side on a simulated panel, with guidance on when to prefer each and a
  note distinguishing the estimation-time `fusion_structure = "event_study"`
  penalty from the reporting-time `family = "event_study"` aggregation (#40).

## Version 1.19.0

### Features

- `fetwfe()` gains a `fusion_structure` argument (`"cohort"` or
  `"event_study"`). Under `"event_study"` the treatment-effect fusion penalty
  fuses effects at the same time since treatment (event time `e = t - g`)
  across cohorts, instead of the default within-cohort / between-cohort two-way
  fusion --- the event-study identifying restriction common in applied DiD.
  This shrinks toward an event-study structure when the data support it,
  recovering event-time-structured effects with smaller error than the default.
  The event-study penalty satisfies the same theoretical guarantees as the
  default (Faletto 2025, event-study singular-value lemma): same singular-value
  constants, so every consistency / selection / asymptotic-normality result
  transfers verbatim. `"cohort"` (the default) is byte-identical to prior
  behavior (#40).
- Note: `fetwfe()` objects now carry a `fusion_structure` slot. Saved fits from
  prior versions lack it and will fail validation when passed to `eventStudy()`,
  `simultaneousCIs()`, or `cohortStudy()`; re-fit under 1.19.0 to restore
  compatibility.

## Version 1.18.1

### User-facing

- `print()` and `summary()` on a `getTes()` (`FETWFE_tes`) object now display a
  "Cohort assignment DGP" section -- the assignment type, strength, and
  propensity-derived cohort weights -- for objects generated under a
  covariate-dependent DGP (`assignment_type = "multinomial"` or `"ordered"`).
  Marginal-DGP objects print exactly as before. To support this, `getTes()`
  output additionally carries the `assignment_type` and `assignment_strength`
  slots, copied from the source `FETWFE_coefs` (#189).
- `fetwfe()` and `betwfe()` now warn when `lambda.max`, `lambda.min`, or
  `nlambda` are supplied under the default `lambda_selection = "cv"`, which
  builds its own cross-validation grid and ignores them. Set
  `lambda_selection = "bic"` to use a custom lambda grid (#185).

## Version 1.18.0

### Inference

- `catt_df` (cohort ATTs) and `eventStudy()` now report a multiplicity-adjusted
  `p_value` that follows the fit's `ci_type`, completing the simultaneous-
  inference default introduced in #197 (#200). Under `ci_type = "simultaneous"`
  (the default) the `p_value` is the single-step max-T (family-wise) adjusted
  p-value -- the exact dual of the simultaneous band, so a coefficient lies
  outside its `1 - alpha` band iff its adjusted `p_value < alpha`. Under
  `ci_type = "pointwise"` the `p_value` stays the per-effect Wald p-value
  (unchanged). The adjustment reuses the joint covariance assembled for the
  band (`mvtnorm::pmvnorm()`; Bonferroni under `se_type = "conservative"`), and
  `simultaneousCIs()` now returns the per-effect `adjusted_p_values`. The
  overall-ATT `att_p_value` is a single scalar (`K = 1`) and is unchanged.

## Version 1.17.5

### Performance

- Vectorized `getBetaBIC()`'s BIC-path coefficient selection (#168, #169):
  the per-`lambda` loop's separate matrix-vector products are now a single
  BLAS-3 matrix multiply, and the redundant per-`lambda` `O(N*T)` scan for
  non-finite fitted values in `sse_bridge()` was dropped, giving roughly a
  1.4x speedup of the selection step when `lambda_selection = "bic"`. The
  selected coefficients and all downstream estimates are byte-identical to
  the previous implementation.

## Version 1.17.4

### Documentation and messages

- Doc/message sweep (#208): corrected stale documentation introduced by the
  #192/#197 simultaneous-CI changes -- `cohortStudy()`'s `ci_low`/`ci_high` are
  now described by the fit's `ci_type` (simultaneous by default) rather than as
  a "Wald CI"; `simultaneousCIs()`'s `se_type` references and the `custom`-family
  variance caveat were clarified (the `Sigma_2 = 0` simplification is
  anti-conservative for a contrast that pools across cohorts in a
  probability-weighted way); and notes were added on `genCoefs()`'s `G = NULL`
  default and `simulateData()`'s deterministic re-seeding. In addition,
  `check_etwfe_core_inputs()` now emits the same friendly single-cohort error
  message as `prep_for_etwfe_core()` instead of a bare `stopifnot` failure.

## Version 1.17.3

### Bug fixes

- Fixed: removed a duplicated `ci_type` entry from the `fetwfe()` and
  `fetwfeWithSimulatedData()` `@return` documentation, so the slot is listed
  once (not twice) in `?fetwfe` / `?fetwfeWithSimulatedData`. The
  `doc-slot-parity` test now also guards against duplicate top-level `@return`
  items. (#206)

## Version 1.17.2

### Bug fixes

- Fixed: `getGramInv()` now degrades gracefully (a warning plus
  `calc_ses = FALSE`) instead of aborting with an uncaught LAPACK
  "system is computationally singular" error when the selected-feature Gram is
  near-singular (reciprocal condition number below `solve()`'s tolerance). The
  singularity guard is now condition-number-aware and `solve()` is wrapped in a
  `tryCatch()` routing to the same fallback. (#205)

## Version 1.17.1

### Bug fixes

- Fixed: `twfeCovs()` now honors its `alpha` argument when building the default
  simultaneous confidence band (previously the band was always computed at
  `alpha = 0.05` regardless of the requested level). The class gains an `alpha`
  slot. (#204)

## Version 1.17.0 (2026-06-01)

### Deprecations

- The user-facing cohort-count is now named **`G`** (the canonical name),
  matching the `did` / Callaway-Sant'Anna staggered-adoption convention. This
  applies to the `genCoefs()` / `genCoefsCore()` / `simulateDataCore()`
  cohort-count argument, the returned-list field on every estimator
  (`fetwfe()`, `etwfe()`, `betwfe()`, `twfeCovs()`) and on the simulation
  objects (`FETWFE_coefs`, `FETWFE_simulated`, `FETWFE_tes`), and the
  `print()` / `summary()` display (`Treated cohorts (G)`). The old name `R`
  is **deprecated**: passing `R =` to those functions still works but emits a
  one-time warning and is mapped to `G`, and the `$R` field is still populated
  (with the same value as `$G`). The deprecated `R` argument and `$R` field
  will be removed in a future release; pass `G` (and read `$G`) to silence the
  warning and adopt the new name (#41).

### Bug fixes

- Replaced three non-ASCII characters (`Ψ`, `′`, `→`) in the roxygen
  documentation of `eventStudy()` and `tidy.eventStudy()` with LaTeX-safe
  `\eqn{}` equivalents, so the PDF reference manual builds. These characters
  were introduced with the simultaneous-CI / tight-Gaussianity documentation
  in 1.15.0–1.16.0 and would have failed `R CMD check`'s PDF-manual build at
  the next CRAN submission (`devtools::check()` runs with `--no-manual`, so
  the per-PR checks did not surface them).

## Version 1.16.0 (2026-06-01)

### Changed defaults

- Behavior change: the default confidence-interval bounds reported for
  cohort-specific ATTs (`catt_df$ci_low` / `catt_df$ci_high`, also surfaced by
  `cohortStudy()` and `broom::tidy()`) and for event-study effects
  (`eventStudy()`, shown in `print` / `summary` / `plot` and surfaced by
  `broom::tidy()`) are now **simultaneous (family-wise, uniform) bands** rather
  than per-effect pointwise Wald intervals (#197). Each family's band now
  covers all of its effects jointly with probability `1 - alpha`, matching the
  default presentation of `did::aggte(cband = TRUE)` (which uses a multiplier
  bootstrap; `fetwfe` uses the parametric `simultaneousCIs()` machinery added
  in 1.15.0). A new fit-time argument
  `ci_type = c("simultaneous", "pointwise")` (default `"simultaneous"`)
  controls this; pass `ci_type = "pointwise"` to `fetwfe()` / `etwfe()` /
  `betwfe()` / `twfeCovs()` to restore the exact pre-1.16.0 bounds everywhere
  (the `broom::tidy()` outputs included). Only the interval bounds change:
  standard errors, per-cohort p-values, selection flags, and the overall-ATT
  confidence interval are identical under both settings.

### Dependencies

- `mvtnorm` moved from `Suggests` to `Imports`. Because simultaneous confidence
  bands are now computed by default at fit time, `mvtnorm` is no longer
  optional (a default whose output depended on whether a suggested package were
  installed would not be reproducible).

## Version 1.15.1 (2026-06-01)

### Documentation

- Expanded the inference vignette (`vignette("inference_vignette")`) with a
  new section explaining FETWFE's restriction selection consistency
  (Faletto 2025, Theorem 6.2) as an implicit test of the zero-effect null:
  a cohort the bridge penalty selects out (`catt_df$selected == FALSE`,
  `p_value == NA`) is the asymptotic conclusion that the cohort's true effect
  is zero --- a stronger statement than a confidence interval on an
  unregularized estimator can provide, where an exactly-zero estimate is a
  measure-zero event. The section attributes the conclusion to the `selected`
  flag (not a Wald statistic, which the paper notes cannot test the zero null),
  states the `q < 1` caveat, and includes a controlled-truth simulated worked
  example in which FETWFE selects out exactly the truly-zero cohort. Also
  corrected a stale theorem cross-reference in the cluster-robust section
  (the bridge oracle property is Theorem 6.3, not 6.2) (#30).

## Version 1.15.0 (2026-05-31)

### New features

- Added `simultaneousCIs()` (#192): a new exported S3 generic that returns
  parametric simultaneous `(1 - alpha)` confidence intervals over user-
  specified families of treatment effects from a fitted `fetwfe()` /
  `etwfe()` / `betwfe()` / `twfeCovs()` object. The simultaneous critical
  value is the `(1 - alpha)` quantile of `max_k |Z_k|` under the joint
  multivariate-normal distribution of the effect family, computed via
  `mvtnorm::qmvnorm(..., algorithm = mvtnorm::GenzBretz())` (mvtnorm's
  default quasi-Monte Carlo integration; deterministic in its inputs via an
  internal fixed seed that does not perturb the caller's random stream),
  under Faletto (2025) Theorem (c') tight Gaussianity and Assumption
  (Psi-IF). Supports four families: `"event_study"` (one effect per post-
  treatment event time), `"cohort"` (one effect per treated cohort),
  `"all_post_treatment"` (one effect per `(r, t)` cell), and `"custom"`
  (a user-supplied `K x num_treats` contrast matrix matching the
  `multcomp::glht()` convention). Returns a `simultaneous_cis` S3 object
  with a `$ci` data frame (estimate / simultaneous / pointwise CIs side by
  side), `$critical_value`, `$pointwise_critical_value`
  (`qnorm(1 - alpha/2)`, for reference), and `$bonferroni_critical_value`
  (`qnorm(1 - alpha/(2K))`, for reference; strictly larger than the
  simultaneous critical value when effects are positively correlated, as is
  typical in DiD). Includes `print` / `tidy` / `plot` methods following the
  existing S3 patterns; `plot()` requires `ggplot2` (already in Suggests).
  `mvtnorm` is added to `Suggests` (used only when this function is called;
  the defensive `requireNamespace()` guard fires with an actionable `stop()`
  if missing, and only when `K > 1` and `se_type != "conservative"`). A new
  vignette section in `inference_vignette.Rmd` shows a worked event-study
  example contrasting simultaneous vs per-point CIs.

### Refactor

- Consolidated three existing inline Jacobian-build sites
  (`getSecondVarTermOLS`, `getSecondVarTermDataApp`,
  `.event_study_var2_fetwfe`) into a single `.build_jacobian()` helper in
  `R/variance_machinery.R` (#192; WORKFLOW_LESSONS §14 Class A: the third
  copy triggers a refactor, ahead of `simultaneousCIs()` becoming a fourth
  call site). Byte-identity of the overall-ATT, per-cohort CATT, and per-
  event-time event-study standard errors was verified pre- and post-
  refactor.

## Version 1.14.1 (2026-05-31)

### New features

- Added optional `assignment_interactions` and
  `assignment_interaction_strength` arguments to `genCoefs()` (#191).
  When non-`NULL`, the covariate-dependent propensity DGPs introduced
  in 1.14.0 gain pairwise multiplicative interaction terms (including
  self-interactions / quadratic terms via `c(j, j)`) in the linear
  predictor — only inside the propensity model; the outcome model
  continues to use the original `X`. Pairs are canonicalized to
  `c(min, max)` and silently deduplicated. The `assignment_coefs` slot
  on `FETWFE_coefs` gains `$interactions` and `$delta` sub-slots; a
  new top-level `$assignment_interaction_strength` slot mirrors
  `$assignment_strength`. `genCoefs()` also gains a `verbose = FALSE`
  argument; setting `verbose = TRUE` emits a `message()` when
  canonicalization removes duplicate or unordered pairs (default is
  silent — users can verify the final canonical list via
  `coefs$assignment_coefs$interactions`). Unlocks adversarial-
  propensity Monte Carlo DGPs for the FETWFE vs Callaway-Sant'Anna
  comparison work tracked at #140. Default
  `assignment_interactions = NULL` preserves v1.14.0 behavior
  byte-identically on every existing test.

## Version 1.14.0 (2026-05-30)

### New features

- Added covariate-dependent cohort-assignment DGPs to `genCoefs()` and
  `simulateData()` (#162). New arguments `assignment_type` (`"marginal"`
  / `"multinomial"` / `"ordered"`) and `assignment_strength` control whether
  cohort membership in the simulated panel is drawn marginally (the current
  default behavior, preserved byte-identically) or from a multinomial-logit
  or ordered-logit (proportional-odds) propensity-score model over the
  covariates. The `FETWFE_coefs` object gains `$assignment_type`,
  `$assignment_strength`, and `$assignment_coefs` slots; the `FETWFE_tes`
  object gains a `$cohort_weights` slot. `getTes()` updated to compute the
  propensity-weighted population ATT truth under the new DGPs via Monte Carlo
  integration over the covariate distribution, matching Eq.
  `att.estimator.weighted` (Faletto 2025, line 837). New
  `simulation_vignette.Rmd` section demonstrates the three DGP modes
  side-by-side with worked examples. Unlocks realistic Monte Carlo studies
  where FETWFE has to handle conditional-on-covariate cohort selection
  (relevant to the FETWFE vs Callaway-Sant'Anna comparison work tracked at
  #140).

## Version 1.13.8 (2026-05-30)

### Internal

- Deleted an unreachable `if (length(cohorts) == 0) stop(...)` block
  in `idCohorts()` (#180 DC1). The block sits AFTER an unconditional
  `cohorts[[as.character(times[1])]] <- character()` insertion, so
  `length(cohorts) >= 1` always by the time it runs. The actual
  "all units treated in first period" catch lives earlier in the
  function at `if (length(units) == 0) stop(...)`. Zero behavior
  change; 5 lines deleted.

- Dropped a vestigial `sel_inds` list and the now-unused `first_inds`
  parameter from `getSecondVarTermOLS()` in `R/variance_machinery.R`
  (#180 DC2). The list was built solely to feed a tautological
  `stopifnot(all.equal(...))` assertion; the Jacobian construction
  below it never read the list. Both call sites
  (`getTeResultsOLS()` and the `eventStudy()` dispatcher) updated to
  stop passing `first_inds`. The sibling function
  `getSecondVarTermDataApp()` keeps its own `sel_inds` — there it
  is genuinely used to slice `d_inv_treat_sel`. Zero behavior
  change.

- Closed the cross-estimator parity gap from PR #144 (#180). FETWFE
  fits now expose `fit$calc_ses` at top level, matching `etwfe()` /
  `betwfe()` / `twfeCovs()`. Value is identical to
  `fit$internal$calc_ses` (already present). `.EXPECTED_SLOTS_FETWFE`
  and the `@return` doc updated. Pure addition; no existing slot
  removed.

### Documentation

- Reworded the misleading "canonical path is `$internal` going
  forward" comment in `R/etwfe_class.R`, `R/betwfe_class.R`,
  `R/twfeCovs_class.R` to match what the code actually does (#180).
  PR #154's `$internal` migration was paused without ever being
  completed — the print extractors in each class file still read
  top-level — so "going forward" was forward-looking but inaccurate.
  The new comment correctly identifies which sub-slots are duplicated
  at top level (the first five: `X_ints`, `y`, `X_final`, `y_final`,
  `calc_ses`) and which live only under `$internal`
  (`variance_components`, `first_year`).

- Reworded 3 inline-source "The five slots" comments at
  `R/fetwfe.R`, `R/betwfe_core.R`, `R/twfeCovs.R` to match the
  post-#179 roxygen wording (#180 carry-over from #179). Both
  reviewers of #179 flagged these as scope-appropriate for #180
  because they contained both the now-stale "five slots" wording
  and the now-retired "canonical going forward" claim.

## Version 1.13.7 (2026-05-30)

### Documentation

- Fixed three roxygen documentation issues in the four estimator
  entry points and their `*WithSimulatedData()` wrappers (#179):
  (a) eight `@return \item{se_type}{}` parenthetical notes across
  `R/fetwfe.R` / `R/betwfe_core.R` / `R/twfeCovs.R` listed only
  `"default"` and `"cluster"`, dropping `"conservative"` (added in
  v1.12.0 / PR #163); they now correctly list all three values;
  (b) six `@return $internal` describe-blocks said "The five
  sub-slots are also duplicated at top level" but the describe-list
  contained seven sub-slots; the framing now correctly notes that
  `variance_components` and `first_year` live only under
  `$internal`; (c) the internal helper `checkFetwfeInputs()` was
  missing `@param` lines for `lambda_selection`, `cv_folds`, and
  `cv_seed`. Zero behavior change; affects rendered help pages
  only.

## Version 1.13.6 (2026-05-29)

### Bug fixes

- Guarded against integer overflow in `getBetaBIC()`'s BIC formula
  (`R/fetwfe_core.R`) at panels with `N * T > .Machine$integer.max ≈
  2.15e9` (#178). Coerces `N * T` to numeric before multiplying,
  mirroring the existing `getBetaCV()` safeguard. No user-visible
  change at realistic panel sizes; the BIC opt-in path on
  >2.15e9-cell panels now computes correctly instead of failing
  with an unhelpful `stopifnot` trace.

### Internal

- Rewrote `test-gls-kronecker-block-apply-165.R` (#178) to call
  `.estimate_variance_and_gls()` directly rather than re-deriving
  both sides of the block-apply identity inline. The previous test
  asserted an algebraic identity that holds by construction; the
  new test pins the production function's `y_gls` / `X_gls` against
  an explicit-Kronecker reference, so a future refactor that
  silently breaks the GLS reshape step is caught here rather than
  via downstream end-to-end tests.

- Sharpened the `tidy.eventStudy(es, conf.level = ...)` regression
  test in `test-broom-methods.R` (#178). The previous parametric
  assertion `(td_99 width) >= (td_95 width)` held even if
  `conf.level` were completely ignored; the new sharp
  `expect_equal((qnorm(0.995) - qnorm(0.975)) * std.error)` form
  pins the exact CI half-width formula.

## Version 1.13.5 (2026-05-29)

### Bug fixes

- Fixed a silent RNG-state regression introduced in v1.13.0 (#177):
  every default-path `fetwfe()` / `betwfe()` call was overwriting
  the caller's `.GlobalEnv$.Random.seed` because `getBetaCV()` in
  `R/fetwfe_core.R` called `set.seed(cv_seed)` without saving and
  restoring the previous state. Wrapped the `set.seed` call in an
  `on.exit({ assign(".Random.seed", old, envir = .GlobalEnv) })`
  guard so the caller's RNG state is preserved across the fit. The
  v1.12.x default BIC path was unaffected (it never set a seed);
  this fix restores RNG semantics to that prior behavior on the
  CV-default path. Estimator output is unchanged.

## Version 1.13.4 (2026-05-29)

### Bug fixes

- Fixed `eventStudy()` on panels with non-consecutive cohort years.
  PR #156 (closing #138) wrapped the `eventStudy(x)` call inside
  `print.<class>()` and `summary.<class>()` in
  `tryCatch(..., error = function(e) NULL)`; on panels where treated
  cohorts adopt at scattered time indices (the canonical
  `bacondecomp::divorce` panel: states adopt no-fault divorce in
  1969-1985 with gaps at 1978, 1979, 1981, 1982, 1983) the underlying
  `eventStudy()` raised `first_inds[R] <= n_vars is not TRUE` inside
  `genInvTwoWayFusionTransformMat()` because `getFirstInds(R, T)` and
  `getNumTreats(R, T)` both assumed cohorts adopt at consecutive
  offsets `2, 3, ..., R + 1`. The tryCatch then silently omitted the
  section, so users on real-world panels saw nothing where the Event
  Study block should have rendered. Two layers of fix:
  - **Layer A (no silent swallow).** The `tryCatch` swallow in
    `R/class_helpers.R::.print_estimator_output()` and
    `.summary_estimator_output()` is removed. `eventStudy()`'s
    contract is "succeeds on any fit produced by `fetwfe()` /
    `betwfe()` / `etwfe()` / `twfeCovs()`"; if it fails on a valid
    fit, that is always a bug and the error propagates.
  - **Layer B (offset-aware D-inverse).** A new internal helper
    `getFirstIndsFromOffsets(cohort_offsets_int, T)` (in
    `R/utility.R`) parallels `getFirstInds(R, T)` but accepts an
    arbitrary strictly-increasing vector of cohort offsets. The
    helper formula `first_inds[r] = 1 + sum_{k < r} (T - offsets[k]
    + 1)` reduces to the existing `1 + (r - 1) * (2 * T - r) / 2`
    when `offsets = c(2, 3, ..., R + 1)`. `eventStudy()` derives the
    actual offsets from `names(x$cohort_probs)` (treatment-start
    years) and `x$internal$first_year` (a new internal slot exposed
    by `prepXints()`); if the cohort labels are not integer-
    coercible (synthetic genCoefs-based fixtures) the path falls
    back to `getFirstInds(R, T)`, preserving byte-identical results
    on the synthetic test suite. The validity set `V_e` at each
    event time is reworked from `seq_len(R) <= T - 1 - e` to the
    offset-aware `cohort_offsets_int <= T - e` for the same reason.
  Issue #174.

### Internal

- New `$internal$first_year` slot on `fetwfe`, `etwfe`, `betwfe`,
  `twfeCovs` results (integer or numeric scalar; the first
  (earliest) `time_var` value in the panel after `idCohorts()`
  processing). Surfaced from `prepXints()`. Used by `eventStudy()` to
  map cohort labels to panel-time-index offsets. Added to all four
  `.EXPECTED_INTERNAL_SLOTS_*` vectors and documented in the
  estimators' `@return` blocks. Pure-additive change.

### Breaking change for v1.13.3-or-earlier serialized fits

- The new `$internal$first_year` slot is added to
  `.EXPECTED_INTERNAL_SLOTS_<CLASS>` on all four estimator classes,
  so `.validate_<class>()` stops with `"Missing slot(s): first_year"`
  when `print()` or `summary()` is called on a fit serialized under
  v1.13.3 or earlier. This matches the established convention from
  #85 (the validator framework) and recurs with every prior
  `.EXPECTED_INTERNAL_SLOTS_*` addition (most recently
  `variance_components` in #141 / #146 / PR #163, v1.12.0). Users
  with saved v1.13.3 fits need to re-fit under v1.13.4.

## Version 1.13.3 (2026-05-28)

### Performance

- Vectorised `my_scale()` in `R/utility.R` (called once per `fetwfe()`
  / `etwfe()` / `betwfe()` / `twfeCovs()` call from
  `prep_for_etwfe_regression()`). Replaced the per-column
  `apply(x, 2, sd)` R-level loop with a single-pass `colSums(x_c^2)`-
  based variance, and the two `sweep(..., FUN = "-" | "/")` calls
  with broadcasted `rep(..., each = nrow(x))` subtract / divide.
  ~2.4x speedup on the function itself (10000 x 50 fixture). Output
  is numerically identical to the prior implementation up to
  floating-point reorder noise (~1e-15) on realistic full-rank
  inputs. The byte-identical-BIC regression test pinned in
  v1.13.0 (#164) passes byte-identically post-this-change,
  confirming the production code paths are unaffected. Issue #167.

## Version 1.13.2 (2026-05-28)

### Performance

- Vectorized `idCohorts()`'s per-unit row-filter loop in `R/utility.R`
  (previously `O(N^2 * T)` due to a `df[df[, unit_var] == s, ]` filter
  inside a `for (s in units)` loop). The post-fix path sorts the data
  frame once by `(unit_var, time_var)`, tabulates per-unit row counts,
  and runs the balance check, absorbing-state check, and cohort
  assignment via column-wise operations on a reshaped `T x N`
  treatment matrix. Same observable behavior: identical balance and
  absorbing-state violation messages (single grouped `stop()`, same
  lexicographically-sorted truncated listing), identical first-period-treated
  `warning()`, identical return shape. `idCohorts()` itself is ~15x
  faster at N = 2000, T = 5; end-to-end `fetwfe()` at the Phase B
  fixture speeds up an additional 1.3x on top of #165, for a combined
  **7.6x vs. pre-#165**. Issue #166.

## Version 1.13.1 (2026-05-28)

### Performance

- Replaced an explicit `kronecker(diag(N), Omega_sqrt_inv) %*% (y | X_mod)`
  GLS whitening step in `.estimate_variance_and_gls()` with the
  block-apply identity `(I_N kron A) %*% vec(M) = vec(A %*% M)`. The
  pre-fix code allocated an `(N*T) x (N*T)` Kronecker matrix on every
  call (50 MB at N = 500, 800 MB at N = 2000) before multiplying it by
  `y` and `X_mod`; the post-fix code allocates nothing larger than the
  inputs and outputs. End-to-end measured speedups on the
  Phase B fixture: **3.6× at N = 500 / 5.7× at N = 2000**. The
  identity is algebraically exact and verified bit-identical (16-digit
  parity) on randomized fixtures; new direct test in
  `tests/testthat/test-gls-kronecker-block-apply-165.R`. The speedup
  carries verbatim to `etwfe()`, `betwfe()`, and `twfeCovs()` (all four
  estimators share the GLS step). Issue #165.

## Version 1.13.0 (2026-05-28)

### Breaking change

The default `lambda_selection` for `fetwfe()` and `betwfe()` is now
`"cv"` (10-fold cross-validation on `cv.grpreg`), replacing the
previous BIC-based selection. Empirical validation across multiple
simulation regimes (#164) showed that the prior BIC selection produced
biased overall-ATT estimates at moderate sample sizes — 95% Wald CI
coverage was as low as 0.00 at N = 2000 with the paper's
second-simulation parameters, and collapsed from 0.83 to 0.33 between
N = 120 and N = 500 in the paper's first-simulation high-dim regime.
The new CV default restores near-nominal coverage in every tested
regime and additionally has the lowest MSE in the high-dim regime
that FETWFE was originally designed for.

To recover the prior behavior:

```r
fetwfe(..., lambda_selection = "bic")
betwfe(..., lambda_selection = "bic")
```

This is a numerical-results-change-for-default-callers release.
Confidence intervals and point estimates from the new default may
differ substantially from prior versions on the same data; see the
inference vignette section "Choosing the bridge penalty parameter"
for guidance.

### New features

- `lambda_selection` parameter on `fetwfe()` and `betwfe()` (and the
  `WithSimulatedData` wrappers).
- `cv_folds` and `cv_seed` parameters control the CV fitting (defaults
  `10L` and `NULL`; the latter defaults internally to `as.integer(N * T)`
  for per-dataset reproducibility without requiring the user to
  specify a seed).
- Result objects include new top-level slots `lambda_selection`,
  `cv_folds`, `cv_seed`. `glance()` exposes all three.
- `inference_vignette.Rmd` updated with a new "Choosing the bridge
  penalty parameter" subsection documenting the default change and
  the BIC opt-out.

Issue #164.

## Version 1.12.0 (2026-05-27)

### Breaking changes

- The default `att_se` value now uses the tight Gaussian variance
  `sqrt(att_var_1 + att_var_2)` from paper Theorem
  `te.asym.norm.thm`(c$'$) under Assumption (Psi-IF). (Psi-IF) is
  satisfied by the package's default cohort sample-proportions
  estimator `hat_pi_r = N_r / N`, by multinomial-logit MLE, by any GLM
  on `W | X`, and by kernel/series regression of `1{W = r}` on `X`, so
  the new default SE is asymptotically exact for the package's default
  estimator out of the box. Previous versions returned the conservative
  Cauchy-Schwarz upper bound
  `sqrt(att_var_1 + att_var_2 + 2 * sqrt(att_var_1 * att_var_2))` as
  the same-data default. The new default SE is strictly narrower; the
  corresponding CIs are strictly tighter. To recover the prior
  conservative behavior, pass `se_type = "conservative"` (intended for
  non-(Psi-IF) propensity-score estimators like
  Robins-Rotnitzky-augmented doubly-robust, which the package does not
  currently implement).

- `se_type = "cluster"` users will also see a numerical shift in
  `att_se` between 1.11.7 and 1.12.0. The cluster path's `V_1` (the
  Liang-Zeger cluster-sandwich quadratic form) and `V_2` (the
  propensity-score contribution) are now combined via the same tight
  formula `sqrt(V_1 + V_2)` by default. To recover the prior
  conservative cluster behavior, pass `se_type = "conservative"` (the
  cluster sandwich is then no longer applied — the conservative path
  uses the model-based V_1 and the Cauchy-Schwarz combination).

### New features

- New named slots on the fitted-object's `internal$variance_components`
  block expose the paper's variance decomposition by name. `V_1` is
  the Kock-2013 regression-coefficient variance contribution
  (`sigma^2 alpha' Sigma^{-1} alpha`); `V_2` is the propensity-score
  variance contribution (`v_psi(beta_0, S)`); both correspond to the
  per-unit paper-notation rescaling `V_1 = N * att_var_1`,
  `V_2 = N * att_var_2`.

- Unit-scaled variance estimators `tilde_v_N := hat_v_N / T` from paper
  line 2006 exposed alongside `V_1` / `V_2` on the variance-components
  block. `tilde_v_N` is the headline value backing `att_se`; the four
  catalogued unit-scaled variance estimators
  (`tilde_v_N_C`, `tilde_v_N_C_pi_hat`, `tilde_v_N_C_pi_hat_cons`,
  `tilde_v_N_cons`) correspond to the paper's equations
  `v.n.r.t.att.const`, `v.n.r.t.att.rand`, `v.n.r.t.att.rand.cons`,
  and `var.est.kock.wooldridge.subgauss.cons` respectively. The Wald
  CI is `[hat_T_N +- qnorm(1 - alpha/2) * sqrt(tilde_v_N / N)]` (paper
  Eq. `conf.int.form`). Issue #146.

- `inference_vignette.Rmd` reframed: same-data inference is now Gaussian
  under (Psi-IF) by default; the conservative subgaussian path is
  documented as the fallback for non-(Psi-IF) propensity-score
  estimators. The vignette also walks through the new
  paper-notation slots (`V_1`, `V_2`, `tilde_v_N`).

Issues #141 and #146.

## Version 1.11.7 (2026-05-26)

### New features

- The `covs` argument on `fetwfe()`, `etwfe()`, `betwfe()`, and
  `twfeCovs()` now accepts a one-sided formula in addition to the
  long-standing character vector. Both forms produce identical fits:

  ```r
  fetwfe(..., covs = c("x1", "x2"))   # existing
  fetwfe(..., covs = ~ x1 + x2)       # new
  ```

  Mirrors the convention used by `did::att_gt(xformla = ...)`. Only
  additive bare variable names are supported in the formula form
  (e.g., `~ x1 + x2 + x3`); interactions (`x1:x2`, `x1 * x2`) and
  `I()` transformations are rejected with a clear pointer to compute
  derived variables in the data frame first. Issue #28.

## Version 1.11.6 (2026-05-26)

### New features

- `plot.fetwfe()`, `plot.etwfe()`, `plot.betwfe()` are new S3 methods
  returning a `ggplot` object for either event-study coefficients
  (default; `type = "event_study"`, preserving the pre-PR behavior of
  the predecessor event-study-only `plot.<class>` methods) or
  per-cohort ATTs (`type = "catt"`). Mirrors the visualization style
  of `did::ggdid()` from the Callaway-Sant'Anna `did` package, giving
  users a single-call route from fitted object to publication-ready
  visualization. For the bridge-penalty estimators (`fetwfe` /
  `betwfe`) in the CATT view, points are shape- and color-coded by
  whether the bridge penalty left the cohort's ATT nonzero
  (`selected = TRUE`) or zeroed it (`selected = FALSE`); for `etwfe`
  (no selection), styling is uniform. `ggplot2` is required (Suggests
  dependency) -- the helper stops with an actionable message if
  missing. `twfeCovs` is not supported by this PR (see #58 for the
  broader treatment of `twfeCovs` class methods). Closes #29.

## Version 1.11.5 (2026-05-26)

### New features

- `print()` and `summary()` for `fetwfe()`, `etwfe()`, and `betwfe()`
  fits now include an event-study preview block alongside the existing
  ATT and CATT blocks, surfacing per-event-time treatment-effect
  estimates in the routine display. The block uses the same snake_case
  column convention as `eventStudy(result)` (`event_time`,
  `n_cohorts`, `estimate`, `se`, `ci_low`, `ci_high`, `p_value`).
  Truncation is controlled by a new `max_event_times` argument
  (defaulting to `getOption("<class>.max_event_times", 10)` for
  `print()` and 20 for `summary()`), parallel to the existing
  `max_cohorts` knob. `summary()` results also carry a new
  `event_study` field between `catt` and `model_info`. The
  event-study is computed on demand (no fit-time cache); if
  `eventStudy()` would error on the fitted object (no realistic
  configuration currently triggers this), the block is silently
  omitted rather than crashing the display. `twfeCovs` is
  unaffected -- `eventStudy()` does not support it. Issue #138.

## Version 1.11.4 (2026-05-26)

### New features

- `etwfe()`, `betwfe()`, and `twfeCovs()` results now expose an
  `$internal` list containing the design-machinery slots (`X_ints`,
  `y`, `X_final`, `y_final`, `calc_ses`), matching the structure
  `fetwfe()` already provided. Downstream consumers can now use a
  single canonical access path (`result$internal$<slot>`) across all
  four estimator classes. The same slots remain available at top level
  (e.g., `result$X_ints`) for backward compat with existing user
  scripts; no scripts are broken by this change. Issue #144.

## Version 1.11.3 (2026-05-26)

### Defensive improvements

- `tidy.eventStudy()` (`R/broom_methods.R`) now reports the missing
  column name when its input is missing a required column, parallel to
  the guard added to `tidy.cohortStudy()` in 1.11.1. Previously a
  user-mutated `eventStudy` frame produced a cryptic "arguments imply
  differing number of rows" error from `data.frame()`; now the failure
  localizes to "tidy.eventStudy(): input is missing required columns:
  ...". Required columns: `event_time`, `n_cohorts`, `estimate`, `se`,
  `p_value`. CI columns (`ci_low`, `ci_high`) are NOT required — this
  method computes its own CIs from `estimate +/- z * se`. Issue #151.

## Version 1.11.2 (2026-05-26)

### Defensive improvements

- The five cluster-sandwich quadratic-form floor sites in the variance
  machinery (PR #111, item 9; the count was six before PR #118
  consolidated `getCohortATTsFinalOLS()` into `getCohortATTsFinal()`)
  now layer a two-tier diagnostic on top of the existing `max(., 0)`
  floor: a `warning()` fires when the quadratic form is more negative
  than `-1e-10` (clearly outside floating-point noise) and a `stop()`
  fires when it is more negative than `-1` (catastrophic; outside any
  plausible SE-squared magnitude in difference-in-differences
  applications). On well-conditioned data the behavior is unchanged.
  The diagnostic surfaces latent bugs that would previously have
  silently produced `SE = 0`. Issue #139.

## Version 1.11.1 (2026-05-25)

### New features

- `cohortStudy(result)` is a new function-style accessor for per-cohort
  average treatment effects, parallel to `eventStudy()`. It returns the
  same information available via `result$catt_df` in a discoverable
  function with its own help page (`?cohortStudy`). Use
  `tidy(cohortStudy(result))` (with `broom` loaded) for broom-shape
  output suitable for `modelsummary` and similar consumers. The returned
  object carries class `c("cohortStudy", "catt_df", "data.frame")`: the
  `cohortStudy` class dispatches `tidy()`; the `catt_df` class preserves
  the helpful-error layer for pre-1.11.0 column-name access; the
  `data.frame` base preserves the standard data-frame methods.

## Version 1.11.0 (2026-05-25)

### Breaking changes

- The columns of `result$catt_df` have been renamed to a snake_case
  convention that matches the columns of `eventStudy(result)`:

  | Old name       | New name   |
  |----------------|------------|
  | `Cohort`       | `cohort`   |
  | `Estimated TE` | `estimate` |
  | `SE`           | `se`       |
  | `ConfIntLow`   | `ci_low`   |
  | `ConfIntHigh`  | `ci_high`  |
  | `P_value`      | `p_value`  |

  The `selected` column (`fetwfe` / `betwfe` only) is unchanged. Broom
  output via `tidy(result)` is also unchanged (still broom-style:
  `estimate`, `std.error`, `statistic`, `p.value`, `conf.low`,
  `conf.high`).

  Update existing code that references `result$catt_df[["Estimated TE"]]`
  etc. The `catt_df` now carries S3 class `c("catt_df", "data.frame")`
  and provides helpful-error methods on `[[`, `$`, `[` that intercept
  the old column names and point to the new name --- so most migration
  friction surfaces as an actionable `stop()` rather than a silent
  `NULL`. Standard data-frame operations (`print()`, `head()`,
  `summary()`, `dim()`, row indexing) all dispatch through to
  `data.frame` methods unchanged. NSE access (dplyr / ggplot / subset)
  strips the `catt_df` class as usual and surfaces the standard
  "object not found" error.

  The cosmetic side of the rename: column headers in the `print()` /
  `summary()` CATT preview blocks now appear in snake_case. The
  snapshot tests in `tests/testthat/_snaps/print-method-snapshot.md`
  have been re-recorded.

## Version 1.10.0 (2026-05-21)

- CRAN release. Version 1.10.0 is the first CRAN release since 1.5.0; it
  consolidates the 1.5.1-1.9.33 development series documented in the
  entries below and introduces no new code beyond 1.9.33. The changes
  most likely to matter when upgrading from the 1.5.0 CRAN release:
  - **Behavior change.** Version 1.9.1 corrected the internal
    variance-component estimator, which had collapsed the unit-level
    variance component to nearly zero and degenerated the GLS weighting
    step. Both standard errors and point estimates shift after upgrading;
    see the 1.9.1 entry for details.
  - New functionality: `eventStudy()` with `plot()` methods (1.7.0);
    broom `tidy()` / `glance()` / `augment()` methods (1.9.0); an
    experimental cluster-robust `se_type = "cluster"` option (1.6.0); the
    `allow_no_never_treated` argument (1.5.6); per-cohort `P_value` and
    `selected` columns (1.5.5).
  - Additional standard-error and correctness fixes (1.8.0, 1.9.2-1.9.13),
    extensive internal hardening, and documentation and vignette updates.

## Version 1.9.33 (2026-05-20)

- Replaced the remaining `bacondecomp::divorce` examples (the four
  estimators' help-page examples and the package's quick-start example)
  with `bacondecomp::castle` examples, matching the main vignette; fixed
  a stale citation year in the package documentation.

## Version 1.9.32 (2026-05-20)

- Replaced the main vignette's flagship empirical example: the
  `bacondecomp::divorce` analysis (on which FETWFE selects a null under
  the corrected variance estimator, on a rank-deficient panel) is
  replaced by a `bacondecomp::castle` analysis of castle-doctrine laws
  and homicide, which yields a non-null, literature-consistent result.
  The vignette's simulated example now also generates and uses
  time-invariant covariates. Vignette content only; no change to package
  behavior.

## Version 1.9.31 (2026-05-20)

- Corrected the documented output column names for `attgtToFetwfeDf()`
  and `etwfeToFetwfeDf()`: their `@return` and `@param` blocks listed
  `time` / `unit` / `y`, but the converters produce `time_var` /
  `unit_var` / `response` (the `out_names` defaults). Documentation fix
  only; no change to either function's behavior.

## Version 1.9.30 (2026-05-20)

- Documentation and cleanup sweep from the 2026-05-19 internal review:
  documented the `cohort_times` slot in the `getTes()` return value;
  documented the `se_type` argument in the four internal `*_core`
  functions; extended the documentation-slot-parity test to cover
  `getTes()`, `genCoefs()`, and `simulateData()`; corrected a misleading
  factor-encoding comment in `processFactors()`; removed a redundant
  cohort-count assertion in `prep_for_etwfe_core()`; and fixed two stale function
  names in `cran-comments.md`. Documentation and cleanup only; no
  user-facing behavior change.

## Version 1.9.29 (2026-05-20)

- Added regression-test coverage for four code paths flagged by the
  2026-05-19 internal review: `processCovs` first-period covariate
  broadcasting when unit IDs are not already in sorted order,
  cluster-robust standard errors on auto-truncated (all-treated) panels
  for `etwfe()` and `twfeCovs()`, `augment()` row-order invariance after
  first-period-treated units are auto-dropped, and the
  `getCohortATTsFinal()` `psi_mat` row-count contract across the OLS and
  bridge calling conventions. Test-only; no user-facing behavior change.

## Version 1.9.28 (2026-05-19)

- Consolidated three small shared helpers across the estimator cores into
  `R/utility.R` (GitHub #119, follow-up to #118): the 4-way
  `grpreg::gBridge` dispatch + lambda-path diagnostic block (formerly
  duplicated across `R/fetwfe_core.R` and `R/betwfe_core.R`), the
  treat-index / treat-interaction-index computation block (formerly
  across all three `*_core.R` files), and the cohort-count rank-condition
  check (formerly across `R/fetwfe.R` and `R/twfeCovs.R`). Pure internal
  refactor with one minor user-visible change: `twfeCovs()`'s warning /
  stop message under the cohort-count rank condition now uses the
  stronger "the design matrix is rank-deficient" wording that `etwfe()`
  always used, converging on the mathematically accurate form.

## Version 1.9.27 (2026-05-19)

- Unified the internal `getCohortATTsFinal()` and previously-separate
  `getCohortATTsFinalOLS()` helpers in `R/variance_machinery.R` into a
  single function controlled by the existing `fused` flag plus a new
  required `include_selected` argument (GitHub #118). Pure internal
  refactor: cohort-ATT point estimates, SEs, p-values, and confidence
  intervals are byte-identical for `fetwfe()`, `etwfe()`, `betwfe()`,
  and `twfeCovs()` on every existing test fixture. Reduces ~140 LOC of
  duplication and consolidates the surface area for future
  variance-estimator fixes.

## Version 1.9.26 (2026-05-19)

- Fixed four roxygen-comment typos surfaced by the 2026-05-19 periodic
  code review (GitHub #115): `#"` instead of `#'` on the
  `simulateData()` and `simulateDataCore()` docstrings in
  `R/gen_data.R`, and an unmatched backtick before `fetwfe()` in the
  `attgtToFetwfeDf()` and `etwfeToFetwfeDf()` descriptions in
  `R/convert_dfs.R`. The latter fix restores correct inline-code
  rendering in `man/attgtToFetwfeDf.Rd` and `man/etwfeToFetwfeDf.Rd`;
  the former is a latent typo that roxygen2 happened to parse around.

## Version 1.9.25 (2026-05-19)

- Reorganized the `R/` source tree so each file's contents match its
  name (GitHub #82). Pure file-move refactor: zero behavior change,
  NAMESPACE unchanged, all 1909 tests pass. Specifically:
  - Renamed `R/ols_calcs.R` to `R/variance_machinery.R` (via `git
    mv`) and consolidated the FETWFE-side variance helpers
    (`getTeResults2`, `getPsiRFused`, `getSecondVarTermDataApp`,
    `getCohortATTsFinal`) that previously lived in
    `R/fetwfe_core.R` so all OLS- and FETWFE-side variance
    machinery now shares one home.
  - Created `R/design_matrix.R` for design-matrix construction
    (`prepXints`, `processCovs`, `processFactors`, `addDummies`,
    `generateFEInts`, `genTreatInts`, `genXintsData`,
    `genTreatVarsRealData`) - moved from `R/etwfe_core.R`.
  - Created `R/fusion_transforms.R` for the five fusion-transform
    matrix generators (`transformXintImproved`,
    `untransformCoefImproved`, `genBackwardsFusionTransformMat`,
    `genBackwardsInvFusionTransformMat`,
    `genInvTwoWayFusionTransformMat`, and
    `genFullInvFusionTransformMat`) - moved from `R/fetwfe_core.R`
    and `R/core_funcs.R`.
  - Split `R/gen_funcs.R` into three smaller files:
    `R/gen_coefs.R` (coefficient generators + truth extraction:
    `genCoefs`, `genCoefsCore`, `getTes`, `getActualCohortTes`),
    `R/gen_data.R` (data simulation pipeline: `simulateData`,
    `simulateDataCore`), and `R/sim_helpers.R` (internal sim
    helpers).
  - Moved `etwfeWithSimulatedData()` from `R/etwfe_core.R` to
    `R/fetwfe.R`, alongside `etwfe()` itself (matching the pattern
    of `betwfeWithSimulatedData` co-located with `betwfe_core`, etc.).
  - Moved `prep_for_etwfe_core()` from `R/etwfe_core.R` to
    `R/core_funcs.R` (alongside its sibling
    `prep_for_etwfe_regression()`).
  - Moved the shared input-validator helpers
    `.collect_etwfe_input_violations()` and
    `.format_input_violations()` from `R/etwfe_core.R` to
    `R/utility.R` (co-located with `.truncate_violations()`).
  - `R/etwfe_core.R` and `R/fetwfe_core.R` now contain only their
    namesake `_core` estimator plus the corresponding input
    validator (and, for fetwfe, `getBetaBIC`).

## Version 1.9.24 (2026-05-19)

- Reconciled three strict-vs-lenient validator inconsistencies (GitHub
  #109, follow-up to #84 item 12):
  - **`lambda.max` (Case 1)**: Tightened the inner `fetwfe_core` check
    from `lambda.max >= 0` to `lambda.max > 0` to match the outer
    `checkFetwfeInputs` invariant. No user-visible behavior change —
    the inner check is downstream of the outer and was effectively
    dead code (any caller already validated `> 0` via the entry
    point).
  - **`sig_eps_c_sq = 0` (Case 2)**: Loosened `simulateData()` (and
    its inner helper `testGenRandomDataInputs()`) to accept
    `sig_eps_c_sq = 0`, matching the validator surface used by
    `fetwfe()` / `etwfe()` / `betwfe()` / `twfeCovs()`.
    Methodologically, `sig_eps_c_sq = 0` yields a panel with no
    unit-level random effects (`rnorm(N, sd = 0) = rep(0, N)`).
    Backward-compatible additive change — previously-rejected inputs
    are now accepted.
  - **`R >= 2` (Case 3)**: Tightened the upstream `R >= 1` checks at
    `R/core_funcs.R::check_etwfe_core_inputs` and
    `genFullInvFusionTransformMat` to `R >= 2`, matching the
    user-facing helpful-error check at `R/etwfe_core.R:1127`
    ("Only one treated cohort detected..."). No user-visible behavior
    change — `R = 1` was already rejected with a clearer message
    downstream. Single-treated-cohort support (R = 1) is filed
    separately as #112.

## Version 1.9.23 (2026-05-19)

- Closes out the remaining items of GitHub #84 (the 2026-05-17 periodic
  code review) in one bundled PR: 11 items across three categories (test
  coverage, defensive guards, and performance). Internal changes only;
  no public-API drift.

  - **Test coverage** (items 1, 3, 5, 6, 7): closes 5 untested-path gaps.
    `eventStudy()` now has direct coverage on auto-truncated panels;
    the q >= 1 ridge-regime `att_se = NA` contract has tightened
    assertions covering `att_p_value`, `catt_ses`, and the CI columns of
    `catt_df`; `augment.<class>()` is now exercised with both
    `se_type = "cluster"` and `indep_counts`-bearing fits; and the
    end-to-end REML + cluster-SE pipeline is verified non-conflicting on
    a single fit.

  - **Defensive guards** (items 8, 9, 10): defense-in-depth additions to
    catch future refactor regressions. (1) `.select_att_branch()`'s
    etwfe/twfeCovs unconditional SE-non-NA guard is now gated on
    `calc_ses` (today only `q >= 1` produces `calc_ses = FALSE`, but the
    gate hardens against a future code path that routes a `calc_ses =
    FALSE` non-bridge fit through the indep-counts branch). (2) The 6
    cluster-sandwich quadratic-form computation sites now floor at zero
    (`max(quad_form, 0)`); the Liang-Zeger sandwich is PSD in exact
    arithmetic, but the `rowsum`-based meat aggregation could in
    principle round a near-zero form to a tiny negative, which would
    NaN downstream `sqrt()`s. (3) `.compute_cluster_robust_sandwich()`
    now `tryCatch`-es the `solve(crossprod(X_S_centered))` call to
    surface the same "Gram matrix corresponding to selected features is
    not invertible" message as `getGramInv()` rather than a raw Lapack
    "system is exactly singular" trace.

  - **Performance** (items 15, 16, 18): (1) `processCovs()` (the
    covariate-cleanup helper called by every estimator entry point) is
    vectorized. The old per-unit nested-loop pattern was O(N^2 * T * d);
    the rewrite uses a single first-period slice + `rep(..., each = T)`
    broadcast, dropping the cost to O(N * T * d). Output is bit-identical
    on the standard test fixtures. (2) `.estimate_variance_and_gls()`
    now computes `Omega^(-1/2)` via a closed-form spectral decomposition
    (`Omega = sig_eps_sq * I_T + sig_eps_c_sq * J_T` has eigenvalue
    `sig_eps_sq + T * sig_eps_c_sq` on `span(1_T)` and `sig_eps_sq` on
    its orthogonal complement) instead of `expm::sqrtm(solve(Omega))`.
    The closed form is exact up to floating-point and eliminates the
    runtime `expm` dependency; `expm` moves from `Imports:` to
    `Suggests:` (it's still used for the numerical-equivalence test).
    (3) `.compute_cluster_robust_sandwich()`'s N-loop is replaced with a
    vectorized `rowsum()` + `crossprod()` assembly, equivalent to the
    pre-rewrite loop bit-for-bit modulo floating-point summation order.

  Items 2 (test coverage already provided), 12 (deferred as GitHub
  #109), and 17 (deferred as GitHub #110, pending benchmark) were
  scoped out of this PR. Items 4, 11, 13, 14 shipped previously (PRs
  #108, #105, #106, #107).

## Version 1.9.22 (2026-05-19)

- `tidy.FETWFE_tes()` now accepts `conf.int` and `conf.level` arguments
  for broom-convention parity with the sibling `tidy.fetwfe()`,
  `tidy.etwfe()`, `tidy.betwfe()`, and `tidy.eventStudy()`
  methods. Defaults are `conf.int = TRUE` (preserves pre-fix output —
  NA-valued `conf.low` / `conf.high` columns included) and
  `conf.level = 0.95`. When `conf.int = FALSE`, the CI columns are
  omitted entirely. `conf.level` is validated but unused (there is no
  sampling distribution for a population truth). GitHub #84 item 14.

## Version 1.9.21 (2026-05-19)

- Added compact `print` methods for `FETWFE_simulated` (output of
  `simulateData()`) and `FETWFE_coefs` (output of `genCoefs()`).
  Pre-fix, `print(simulateData(...))` fell through to `print.list`
  and dumped the full `N*T x p` design matrix — potentially hundreds
  of MB for medium-sized panels. The new methods print a compact
  dimensions summary (N, T, R, d, p, cohort sizes, noise variances
  for `FETWFE_simulated`; R, T, d, beta length, theta sparsity, seed
  for `FETWFE_coefs`). Live behavior is unchanged for everything
  except the printed output. GitHub #84 item 13.

## Version 1.9.20 (2026-05-19)

- User-facing input-validation error messages for the four estimator
  entry points (`fetwfe()`, `etwfe()`, `betwfe()`, `twfeCovs()`) are
  now actionable instead of cryptic (GitHub #84, item 11). The legacy
  `Error: is.character(time_var) is not TRUE` family produced by
  `stopifnot()` is replaced with a named, multi-line `stop()` message
  of the form
  `Invalid inputs:\n  - <arg>: <expected>; got <received>`, where
  ALL malformed args are collected on a single call so the user sees
  every violation at once instead of one error round-trip per arg.
  Two new internal helpers in `R/etwfe_core.R` ---
  `.collect_etwfe_input_violations()` and `.format_input_violations()`
  --- absorb the collect-and-render plumbing; `checkEtwfeInputs()` and
  `checkFetwfeInputs()` (the latter at `R/fetwfe_core.R`) became thin
  wrappers around them. Per-arg failures cascade within an arg (the
  first failure in an arg's check chain short-circuits — e.g.,
  `time_var %in% colnames(pdata)` isn't run if `time_var` isn't a
  length-1 character) but collect independently across args. NULL
  inputs are handled cleanly (`class(NULL) = "NULL"`, `length(NULL) = 0`
  interpolate fine). The PR #103 snapshot guardrail (which locked the
  legacy `is.foo() is not TRUE` strings) was deliberately re-blessed
  with the new helpful messages; the snapshot mechanism continues to
  lock the new contract against accidental drift. Internal validators
  (`check_etwfe_core_inputs()`, `prepXints()`) were intentionally left
  as terse `stopifnot()`s — they check post-prep invariants, not user
  input. Item 12 of #84 (strict-vs-lenient reconciliations:
  `lambda.max >= 0` vs `> 0`; `sig_eps_c_sq = 0`; `R >= 1` vs `R >= 2`)
  is deferred to a follow-up since it involves API-behavior decisions
  rather than UX wording.

## Version 1.9.19 (2026-05-19)

- Internal consolidation of the input/output pipeline shared by the
  four public estimator entry points `fetwfe()`, `etwfe()`,
  `betwfe()`, and `twfeCovs()` (GitHub #79). Two new internal helpers
  in `R/utility.R` --- `.run_estimator_input_prep()` and
  `.select_att_branch()` --- absorb the verbatim Steps 3-5 (input
  validation + auto-truncation + design-matrix prep) and Step 8
  ("pick from in-sample vs indep ATT/SE/cohort-probs") that were
  previously copy-pasted across the four entry points. The four
  caller-side rewrites trim ~250 LOC of duplication while leaving
  the eleven other steps (match.arg, original-covs capture, core
  call, rank-deficiency warning loop, p-value, output-list build,
  validator, class assignment) inline. No public API change; the
  snapshot guardrail from v1.9.18 (PR A, #103) confirms all four
  validators' error messages remain byte-identical. The new helpers
  carry `@noRd` / `@keywords internal` and do not appear in the
  exported docs surface.

## Version 1.9.18 (2026-05-19)

- Internal refactor of `R/core_funcs.R::prep_for_etwfe_regression()`
  (GitHub #81). The 237-line god function is split into four named
  helpers — `.estimate_variance_and_gls()`, `.collapse_design_for_twfe_covs()`,
  `.append_ridge_rows()`, `.compute_cohort_probs()` — with the
  orchestrator reduced to an ~100-line sequence of helper calls. No
  public API change; all four caller sites (`fetwfe_core()`,
  `etwfe_core()`, `betwfe_core()`, `twfeCovs_core()`) continue to
  unpack the same 12-field result list. Four direct return-shape
  tests added in `tests/testthat/test-prep-helpers-81.R` to lock the
  helpers' contracts against future drift.

## Version 1.9.17 (2026-05-19)

- Small consistency-and-cleanup bundle (GitHub #76). Eight items:
  (1) `twfeCovs()` now returns an object with class `"twfeCovs"`, with
  minimal-stub `print.twfeCovs()` and `coef.twfeCovs()` methods to bring
  it into parity with the three sibling estimators (`fetwfe`, `etwfe`,
  `betwfe`). A new `R/twfeCovs_class.R` file holds these alongside a
  NULL doc stub. The class dispatcher in `.assert_estimator_object()`
  was extended to recognize `twfeCovs`. Live `print()` and `coef()`
  behavior is preserved exactly.
  (2) Three `augment.*` docstrings (for `fetwfe`, `etwfe`, `betwfe`)
  were rewritten to reflect the actual auto-trim + auto-sort behavior
  of `.augment_estimator_output()`. The pre-fix docstrings claimed the
  user must pre-trim and pre-sort; that has been wrong since #27.
  (3a) `checkEtwfeInputs()` `@title` now correctly references the
  `etwfe` / `betwfe` / `twfeCovs` functions (was: `fetwfe`).
  (6) Deleted the unused, internal `genInvFusionTransformMat()` helper.
  Confirmed zero callers before removal.
  (7a) Removed Test 4 from `tests/testthat/test-fetwfe-var2-fix.R`:
  its "hand-computed" reference reproduced the helper's own loop
  structure step-for-step, so passing only confirmed the helper was
  self-consistent, not that it matched an independent target. Tests 2
  (Monte Carlo validation) and 3 (anti-regression) are the
  load-bearing var2 coverage.
  (7b) Demoted the tautological `.fitted + .resid == y` round-trip
  expectation in `tests/testthat/test-augment-treatment-parity.R` to
  a comment. The load-bearing row-order check is preserved.
  (8a) Dropped unused `p` argument from `getCohortATTsFinal()`.
  (8b) Dropped unused `psi_mat` and `tes` arguments from
  `getSecondVarTermDataApp()`, with corresponding test call-site
  updates.
  (8c) Dropped the dead `fused = FALSE` branch (and the `fused`
  argument) from `getSecondVarTermDataApp()`, along with the
  corresponding roxygen `@param fused` block and stale docstring
  sentence about the fused-only restriction.
- The issue's Item 5 (NULL-alpha defensive fix in `tidy.*`) was
  surveyed and skipped: the post-#85 constructor validators catch
  NULL `alpha` before any `tidy.*()` method runs, so the proposed
  defensive default would be unreachable. See `.plans/cleanup-bundle-76/PLAN.md`
  Decision D7 for the full trace.

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

- Fixed a cross-method consistency bug in `eventStudy()` (GitHub #73):
  the function previously reported finite SEs and p-values for `q >= 1`
  fits even though the corresponding `fetwfe()` / `betwfe()` objects
  correctly return `att_se = NA` in that regime (the bridge oracle
  property is required for the model-based SEs). After this fix,
  `eventStudy()` propagates the fit's `calc_ses` status: when the fit
  has no valid SEs, `eventStudy()` returns NA SEs and NA p-values,
  matching the parent object's contract. The point estimates remain
  valid (and finite) in both regimes. Same fix applies under
  `se_type = "cluster"`.
- Added internal method-entry preconditions for cross-class consumers
  (`eventStudy`, `augment`, `tidy`, `glance`, `plot`, `coef`). Each
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
  outputs of `eventStudy()` and `getTes()`, so users can pipe estimator
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
  to the estimator or a pre-trimmed panel. `tidy.eventStudy()`
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
  event-study output (`eventStudy()` and `plot.fetwfe()` columns
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

- Added an exported `eventStudy(x, alpha)` function and S3 methods
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
  `eventStudy()` machinery and are useful for downstream tooling.

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
