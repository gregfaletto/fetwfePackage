# Debiased overall-ATT with a uniformly-valid standard error

Computes a *debiased* estimate of the overall average treatment effect
on the treated (ATT) from a fitted
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
object, together with a uniformly-valid standard error. This is a
complement to the fused plug-in interval reported on the fit
(`fit$att_hat` / `fit$att_se`): the fused SE is the within-selection
variance and tends to *under-cover* the aggregated overall ATT (it omits
the between-selection term), whereas `debiasedATT()` returns an interval
with asymptotically nominal coverage at roughly the efficiency of
unrestricted ETWFE.

The debiased point estimate **differs** from the fused `fit$att_hat`:
when `p < NT`, by the OLS identity (paper eq. `debiased.ols.identity`)
it equals the unrestricted ETWFE/OLS estimate in the ATT direction. You
therefore get a dual offering — `fit$att_hat` / `fit$att_se` (fused,
efficient, pointwise) and `debiasedATT(fit)` (debiased, ETWFE-efficient,
uniformly valid). It is exposed as a separate accessor (rather than a
`se_type` option) precisely because the debiased SE accompanies a
different point estimate than the fused one.

When `p >= NT` the unrestricted ETWFE is infeasible, so the OLS identity
no longer applies; the accessor instead builds the debiasing direction
by a nodewise (desparsified-lasso) relaxed inverse (**Theorem
`debiased.highdim.thm`** of Faletto (2025)). The point estimate, SE
formula, and return value are otherwise identical — "one estimator, two
regimes." The high-dimensional overall-ATT interval is **uniformly
valid** under that theorem, and its coverage is validated near-nominally
in simulation (0.94 at the `p >= NT` anchor of Faletto (2025), robust
across the nodewise penalty scales studied); see Assumptions.

## Usage

``` r
debiasedATT(
  fit,
  alpha = NULL,
  method = c("analytic", "bootstrap"),
  B = 1000L,
  seed = NULL,
  multiplier = c("webb", "rademacher", "mammen"),
  lambda_c = 1,
  riesz_max_iter = 5000L,
  riesz_tol = 1e-09,
  cv_time_budget = Inf
)
```

## Arguments

- fit:

  A fitted object from
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).
  In the fixed-p (`p < NT`) regime the SE needs either bridge selection
  (`q < 1`, `gls = TRUE`) or an un-whitened fit (`gls = FALSE`, any `q`,
  the Omega-free cluster-robust path; \#312); a GLS-whitened `q >= 1`
  fit is rejected. The high-dimensional (`p >= NT`) path accepts any
  fit.
  [`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
  /
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
  /
  [`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
  are not supported (the construction is specific to the FETWFE
  transformed-coefficient space).

- alpha:

  Numeric in `(0, 1)`; the confidence level is `1 - alpha`. Defaults to
  the `alpha` stored on the fit.

- method:

  How to form the confidence interval. `"analytic"` (default) uses the
  two-channel unit-clustered sandwich SE with a Gaussian critical value.
  `"bootstrap"` replaces the Gaussian critical value with a studentized
  score / influence-function **wild cluster bootstrap** (#360),
  **floored at the Gaussian quantile** so the interval is never narrower
  than the analytic Wald interval. The point estimate and the reported
  `se` are unchanged; only the critical value changes.

  **This is not a few-clusters remedy.** As an *unrestricted, no-refit*
  score bootstrap it does not reproduce the tail inflation of the
  *restricted* wild-cluster bootstrap-t, and under heterogeneous cluster
  influence its critical value (before the floor) falls *below* the
  Gaussian — so under the floor it reduces to the analytic interval
  exactly in the small-`N` regime; it only ever *widens* the interval
  when cluster influence is near-homogeneous. For a genuine few-clusters
  correction the restricted bootstrap-t or a CR2-type analytic
  adjustment is required (tracked as future work, \#361). Not supported
  for `indep_counts` (two-sample) fits.

- B:

  Integer; the number of wild-bootstrap replicates (default `1000`).
  Ignored unless `method = "bootstrap"`.

- seed:

  `NULL` (draw from the ambient RNG, the default) or a single integer
  for a reproducible bootstrap (the RNG state is saved and restored).
  Ignored unless `method = "bootstrap"`.

- multiplier:

  The wild-bootstrap weight distribution: `"webb"` (default),
  `"rademacher"`, or `"mammen"`. `"rademacher"` (`±1`) gives a constant
  studentization denominator (a *percentile* bootstrap); `"webb"` and
  `"mammen"` vary it (a *studentized* statistic). Because the critical
  value is floored at the Gaussian quantile, the multiplier only affects
  the near-homogeneous case where the bootstrap widens above the floor;
  in the small-`N` / heterogeneous regime the interval is the analytic
  one regardless of the weight. Ignored unless `method = "bootstrap"`.

- lambda_c:

  The leading constant of the high-dimensional (`p >= NT`) nodewise
  penalty `lambda_node = lambda_c * max(|a|) * sqrt(log(p) / N)` (`N` =
  number of units). Either a single positive number (a fixed constant;
  default `1.0` = the theory scale) or the string `"cv"`, which selects
  the constant **per fit** by cross-validating the unit-level Riesz loss
  `0.5 v'Sigma v - a'v` over the KKT-feasible region (the
  desparsified-lasso / auto-DML standard; van de Geer;
  Chernozhukov-Newey-Singh), falling back to the theory scale `1.0` when
  no grid penalty is feasible (#295). Larger values shrink the debiasing
  direction more. **Ignored when `p < NT`.** The default stays the fixed
  `1.0` (the theory scale); the high-dimensional coverage validated in
  Faletto (2025) is robust across the penalty values studied — a
  cross-validated selector below the theory scale and a fixed constant
  above it, which bracket the default `1.0` — so it needs no tuning for
  valid inference.

- riesz_max_iter, riesz_tol:

  Integer / numeric; coordinate-descent controls for the
  high-dimensional nodewise solver. **Ignored when `p < NT`.**

- cv_time_budget:

  Numeric; a wall-clock backstop (in seconds) for the `lambda_c = "cv"`
  cross-validation. `Inf` (the default) leaves selection fully
  deterministic / reproducible; a finite value stops the CV early and
  falls back to the theory scale (`lambda_c = 1.0`) with a warning, so
  an adversarial high-dimensional draw cannot spin indefinitely (#384).
  Whether a finite budget fires depends on machine speed, so the result
  can too. Ignored unless `lambda_c = "cv"`.

## Value

An object of S3 class `"debiased_att"` (a named list, with a `print` and
a `tidy` method) with elements:

- att:

  Numeric scalar; the debiased overall-ATT point estimate.

- se:

  Numeric scalar; the uniformly-valid standard error,
  `sqrt(var_reg + var_weight)`.

- ci_low, ci_high:

  Numeric; the `1 - alpha` interval. With `method = "analytic"`
  (default) this is the Wald interval `att +/- qnorm(1 - alpha/2) * se`;
  with `method = "bootstrap"` the critical value is the bootstrap
  `crit_value` (floored at `qnorm(1 - alpha/2)`, so never narrower than
  the Wald interval) in place of `qnorm(1 - alpha/2)`.

- var_reg:

  Numeric; the regression (outcome) channel's contribution to `se^2`,
  per-unit-clustered.

- var_weight:

  Numeric; the cohort-weight channel's contribution to `se^2` — the
  fit's plug-in propensity variance (`att_var_2`), which carries the
  correct single- or two-sample (`indep_counts`) formula.

`var_reg` and `var_weight` are the two **`se^2`-scale** (unit-scale)
components, so `se = sqrt(var_reg + var_weight)`. They are *not* the
paper's CLT-scale variances `V_1^full` / `V_2` from Theorem
`debiased.att.thm`, which are larger by a factor on the order of `NT`.

For high-dimensional (`p >= NT`) fits the list additionally contains
`feasibility` (`= ||Sigma v - a||_inf`, the nodewise KKT certificate,
which should be `<= lambda_node`), `converged` (the coordinate-descent
flag), `lambda_node` (the penalty used), `lambda_c` (the leading
constant used), and `lambda_c_selection` (`"fixed"` or `"cv"`). When
`lambda_c = "cv"` it also contains `lambda_cv`, the cross-validation
diagnostics (the grid, the per-grid feasibility flags and CV losses, and
whether the theory-scale fallback fired). These are absent for `p < NT`
fits.

With `method = "bootstrap"` the list additionally contains `method`,
`crit_value` (the studentized wild-bootstrap critical value, floored at
`qnorm(1 - alpha/2)`), `B`, `multiplier`, and `alpha`. The analytic
default adds none of these, so its
[`names()`](https://rdrr.io/r/base/names.html) are unchanged.

## Details

The standard error has two channels that **add** under marginal cohort
assignment / Assumption (Psi-IF): a regression (outcome) channel,
per-unit-clustered, and a cohort-weight channel for the sampling noise
in the estimated cohort weights. See *Assumptions* below.

The cohort-weight channel is the fit's own plug-in propensity variance
(`att_var_2`), the same quantity `fit$att_se` uses — so it automatically
carries the correct single-sample formula and, for `indep_counts` fits,
the two-sample formula (the weight channel concerns the estimated cohort
weights, not the selection, so it is shared by the fused and debiased
estimators). The fit's `se_type` does **not** affect the regression
channel — `debiasedATT()` always computes its own per-unit-clustered
version. Fits made with `add_ridge = TRUE` are **not** supported (the
ridge-augmented design and post-fit un-shrink are outside the validated
construction); the accessor errors on them.

## Assumptions

The uniformly-valid SE relies on the hypotheses of paper Theorem
`debiased.att.thm`:

- **Marginal cohort assignment / (Psi-IF).** Cohort membership is
  independent of the outcome noise, so the regression and cohort-weight
  variance channels are uncorrelated and add. The package's default
  cohort weight estimator (sample proportions) satisfies this. Under
  *confounded* assignment (cohort membership correlated with
  treatment-effect heterogeneity) the two channels correlate and the
  additive SE can mis-state uncertainty; a combined per-unit score is
  then required (not yet implemented).

- **Correctly specified conditional mean** of the outcome.
  Neyman-orthogonality removes first-order sensitivity to the selected
  nuisance coefficients, but not to mean-misspecification. The
  regression channel is the unit-clustered sandwich, valid as the number
  of units `N -> infinity` with **no** model for within-unit dependence
  (so a `gls = FALSE` fit, which skips GLS whitening, yields a valid
  cluster-robust SE; \#307, paper Decision D1). GLS whitening
  (`gls = TRUE`) buys *efficiency* (an asymptotically smaller
  `var_reg`), not validity.

- **(Low-dimensional `p < NT` only) Asymptotically negligible ridge,**
  `lambda = o((NT)^(-1/2))` (the theoretical condition; the leading case
  is the exact inverse `lambda = 0`). The implementation does **not**
  use a vanishing schedule — it adds a fixed numerical stabilizer
  `1e-6 * mean(diag(Sigma))` to the Gram before solving, which is
  negligible and cancels in the OLS-identity case.

- **Growing number of clusters,** `N -> infinity`. The CLT is over the
  `N` independent units, not the `NT` rows. With few treated units the
  cluster approximation is poor and the interval can under-cover; a
  genuine few-clusters correction (a restricted wild-cluster bootstrap-t
  or a CR2-type analytic adjustment, \#361) is future work. The
  `method = "bootstrap"` option does *not* remedy this (see its
  documentation).

- **Regularity / two regimes.** When `p < NT` (Theorem
  `debiased.att.thm`) the full design Gram is nonsingular and the
  debiasing direction is the exact inverse; the accessor reduces to
  debiased ETWFE. When `p >= NT` (the high-dimensional FETWFE theory)
  the Gram is singular and the direction is the **nodewise
  (desparsified-lasso) relaxed inverse** of equation
  `debiased.highdim.v` (the same estimate and SE, only `v` differs —
  "one estimator, two regimes"); it requires sparsity
  (`s_N log(p_N) / sqrt(N) -> 0`), a restricted-eigenvalue condition
  over sparse cones, `||v*||_1 = O(1)`, and a bounded limiting variance
  `a' Sigma_theta^(-1) a`. The high-dimensional branch realizes
  **Theorem `debiased.highdim.thm`** (Faletto 2025) with the
  theory-scaled penalty; its overall-ATT interval is **uniformly
  valid**, with coverage validated near-nominally in simulation (0.94 at
  the `p >= NT` anchor, robust across the nodewise penalty scales
  studied). Inspect the returned `feasibility` / `converged` diagnostics
  to confirm the per-fit nodewise directions are well-behaved.

The two regimes use different nuisances. The **fixed-`p`** branch uses
the fit's bridge (`q < 1`) `theta_hat`: Theorem `debiased.att.thm` needs
only *consistency* of the nuisance, and with the exact inverse this is
the OLS identity. The **high-dimensional** branch instead fits an
internal `q = 1` fused lasso
(`grpreg::cv.grpreg(penalty = "gBridge", gamma = 1)`, CV-selected
`lambda`) as the nuisance, because the high-dimensional FETWFE theory
controls the orthogonalization remainder through the nuisance's `l1`
*rate* (`||theta_hat - theta*||_1 = O_p(s_N lambda_theta)`) — the rate
the `q = 1` fused lasso enjoys under a restricted-eigenvalue condition
(the standard desparsified-lasso rate; cf. Negahban et al. 2012), the
`p >= NT` extension the paper points to; the `q < 1` bridge is
super-efficient / non-uniform — the wrong nuisance for a uniformly-valid
CI (#303). The *input* fit must still be `q < 1` (it supplies the
cohort-weight variance channel `att_var_2`); only the high-dimensional
debiasing nuisance is `q = 1`. That nuisance's cross-validation uses a
fixed, data-derived seed, so it is deterministic across calls
(`debiasedATT()` and the high-dimensional
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
band center agree exactly) — but it is not tunable via the fit's
`cv_seed`.

## References

Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions. *arXiv preprint
arXiv:2312.05985*. <https://arxiv.org/abs/2312.05985>.

## See also

[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
for the fused fit;
[`cohortTimeATTs()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortTimeATTs.md)
/
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
/
[`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md)
for the disaggregated fused effects.

## Examples

``` r
if (FALSE) { # \dontrun{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2, seed = 1)
  dat <- simulateData(coefs, N = 200, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 1)
  fit <- fetwfeWithSimulatedData(dat, q = 0.5)
  fit$att_hat # fused (pointwise)
  debiasedATT(fit) # debiased (uniformly valid)
} # }
```
