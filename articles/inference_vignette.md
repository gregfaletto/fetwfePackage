# Inference in fetwfe: standard errors, simultaneous bands, and high-dimensional debiasing

``` r

library(fetwfe)
```

This vignette documents inference in the package: how
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md),
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
and
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
compute their standard errors (and the three values `se_type` can take),
the family-wise simultaneous confidence bands from
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md),
and the high-dimensional ($`p \geq NT`$) debiased path of
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md).
The four estimators share the same inferential machinery, so the
discussion applies to all of them; we use
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
in the running example.

## 1. The headline: tight Gaussian CIs by default under (Ψ-IF)

By default, the package’s standard errors — the `att_se` slot on the
returned object and the entries of `catt_ses` — are the *tight Gaussian
variance* from the paper’s Theorem `te.asym.norm.thm`(c$`'`$), valid
under the influence-function condition (Ψ-IF) (Faletto 2025,
[arXiv:2312.05985](https://arxiv.org/abs/2312.05985), Assumption (Ψ-IF)
at paper Eq. `psi.if.assum`). The combined ATT variance has two pieces,

``` math
V_1 \;=\; \sigma_\varepsilon^2 \, \boldsymbol{\alpha}^\top \boldsymbol{\Sigma}^{-1} \boldsymbol{\alpha},
\qquad
V_2 \;=\; v_\psi(\boldsymbol{\beta}_0, \mathcal{S}),
```

where $`V_1`$ is the Kock-2013 regression-coefficient variance
contribution and $`V_2`$ is the propensity-score variance contribution;
the names $`V_1`$ and $`V_2`$ are the paper’s, and the package exposes
both of them by name on the fitted object under
`result$internal$variance_components`. Theorem (c$`'`$) shows that under
(Ψ-IF) the same-data confidence interval is asymptotically Gaussian with
the **tight** variance $`V_1 + V_2`$, *identical to the independent-data
variance* in part (b) of the theorem. Concretely, the package returns

``` math
\widehat{\text{SE}}(\hat{\tau}_{\text{ATT}})
\;=\;
\sqrt{\widehat{\text{Var}}_1 + \widehat{\text{Var}}_2}
\;=\;
\sqrt{\tilde{v}_N \,/\, N},
```

where $`\tilde{v}_N := \hat{v}_N / T`$ is the unit-scaled paper-notation
variance estimator (paper line 2004), and the Wald confidence interval
is

``` math
\text{CI}_N(\tilde{v}_N;\,\alpha)
\;=\;
\left[\,
\hat{T}_N \;\pm\; z_{1 - \alpha/2}\,\sqrt{\tilde{v}_N \,/\, N}
\,\right]
```

per paper Eq. `conf.int.form` (line 1259). The same combined variance
backs the per-cohort entries of `catt_ses` and the `catt_df`
confidence-interval columns.

### 1.1 Which propensity-score estimators satisfy (Ψ-IF)?

Assumption (Ψ-IF) requires the cohort-weight estimator
$`\widehat{\boldsymbol{\psi}}_N`$ to admit an i.i.d.-across-$`i`$
unit-level influence-function representation. It is satisfied by

- **Cohort sample proportions** $`\widehat{\pi}_g = N_g / N`$ — the
  package’s *default* propensity-score estimator. So out of the box the
  tight Gaussian CI is asymptotically exact.
- **Multinomial-logit MLE** for $`P(W = g \mid \boldsymbol{X})`$.
- **Any GLM** fit by maximum likelihood on $`W \mid \boldsymbol{X}`$.
- **Kernel or series regression** of $`\mathbf{1}\{W = g\}`$ on
  $`\boldsymbol{X}`$.

(See the discussion at paper line 1268: “When (Ψ-IF) holds — which is
the case for every standard parametric or semi-parametric
propensity-score estimator listed in Remark `psi.if.coverage` (cohort
sample proportions, multinomial logit, any GLM, kernel/series
regression) — Theorem `te.asym.norm.thm`(c$`'`$) gives asymptotic
Gaussianity with the tight variance $`V_1 + V_2`$.”)

The (Ψ-IF) condition fails for **Robins-Rotnitzky-style augmented
doubly-robust estimators** that augment the propensity score with
outcome residuals. The package does not currently implement these; when
it does, users will be able to opt out of the tight CI via
`se_type = "conservative"` (see Section 3).

### 1.2 What changed in v1.12.0

Versions $`\le`$ 1.11.7 of the package returned the *conservative*
Cauchy-Schwarz upper bound

``` math
\widehat{\text{SE}}(\hat{\tau}_{\text{ATT}})
\;=\;
\sqrt{\widehat{\text{Var}}_1 + \widehat{\text{Var}}_2 + 2\sqrt{\widehat{\text{Var}}_1\,\widehat{\text{Var}}_2}}
```

as the same-data default (from paper Theorem (c)). Theorem (c$`'`$) was
added in the paper’s same-data Gaussianity upgrade, and v1.12.0 switches
the package default to match: the tight Gaussian variance is now the
headline number. The conservative formula remains available as an opt-in
fallback (Section 3) for any users whose propensity-score estimator
violates (Ψ-IF).

## 2. The two variance pieces and how the package exposes them

The combined ATT variance decomposes as $`V_1 + V_2`$, with both pieces
estimated by the package and exposed as named slots.

### 2.1 $`V_1`$: regression-coefficient variance

The regression-coefficient piece is

``` math
\widehat{\text{Var}}_1
\;=\;
\frac{\sigma_\varepsilon^2}{NT}\;\psi_{\text{att}}^\top \widehat G^{-1} \psi_{\text{att}},
```

where $`\widehat G^{-1}`$ is the Gram inverse on the bridge-selected
support and $`\psi_{\text{att}}`$ encodes the cohort-weighted ATT
contrast. The per-cohort SEs (`catt_ses`) have the same form with a
cohort selector $`\psi_g`$ in place of $`\psi_{\text{att}}`$. This piece
is stored on the fitted object as
`result$internal$variance_components$att_var_1`; the paper-notation
per-unit-scaled version is `result$internal$variance_components$V_1`
(with `V_1 = N * att_var_1`).

### 2.2 $`V_2`$: propensity-score variance

The propensity-score piece accounts for variability from estimating the
cohort-membership probabilities $`\widehat{\pi}_g`$. It scales like
$`1/N`$ and is unrelated to the regression residuals. This piece is
stored as `result$internal$variance_components$att_var_2`, with
paper-notation per-unit-scaled `V_2 = N * att_var_2`.

### 2.3 The combined variance: `tilde_v_N` and the paper’s catalogue

The unit-scaled total variance is

``` math
\tilde{v}_N \;=\; N \cdot \widehat{\text{Var}}(\hat{T}_N),
```

stored as `result$internal$variance_components$tilde_v_N`. The original
(un-scaled) paper notation `hat_v_N = T * tilde_v_N` is also exposed.

The paper catalogues six variance estimators at line 2006; the package
exposes the four that apply to the ATT (the other two apply to the
CATT(x), which connects to the on-hold
[`predict()`](https://rdrr.io/r/stats/predict.html) work in issue \#33):

| Slot | Paper Eq. | Meaning |
|----|----|----|
| `tilde_v_N_C` | `v.n.r.t.att.const` | Fixed-π exact (no propensity-score noise) |
| `tilde_v_N_C_pi_hat` | `v.n.r.t.att.rand` | Random-π exact (tight Gaussian under (Ψ-IF)) |
| `tilde_v_N_C_pi_hat_cons` | `v.n.r.t.att.rand.cons` | Random-π conservative (Cauchy-Schwarz) |
| `tilde_v_N_cons` | `var.est.kock.wooldridge.subgauss.cons` | Subgaussian conservative bound |

The default `att_se` value is computed from `tilde_v_N_C_pi_hat` (via
`sqrt(tilde_v_N_C_pi_hat / N)`); `se_type = "conservative"` switches the
headline to `tilde_v_N_C_pi_hat_cons`. When `indep_counts` is supplied
(two-sample regime), the two-sample-exact formula `tilde_v_N_C_pi_hat`
is in force regardless of `se_type`.

## 3. The conservative fallback: `se_type = "conservative"`

If a user’s propensity-score estimator violates (Ψ-IF) — the canonical
example being a Robins-Rotnitzky-style augmented doubly-robust estimator
that injects outcome-residual information into the cohort weights —
Theorem (c$`'`$) no longer applies. The paper’s fallback is Theorem
(c)’s subgaussian result with the **conservative** variance bound

``` math
\widehat{\text{Var}}^{\text{(cons)}}(\hat{\tau}_{\text{ATT}})
\;=\;
\widehat{\text{Var}}_1 + \widehat{\text{Var}}_2 + 2\sqrt{\widehat{\text{Var}}_1\,\widehat{\text{Var}}_2}.
```

This is the Cauchy-Schwarz upper bound on $`V_1 + V_2`$ that holds
without any influence-function structure on the cohort-weight estimator.
Until v1.11.7 the package returned this bound as the default; v1.12.0
demoted it to an opt-in fallback via

``` r

fetwfe(..., se_type = "conservative")
```

The package does not currently implement any non-(Ψ-IF) propensity-score
estimator, so for users running the default cohort-sample-proportions
estimator there is no practical reason to pass
`se_type = "conservative"`. The argument exists for users who supply
their own propensity-score estimates (see issue \#32) or who want to
reproduce the pre-v1.12.0 conservative default for backward
compatibility.

A worked comparison on a clean F1-conforming panel:

``` r

set.seed(2026)
sim_coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
# sim_coefs was built without a seed, so draw from the ambient generator
# (seeded just above); seed = NA selects the ambient generator silently.
sim_data  <- simulateData(
  sim_coefs,
  N = 120,
  sig_eps_sq = 1,
  sig_eps_c_sq = 0.5,
  seed = NA
)

res_default <- fetwfe(
  pdata = sim_data$pdata,
  time_var = sim_data$time_var,
  unit_var = sim_data$unit_var,
  treatment = sim_data$treatment,
  covs = sim_data$covs,
  response = sim_data$response,
  sig_eps_sq = sim_data$sig_eps_sq,
  sig_eps_c_sq = sim_data$sig_eps_c_sq
)
res_cons <- fetwfe(
  pdata = sim_data$pdata,
  time_var = sim_data$time_var,
  unit_var = sim_data$unit_var,
  treatment = sim_data$treatment,
  covs = sim_data$covs,
  response = sim_data$response,
  sig_eps_sq = sim_data$sig_eps_sq,
  sig_eps_c_sq = sim_data$sig_eps_c_sq,
  se_type = "conservative"
)

c(
  default      = res_default$att_se,
  conservative = res_cons$att_se,
  ratio_cons_over_default = res_cons$att_se / res_default$att_se
)
#>                 default            conservative ratio_cons_over_default 
#>               0.1050007               0.1469506               1.3995195
```

The conservative SE is strictly wider; the tight Gaussian SE is the
asymptotically-exact one under the package’s default
cohort-sample-proportions estimator.

## 4. Experimental: cluster-robust standard errors via `se_type = "cluster"`

The package’s default standard errors rely on the model-based covariance
structure (paper Assumption F1: i.i.d. units, idiosyncratic noise plus a
unit random effect). In applied DiD work it is not unusual to suspect
one or more violations of F1:

- **Within-unit serial correlation beyond the random effect.** The unit
  random effect $`c_i`$ in F1 absorbs a single, time-constant deviation
  per unit. It does *not* model serial correlation in the idiosyncratic
  shocks. Bertrand, Duflo, and Mullainathan (2004) is the classic
  warning that ignoring within-unit serial correlation can drastically
  understate DiD standard errors.
- **Heteroskedasticity across units.** F1 imposes a single
  $`\sigma_\varepsilon^2`$ across all units.
- **Higher-level clustering.** F1 treats units as i.i.d.; if
  observations within a state-year or industry-year share unobserved
  shocks, the unit-level F1 variance will understate the true sampling
  variability.

Starting in version 1.6.0, all four estimators accept an experimental
`se_type` argument with value `"cluster"`:

``` r

fetwfe(..., se_type = "cluster")
```

Setting `se_type = "cluster"` swaps the model-based
regression-coefficient variance $`V_1`$ for a **unit-clustered
Liang-Zeger CR1 sandwich** computed on the bridge-selected support. The
$`V_2`$ piece and the (Ψ-IF)-tight combination $`V_1 + V_2`$ are
unchanged — only the $`V_1`$ formula switches.

### 4.1 The formula

Let $`\widehat{S}`$ be the support selected by the bridge regression,
$`X_{\widehat{S}}`$ the corresponding design matrix in the coordinate
system the regression was solved in (GLS-transformed for ETWFE/twfeCovs,
fusion-then-GLS-transformed for FETWFE/BETWFE), and
$`\widehat\varepsilon`$ the residuals from OLS on that selected support.
The cluster-robust variance is

``` math
V_{\text{CR}}
\;=\;
\frac{N}{N-1}\;
(X_{\widehat{S}}^\top X_{\widehat{S}})^{-1}\;
\left(\sum_{i=1}^N X_{i\cdot\widehat{S}}^\top \widehat\varepsilon_{i\cdot} \widehat\varepsilon_{i\cdot}^\top X_{i\cdot\widehat{S}}\right)\;
(X_{\widehat{S}}^\top X_{\widehat{S}})^{-1},
```

with units $`i = 1, \dots, N`$ as clusters and an $`N/(N-1)`$
small-sample adjustment (matching
`sandwich::vcovCL(cadjust = TRUE, type = "HC0")`). The CATT SE for
cohort $`g`$ is $`\sqrt{\psi_g^\top V_{\text{CR}} \psi_g}`$ (using a
zero-padded $`\psi_g`$ on the full selected support); the ATT
regression-coefficient variance is
$`\psi_{\text{att}}^\top V_{\text{CR}} \psi_{\text{att}}`$, replacing
$`V_1`$ above. The propensity-score piece $`V_2`$ is unchanged because
it depends on empirical cohort proportions, not regression residuals.

For FETWFE and BETWFE, `se_type = "cluster"` is only meaningful when
`q < 1` (the bridge oracle property is required); for `q >= 1` the
cluster path returns `NA` just like the default. ETWFE and `twfeCovs`
have no `q` argument, so the cluster path always runs when the Gram
matrix is invertible.

### 4.2 Why we call this experimental

The CR1 sandwich is a textbook estimator and the package’s
implementation matches `sandwich::vcovCL()` to numerical precision on a
clean panel without selection. What is *not* yet covered by the paper’s
theory is:

- Verification that the bridge oracle property of Theorem 6.2 still
  holds under the relaxed covariance structure that motivates
  cluster-robust SEs in the first place.
- Sandwich consistency *after model selection*, i.e., that the CR1
  sandwich evaluated at the bridge-selected support is a consistent
  estimator of the true asymptotic variance.

These extensions are mechanically routine but conceptually non-trivial,
and they are explicitly outside the package’s current scope. Until they
land, `se_type = "cluster"` is exposed as an opt-in, clearly-labelled
experimental feature.

**Recommendation.** Until the theory lands, treat `se_type = "cluster"`
as a sensitivity check: report both `se_type = "default"` and
`se_type = "cluster"` in applied work, comment on the gap, and lean on
the default for headline numbers.

### 4.3 Worked cluster comparison

``` r

res_cluster <- fetwfeWithSimulatedData(sim_data, se_type = "cluster")

c(
  default = res_default$att_se,
  cluster = res_cluster$att_se
)
#>   default   cluster 
#> 0.1050007 0.1000234
```

On this F1-conforming simulated panel the two SEs are similar by
construction: the data-generating process satisfies F1, so the
model-based SE is already valid and the cluster-robust SE estimates the
same underlying variance. Under a deliberately serially-correlated DGP
(or under heteroskedasticity, or higher-level clustering) the
cluster-robust SE will typically be larger.

The [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods label the SE
so it is clear which one was used:

``` r

print(res_cluster)
#> Fused Extended Two-Way Fixed Effects Results
#> ===========================================
#> 
#> Overall Average Treatment Effect (ATT):
#>   Estimate:   -0.7798
#>   Std. Error (cluster-robust): 0.1000
#>   P-value:    6.38e-15
#>   Selected:   TRUE
#>   95% CI:    [-0.9758, -0.5838]
#> 
#> Cohort Average Treatment Effects (CATT) [simultaneous 95% CI]:
#>  cohort   estimate         se     ci_low    ci_high     p_value selected
#>       2 -0.4137526 0.06057451 -0.5489811 -0.2785242 1.69269e-11     TRUE
#>       3 -1.6817852 0.13059071 -1.9733200 -1.3902505 0.00000e+00     TRUE
#>       4  0.0000000 0.00000000  0.0000000  0.0000000          NA    FALSE
#> 
#> Event-Study Average Treatment Effects (per event time) [simultaneous 95% CI]:
#>  event_time n_cohorts   estimate        se     ci_low    ci_high      p_value
#>           0         3  0.0000000 0.0000000  0.0000000  0.0000000           NA
#>           1         3 -0.6359733 0.1104477 -0.8975561 -0.3743905 2.162509e-08
#>           2         3 -0.6359733 0.1104477 -0.8975561 -0.3743905 2.529819e-08
#>           3         2 -3.0326244 0.1919425 -3.4872181 -2.5780307 0.000000e+00
#>           4         1  0.1684927 0.2513507 -0.4268026  0.7637881 8.581460e-01
#> 
#> Model Details:
#>   Units (N)           : 120
#>   Time periods (T)    : 6
#>   Treated cohorts (G) : 3
#>   Covariates (d)      : 2
#>   Features (p)        : 62
#>   Selected size       : 26
#>   Lambda*             : 0.0114
```

The CATT SEs and confidence intervals in `catt_df` are recomputed from
the same cluster-robust sandwich; the CATT p-values follow accordingly.

## 5. Choosing the bridge penalty parameter `lambda`

As of v1.13.0,
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
and
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
select the bridge penalty `lambda` via 10-fold cross-validation by
default (using
[`grpreg::cv.grpreg`](https://pbreheny.github.io/grpreg/reference/cv.grpreg.html)).
The prior default was BIC over the same `grpreg` lambda grid. The change
addresses a finite-sample bias issue documented in simulation studies
(see issue \#164 in the package repository): under the prior BIC
default, the overall ATT estimator was biased toward zero at moderate
sample sizes, producing 95% confidence intervals whose empirical
coverage was as low as 0.00 in some regimes. Cross-validation restores
near-nominal coverage in every regime tested and additionally has the
lowest MSE in the high-dim regime FETWFE was originally designed for.

The `lambda_star` and `lambda_star_model_size` result slots are
populated identically under both selection paths; only the underlying
selection criterion differs. The CV path stores its provenance in three
new top-level slots: `lambda_selection`, `cv_folds`, and `cv_seed`. The
seed defaults to `as.integer(N * T)` so consecutive calls on the same
dataset are reproducible without the user having to specify a seed.

To recover the prior BIC behavior — for reproducing results from
analyses run against v1.12.0 or earlier, for example — pass
`lambda_selection = "bic"`:

``` r

fit_bic <- fetwfe(..., lambda_selection = "bic")
```

To control the CV fold assignment explicitly, pass `cv_folds` and/or
`cv_seed`:

``` r

fit_cv <- fetwfe(..., cv_folds = 20L, cv_seed = 42L)
```

## 6. Simultaneous confidence intervals via `simultaneousCIs()`

The package’s per-point standard errors and Wald intervals deliver
`(1 - alpha)` coverage *per effect*: for any single cohort ATT, or any
single event-time coefficient, the CI contains the true value with
probability approximately `1 - alpha`. But when you make a claim
involving multiple effects jointly — for example, “no event-time
coefficient is significantly positive at the 5% level” — per-effect
coverage does not imply family-wise coverage. The probability that
*every* coefficient in a family of `K` effects lies in its per-point CI
can be far below `1 - alpha` when `K` is large.

[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
returns intervals that maintain family-wise coverage at the requested
`1 - alpha`, using the exact joint multivariate-normal critical value
$`c_{1-\alpha}`$ defined as the $`(1 - \alpha)`$ quantile of
$`\max_k |Z_k|`$ where $`Z \sim N(0, \rho)`$ and $`\rho`$ is the
correlation matrix of the effect family’s joint covariance. Under
Theorem (c$`'`$) (`paper_arxiv.tex:1233`) and Assumption
$`(\Psi\text{-IF})`$ (assumption equation `paper_arxiv.tex:2013`;
in-prose discussion of when it implies tight Gaussianity at paper line
1268), this joint Gaussianity is exact asymptotically; under the paper’s
fixed-dim framing (AE 1(d)), no high-dimensional correction is needed.

For an event-study family on the fitted `res_default` object from the
previous sections:

``` r

sci <- simultaneousCIs(res_default, family = "event_study", alpha = 0.05)
print(sci)
#> Parametric simultaneous (1 - alpha) confidence intervals
#> Family: event_study, K = 5, alpha = 0.05
#> Simultaneous critical value: 2.3693 (pointwise 1.9600, Bonferroni 2.5758)
#> 
#>  effect   estimate simultaneous_ci_low simultaneous_ci_high pointwise_ci_low
#>      e0  0.0000000           0.0000000            0.0000000        0.0000000
#>      e1 -0.6726640          -0.9531319           -0.3921962       -0.9046805
#>      e2 -0.6726640          -0.9531319           -0.3921962       -0.9046805
#>      e3 -3.0092312          -3.4694847           -2.5489777       -3.3899750
#>      e4  0.1684927          -0.4085822            0.7455677       -0.3088914
#>  pointwise_ci_high
#>          0.0000000
#>         -0.4406476
#>         -0.4406476
#>         -2.6284874
#>          0.6458769
```

The function returns a `simultaneous_cis` S3 object with `$ci` (the
interval data frame), `$critical_value` (the simultaneous
$`c_{1-\alpha}`$), `$pointwise_critical_value` ($`z_{1-\alpha/2}`$, for
reference), and `$bonferroni_critical_value` ($`z_{1-\alpha/(2K)}`$, the
conservative family-wise upper bound). When effects are positively
correlated — as they typically are in DiD, because event-study
coefficients share the regression-coefficient variance piece — the
simultaneous critical value sits strictly between the pointwise and
Bonferroni values, giving tighter intervals than Bonferroni while
preserving family-wise coverage.

The four available family options are `"event_study"` (one effect per
post-treatment event time), `"cohort"` (one effect per treated cohort),
`"all_post_treatment"` (one effect per `(g, t)` cell), and `"custom"` (a
user-supplied `K x num_treats` contrast matrix following the
`multcomp::glht()` convention; see
[`?simultaneousCIs`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)).

The critical value is computed via
`mvtnorm::qmvnorm(..., algorithm = mvtnorm::GenzBretz())` (mvtnorm’s
default quasi-Monte Carlo integrator), sub-second on a typical fit for
any realistic FETWFE family (no `K` cap of practical concern). The
function is deterministic in its inputs: the same fit plus the same
family plus the same `alpha` plus the same `contrasts` produces the same
critical value across calls. This is achieved by save/restoring the
caller’s `.Random.seed` and using a fixed internal seed immediately
before the `qmvnorm()` call; the function does not mutate the caller’s
RNG state.

As of version 1.16.0 the `mvtnorm` package is an `Imports` dependency
(it ships with the package), because simultaneous bands are now computed
by default at fit time — see the `ci_type` subsection below. (In
versions 1.15.x `mvtnorm` was in `Suggests` and only loaded when
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
was called.) The `plot.simultaneous_cis` method additionally requires
`ggplot2`, which is in `Suggests`; install it with
`install.packages("ggplot2")` if you have not already.

The default analytic path matches `did::aggte(cband = TRUE)` from the
Callaway-Sant’Anna `did` package, which computes simultaneous bands via
the CCK (2013) multiplier bootstrap. Asymptotically the two critical
values are identical under fixed-dim Gaussianity; the analytic
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
path is deterministic and sub-second, whereas a multiplier bootstrap is
stochastic and resampling-based.

### The multiplier-bootstrap method: `method = "bootstrap"`

`simultaneousCIs(..., method = "bootstrap")` computes the same sup-t
critical value by the CCK (2013) multiplier bootstrap directly — the
procedure `did` uses. It perturbs the per-unit influence functions
$`\hat f_i^{(k)}`$ with random multipliers $`\xi_i`$
(`multiplier = "rademacher"` or `"mammen"`), forms
$`T^{(b)} = \max_k |\frac{1}{N}\sum_i \xi_i \hat f_i^{(k)}| / \widehat{\mathrm{sd}}_k`$
over $`B`$ replicates, and takes the $`(1-\alpha)`$ quantile. The
analytic and bootstrap critical values agree up to Monte-Carlo error,
and both sit between the pointwise and Bonferroni values:

``` r

sci_boot <- simultaneousCIs(res_default, family = "all_post_treatment",
                            method = "bootstrap", B = 1000, seed = 1)
c(pointwise  = sci_boot$pointwise_critical_value,
  analytic   = simultaneousCIs(res_default, family = "all_post_treatment")$critical_value,
  bootstrap  = sci_boot$critical_value,
  bonferroni = sci_boot$bonferroni_critical_value)
#>  pointwise   analytic  bootstrap bonferroni 
#>   1.959964   2.479436   2.374132   2.865260
```

Two reasons to prefer the bootstrap. First, it **scales to large effect
families**: `qmvnorm` integration strains at hundreds of effects
(e.g. simultaneous bands over the full disaggregated cohort $`\times`$
event-time grid), whereas the bootstrap stays stable and linear in
$`B \cdot K`$. Second, it is **heteroskedasticity/cluster-robust by
construction** — its per-unit influence matrix $`F`$ reproduces the
cluster-robust joint covariance ($`\mathrm{crossprod}(F)/(NT)^2`$ equals
it up to the $`N/(N-1)`$ finite-sample factor), so on a
`se_type = "cluster"` fit the bootstrap per-effect standard errors equal
the analytic ones exactly (and on a default fit the bootstrap reports
the cluster-robust band regardless). Pass an integer `seed` for
reproducible draws; with `seed = NULL` the draws come from the ambient
RNG.

`method = "bootstrap"` supports all four families. For
`family = "event_study"` the bootstrap additionally perturbs a **second,
independent** multiplier stream over a per-unit cohort-probability
(propensity) influence function, alongside the regression-channel stream
— a two-channel bootstrap reproducing the analytic variance
decomposition $`\Sigma = \Sigma_1 + \Sigma_2`$. The event-study family
is the one whose propensity term $`\Sigma_2`$ is non-zero, because each
event-time effect pools across the cohorts treated by that event time,
weighted by the estimated cohort probabilities $`\hat\pi_g = N_g / N`$;
its per-unit propensity influence function is the
cohort-sample-proportion summand $`\xi_i = e_{W_i} - \hat\pi`$ (paper
assumption $`(\Psi\text{-IF})`$). The two channels are asymptotically
uncorrelated (paper Theorem C.1, Step 5), so their variances add with no
cross-term — which is exactly why the multipliers are drawn
independently. On a `se_type = "cluster"` fit the event-study bootstrap
per-effect standard errors equal the analytic ones to machine precision,
and the critical value matches the analytic `qmvnorm` one up to
Monte-Carlo error:

``` r

es_boot <- simultaneousCIs(res_default, family = "event_study",
                           method = "bootstrap", B = 1000, seed = 1)
c(pointwise  = es_boot$pointwise_critical_value,
  analytic   = simultaneousCIs(res_default, family = "event_study")$critical_value,
  bootstrap  = es_boot$critical_value,
  bonferroni = es_boot$bonferroni_critical_value)
#>  pointwise   analytic  bootstrap bonferroni 
#>   1.959964   2.369258   2.355859   2.575829
```

**High-dimensional $`p \geq NT`$.** When the (full) design is
high-dimensional — where the analytic Gram inverse need not exist, so
there is no analytic sup-t band — the bootstrap automatically switches
to the full-design **desparsified** construction of
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
(the high-dimensional FETWFE theory): each effect’s debiasing direction
becomes a nodewise (`riesz_lasso`) relaxed inverse over the singular
Gram, with the bridge residuals. This is *uniformly valid* (not
post-selection) and is the only route to simultaneous bands in that
regime. It is **experimental**
([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
fits only): its underlying overall-ATT coverage is validated
near-nominally at the `p > NT` anchor of Faletto (2025), but the
family-wise band coverage is not itself simulation-validated, so the
returned object carries a `regime` field and per-effect `feasibility` /
`converged` / `lambda_node` diagnostics, and `lambda_c` tunes the
nodewise penalty `lambda_node = lambda_c * max(|a|) * sqrt(log(p) / N)`
— too-small values leave the directions infeasible and trigger a
warning. `lambda_c` accepts either a fixed positive number (default the
theory scale `1.0`) or the string `"cv"`, which selects the constant
**per fit** by cross-validating the Riesz loss over the KKT-feasible
region (the desparsified-lasso / auto-DML standard), falling back to
`1.0` when no grid penalty is feasible; one constant then serves both
the
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
point estimate and every bootstrap band effect (issue \#295). The
default stays the fixed `1.0` (the opt-in CV selector picks the
feasibility-appropriate constant; \#88 validated overall-ATT coverage
with it at the studied anchor). The simulator can generate the `p > NT`
panels this needs (see the
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)
documentation). `family = "event_study"` uses this desparsified path too
in the high-dimensional regime, additionally carrying the per-unit
propensity channel (the two-channel $`\Sigma_1 + \Sigma_2`$ bootstrap).
In the high-dimensional regime the band is **centered on the debiased
estimate** (the high-dimensional FETWFE theory correction, matching
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md))
rather than the post-selection bridge estimate, which would otherwise
carry the selection bias. The desparsified path is
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)-only;
a
non-[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
$`p \geq NT`$ fit
(e.g. [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md))
has no debiased band, so it falls back to the fixed-$`p`$
selected-support band (rather than erroring) but emits a
[`warning()`](https://rdrr.io/r/base/warning.html). A coverage study
(issue \#308) shows that fallback band substantially **under-covers** in
the $`p \geq NT`$ regime — even when the selected support is
low-dimensional — because the bridge shrinks the band center toward zero
(so it is biased away from the truth) and the post-selection standard
error understates the sampling variability; treat such a band as
unreliable and use a
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
fit for valid high-dimensional simultaneous bands.

### Simultaneous bands are the default: the `ci_type` argument

As of version 1.16.0, the confidence-interval bounds that
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
/
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
/
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)
/
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md)
*return and display* are simultaneous by default — matching `did`’s
`cband = TRUE` default. A new fit-time argument
`ci_type = c("simultaneous", "pointwise")` controls this:

- `ci_type = "simultaneous"` (the default) makes the returned
  `catt_df$ci_low` / `catt_df$ci_high` the **cohort-family**
  simultaneous bounds, and
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)’s
  `ci_low` / `ci_high` the **event-study-family** simultaneous bounds —
  the very bounds `simultaneousCIs(fit, family = "cohort")` /
  `simultaneousCIs(fit, family = "event_study")` return.
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) display these
  (labeled `[simultaneous 95% CI]`), and
  [`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) (on
  the fit, on
  [`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md),
  and on
  [`cohortStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/cohortStudy.md))
  surfaces them too.
- `ci_type = "pointwise"` restores the pre-1.16.0 per-effect Wald
  intervals everywhere (the bounds you would get from
  `estimate ± qnorm(1 - alpha/2) * se`).

The interval bounds **and** the per-cohort p-values (`p_value`) follow
`ci_type`. Under `"simultaneous"` the `p_value` is the single-step max-T
multiplicity-adjusted (family-wise) p-value matching the band; under
`"pointwise"` it is the per-effect Wald p-value. The two correspond
exactly: a coefficient lies outside its `1 - alpha` simultaneous band
iff its adjusted `p_value` is `< alpha` (version 1.18.0, \#200). The
standard errors (`se`), selection flags (`selected`), and the
overall-ATT confidence interval and p-value (a single scalar, `K = 1`,
with no family and hence no widening) are **identical** under both
settings — for the bounds the only difference is the critical-value
multiplier (`qnorm(1 - alpha/2)` for pointwise vs the joint
$`c_{1-\alpha}`$ for simultaneous).

``` r

fit_sim <- fetwfeWithSimulatedData(sim_data) # ci_type = "simultaneous" (default)
fit_pw <- fetwfeWithSimulatedData(sim_data, ci_type = "pointwise")

# Same point estimates and SEs; only the bounds differ (simultaneous wider).
cbind(
  estimate = fit_sim$catt_df$estimate,
  se = fit_sim$catt_df$se,
  simul_width = fit_sim$catt_df$ci_high - fit_sim$catt_df$ci_low,
  pointwise_width = fit_pw$catt_df$ci_high - fit_pw$catt_df$ci_low
)
#>        estimate         se simul_width pointwise_width
#> [1,] -0.4137526 0.07754859   0.3463267       0.3039849
#> [2,] -1.6817852 0.12781952   0.5708332       0.5010433
#> [3,]  0.0000000 0.00000000   0.0000000       0.0000000
```

Because the tidiers pass the object’s bounds straight through, every CI
surface in the package agrees under the default: `broom::tidy(fit)`’s
cohort rows, `summary(fit)`, `cohortStudy(fit)`, and
`plot(fit, type = "catt")` all report the same simultaneous bounds. (The
overall-ATT row of `broom::tidy(fit)` stays the scalar pointwise
interval.) When standard errors are unavailable (`q >= 1` for the bridge
estimators, or a rank-deficient design) the bounds are `NA` under both
settings, and the fit still succeeds — the simultaneous path degrades
cleanly to the (`NA`) pointwise bounds rather than throwing an error.

## 7. Selection consistency as an implicit zero-effect test

A question applied users ask constantly is: *is this cohort’s treatment
effect distinguishable from zero — and can I conclude that it actually
**is** zero?* For
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
(and
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md)),
the answer is built into the estimator. The bridge penalty’s
*restriction selection consistency* (Faletto 2025, Theorem 6.2) makes
the `selected` flag on `catt_df` an **implicit hypothesis test** of
$`H_0: \tau_{\text{ATT}}(g) = 0`$, for free, with no extra computation.
The mapping is:

- `selected = TRUE` (estimate $`\neq 0`$) is the asymptotic analogue of
  **reject $`H_0`$**: the effect is nonzero.
- `selected = FALSE` (estimate exactly $`0`$) is the asymptotic analogue
  of **fail to reject $`H_0`$** — and, under Theorem 6.2, it is the
  *stronger* statement that the truth is zero with probability tending
  to one.

The conclusion for a selected-out cohort is delivered by the `selected`
flag, **not** by a p-value. A cohort the penalty zeros out has
`estimate == 0`, `se == 0`, and `p_value == NA` by construction — there
is no Wald p-value to read, and the inferential content lives entirely
in `selected`. (A condensed, first-time-user-oriented version of this
discussion appears under “Testing the zero-effect null” in the main
package vignette; the treatment here is the deeper, inference-focused
one.)

### Why this is stronger than a confidence interval on an unregularized estimator

For an unregularized estimator — OLS, or plain
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md)
— the coefficient estimate $`\hat\beta`$ is a continuous random
variable, so $`\hat\beta = 0`$*exactly* is a **measure-zero event**: it
essentially never happens, and would carry no information about the
truth if it did. The only zero-effect inference such an estimator offers
is “the confidence interval contains $`0`$,” i.e. *fail to reject* —
never *accept*. Selection consistency changes the picture qualitatively.
Under Theorem 6.2, $`\hat\beta_j = 0`$ exactly with probability
$`\to 1`$ when the truth is zero, and $`\hat\beta_j \neq 0`$ with
probability $`\to 1`$ when the truth is nonzero, so observing
$`\hat\beta_j = 0`$*is* informative. The paper states the consequence
directly: if $`\tau_{\text{ATT}}(g,t) = 0`$ then
$`\lim_{N\to\infty} P(\hat\tau_{\text{ATT}}(g,t) = 0) = 1`$, which it
calls “stronger than convergence in probability to $`0`$” (Faletto 2025,
around the statement of Theorem 6.2).

It is worth being precise about *which* part of the theory delivers the
test, because the paper is careful here. A conventional Wald statistic
**cannot** test the zero null in this setting: when the relevant effect
is truly zero the statistic is not asymptotically normal (it converges
in probability to $`0`$ even after scaling by $`\sqrt{N}`$), so there is
no valid Wald test of $`H_0: \tau = 0`$ to construct. That is *precisely
why* the right tool is selection consistency and the `selected` flag,
rather than a Wald statistic on the point estimate. The complementary
half of the theory backs the effects the estimator *retains*: when an
effect is truly nonzero, FETWFE is $`\sqrt{N}`$-consistent and
asymptotically Gaussian on the selected support (the oracle property,
Theorem 6.3), with the practical Wald confidence intervals and p-values
of Theorem 6.4 — the same machinery documented in Sections 1–2. So
`selected = TRUE` cohorts come with the usual valid CIs and finite
p-values, while `selected = FALSE` cohorts carry the stronger,
selection-based zero-effect conclusion.

### The `q < 1` caveat

This interpretation rests on selection consistency, which the paper
proves **only for the bridge exponent $`q \in (0, 1)`$**. The package
default $`q = 0.5`$ satisfies it. For $`q = 1`$ (the lasso) selection
consistency of FETWFE is *not* proven, so the implicit-test framing does
not apply; for $`q = 2`$ (ridge) the solver never zeros coefficients at
all, so `selected` is trivially `TRUE` everywhere and conveys nothing.
Read the `selected` flag as a zero-effect test only under the default
(or any $`q < 1`$) bridge penalty.

### A controlled-truth worked example

The cleanest way to see the implicit test is on simulated data where the
truth is known. We draw coefficients in which the *middle* treated
cohort has a true average treatment effect of exactly zero while its two
neighbors are nonzero, then check whether FETWFE selects out exactly the
zero cohort.
[`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md)
reports the known per-cohort truth for comparison.

``` r

# Known truth: cohort 3 (the middle treated cohort) has a true effect of
# exactly zero; cohorts 2 and 4 are nonzero.
coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.2, eff_size = 2, seed = 3)
getTes(coefs)$actual_cohort_tes
#> [1] -1.600000  0.000000  1.333333

dat <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 3)
res <- fetwfeWithSimulatedData(dat, verbose = FALSE)

# The implicit test: the `selected` and `p_value` columns of catt_df.
res$catt_df[, c("cohort", "estimate", "se", "p_value", "selected")]
#>   cohort  estimate         se p_value selected
#> 1      2 -1.425100 0.12554417       0     TRUE
#> 2      3  0.000000 0.00000000      NA    FALSE
#> 3      4  1.426826 0.09635078       0     TRUE
```

The true effects are $`(-1.6, 0, 1.333)`$: cohort 3’s effect is exactly
zero. FETWFE selects out precisely that cohort — cohort 3 returns
`selected = FALSE` with an estimate of exactly $`0`$ and `p_value = NA`
— while it recovers the two truly-nonzero cohorts with `selected = TRUE`
and finite standard errors and p-values. Reading the `selected` column
through the test mapping above: cohort 3 is *“fail to reject (and
conclude zero),”* cohorts 2 and 4 are *“reject.”* The estimator
performed the zero-effect test automatically, and its verdict matches
the known truth.

### Scope: this is the *post-treatment* zero-effect test

The package currently estimates only post-treatment treatment effects
(event times $`e \geq 0`$), so the demonstration above is specifically
the **post-treatment** zero-effect test, surfaced through `catt_df`. The
same selection-consistency logic would support two further applications
once the corresponding extended regression specifications exist: (a)
selecting out pre-treatment placebo coefficients would give a
**parallel-trends test**, and (b) selecting out heterogeneous-trend
augmenting parameters would give a **homogeneous-trends test**. The
package does not estimate those terms today, so this vignette makes no
claim to perform either test; they are noted only as the natural
extensions the same `selected`-as-implicit-test reasoning will cover
when the specifications land.

## 8. High-dimensional inference ($`p \geq NT`$): `debiasedATT()`

Everything above assumes a low-dimensional design ($`p < NT`$), where
the analytic Gram inverse exists and the fit’s standard errors are the
asymptotically valid ones from Section 1. When the design is
**high-dimensional** — the number of design columns $`p`$ (cohort and
time fixed effects, covariates, and all their treatment interactions)
meets or exceeds the number of observations $`NT`$ — two things change:

1.  The fused fit is a *regularized* (bridge) estimate, and its
    pointwise standard errors (and the simultaneous bands of Section 6)
    are **post-selection**: valid conditional on the selected support,
    not *uniformly* valid. They can under-cover, because the bridge’s
    regularization biases the point estimate (it shrinks and prunes
    effects) and the post-selection standard error does not account for
    the selection.

2.  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
    returns a **uniformly valid** confidence interval for the overall
    ATT by *desparsifying* the fit (the nodewise relaxed-inverse
    construction of the high-dimensional FETWFE theory in Faletto 2025),
    which removes the regularization bias. This is “one estimator, two
    regimes”: the call and the returned object are identical to the
    $`p < NT`$ case; only the internal nuisance changes.

The $`p \geq NT`$ regime arises with many covariates, short panels, or
many cohorts.
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)
can generate such panels.

### A worked example

``` r

# d = 12 covariates, N = 40 units, T = 5 periods, G = 3 cohorts.
coefs_hd <- genCoefs(G = 3, T = 5, d = 12, density = 0.2, eff_size = 2, seed = 2)
sim_hd <- simulateData(
  coefs_hd, N = 40, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 2
)
c(p = sim_hd$p, NT = sim_hd$N * sim_hd$T) # p >= NT: high-dimensional
#>   p  NT 
#> 220 200
```

``` r

fit_hd <- fetwfeWithSimulatedData(sim_hd, q = 0.5)
# The fused (post-selection, pointwise) overall ATT:
round(c(att = fit_hd$att_hat, se = fit_hd$att_se), 3)
#>    att     se 
#> -2.416  0.186
```

That fused estimate is post-selection. For uniformly valid inference,
pass the fit to
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md),
selecting the nodewise penalty by cross-validation (`lambda_c = "cv"`) —
the configuration whose overall-ATT coverage was validated at the anchor
below. (The default fixed `lambda_c = 1.0` is a starting point, not the
validated penalty.)

``` r

db <- debiasedATT(fit_hd, lambda_c = "cv")
db
#> Debiased overall ATT (high-dimensional regime)
#>   Estimate:   -2.3684
#>   Std. Error: 0.2507
#>   95% CI: [-2.8597, -1.8771]
#> High-dimensional (p >= NT) desparsified interval [EXPERIMENTAL:
#>   overall-ATT coverage was validated near-nominally at a p >= NT anchor
#>   with the CV-selected penalty (Faletto 2025), but
#>   coverage at other penalties is not separately validated -- inspect the diagnostics below]
#>   nodewise penalty: lambda_c = 0.5 (CV-selected); lambda_node in [0.1836, 0.1836]
#>   nodewise directions: 1/1 converged, 1/1 KKT-feasible (worst feasibility/lambda_node = 1)
```

The print surfaces the high-dimensional diagnostics worth inspecting:
the resolved nodewise penalty (`lambda_c` / `lambda_node`) and, per
debiasing direction, whether it **converged** and met its **KKT
feasibility** certificate. A direction that fails to converge or stay
feasible flags a potentially unreliable standard error — raise
`riesz_max_iter` and/or `lambda_c` and re-inspect.

### Why debiased, and how well it covers

In any single draw the fused and debiased point estimates can be close,
so the value of
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
is not visible from one fit — it is a **repeated-sampling** guarantee. A
coverage study at a $`p > NT`$ anchor ([issue
\#88](https://github.com/gregfaletto/fetwfe/issues/88)) found:

| Overall-ATT interval | Coverage (nominal 95%) |
|----|----|
| Debiased ([`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md) with `lambda_c = "cv"`) | **≈ 93%** |
| Fused plug-in (pointwise) | ≈ 77% |

The fused interval under-covers because of the post-selection
non-uniformity; the debiased interval is near-nominal. This is validated
at the studied anchor (with the feasibility-appropriate penalty); treat
the $`p \geq NT`$ path as **experimental** outside that regime — which
is why the diagnostics above matter.

### An asymptotic alternative: the wild bootstrap

`method = "bootstrap"` offers an alternative, asymptotically-valid
reference distribution for the overall-ATT interval: a studentized score
/ influence-function **wild cluster bootstrap** that re-signs the
per-unit influence summands (no refit per replicate), leaving the point
estimate and the reported `se` unchanged and changing only the critical
value. It is **floored at the Gaussian quantile**, so the interval is
never narrower than the analytic Wald interval.

**It is not a few-clusters remedy.** The analytic sandwich SE is
downward-biased with few clusters, but this bootstrap does not fix that:
as an *unrestricted, no-refit* score bootstrap it does not reproduce the
tail inflation the *restricted* wild-cluster bootstrap-t uses for
few-cluster coverage, and under heterogeneous cluster influence its
critical value (before the floor) actually falls *below* the Gaussian —
so the floor reduces it to the analytic interval in exactly the
small-`N` regime. Its only effect is to *widen* the interval when
cluster influence is near-homogeneous. For genuine few-clusters
inference a restricted bootstrap-t or a CR2-type analytic correction is
required (future work).

``` r

# `fit` must be a single-sample fetwfe() fit -- not an `indep_counts` fit, which
# the worked example above (via fetwfeWithSimulatedData) happens to be.
debiasedATT(fit, method = "bootstrap", multiplier = "webb", seed = 1)
```

The `multiplier` (`"webb"`, the default; `"rademacher"`; or `"mammen"`)
only affects the near-homogeneous case where the bootstrap widens above
the floor — in the small-`N` / heterogeneous regime the interval is the
analytic one regardless of the weight. The bootstrap is defined only for
**single-sample** fits (real observational panels are single-sample),
not `indep_counts` two-sample fits.

### Simultaneous bands at $`p \geq NT`$

[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
also switches to the desparsified construction in the high-dimensional
regime (via `method = "bootstrap"`; the analytic Gram-inverse band does
not exist when $`p \geq NT`$):

``` r

sc_hd <- simultaneousCIs(
  fit_hd,
  family = "cohort",
  method = "bootstrap",
  B = 500,
  seed = 1
)
sc_hd$regime
#> [1] "high-dimensional"
```

One scope caveat: \#88 validated the **scalar overall-ATT**
([`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md));
the **family-wise band** coverage at $`p \geq NT`$ is not itself
simulation-validated, so the band is the more experimental of the two.

### Controlled high-dimensional truth

The default `density` placement spreads signal across all coordinates,
which when $`p`$ greatly exceeds the number of treatment-effect
parameters rarely lands signal on the treatment-effect block. For a
controlled, sparse, non-degenerate high-dimensional data-generating
process,
[`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
offers a **targeted-sparsity** mode (`n_signal_cohorts` or
`treat_base_levels`) that places a heterogeneous treatment-effect signal
directly on the per-cohort base levels — the DGP the \#88 coverage study
uses. See
[`?genCoefs`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md).

## References

Bertrand, M., Duflo, E., & Mullainathan, S. (2004). “How much should we
trust differences-in-differences estimates?” *Quarterly Journal of
Economics* 119(1), 249–275.

Faletto, G. (2025). “Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions.” *arXiv preprint*
[arXiv:2312.05985](https://arxiv.org/abs/2312.05985).

Kock, A. B. (2013). “Oracle Efficient Variable Selection in Random and
Fixed Effects Panel Data Models.” *Econometric Theory* 29(1), 115–152.

Liang, K.-Y., & Zeger, S. L. (1986). “Longitudinal data analysis using
generalized linear models.” *Biometrika* 73(1), 13–22.
