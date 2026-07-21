# Parametric simultaneous (1 - alpha) confidence intervals over a family of treatment effects

Computes simultaneous (family-wise) confidence intervals for a user-
specified family of treatment effects from a fitted FETWFE / ETWFE /
BETWFE / twfeCovs object. The simultaneous critical value
`c_{1 - alpha}` is the `(1 - alpha)` quantile of `max_k |Z_k|` where `Z`
follows a multivariate normal with correlation matrix `cov2cor(Sigma)`;
it is computed deterministically via
[`mvtnorm::qmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/qmvnorm.html).
Under Faletto (2025) Theorem (c') tight Gaussianity and Assumption
(Psi-IF), the family of psi-linear effects is asymptotically
multivariate normal with a covariance that is estimable from the
package's existing variance machinery; under the paper's fixed-dim
framing no high-dimensional correction is needed.

The pointwise critical value `qnorm(1 - alpha/2)` (per-effect coverage)
and the Bonferroni-conservative critical value `qnorm(1 - alpha/(2K))`
(family- wise coverage with no correlation assumption) are returned for
side-by-side comparison; the simultaneous critical value is always
between them when the effects are positively correlated (as is typical
in difference-in- differences, where effects share the
regression-coefficient variance piece).

## Usage

``` r
simultaneousCIs(
  result,
  family = c("event_study", "cohort", "all_post_treatment", "custom"),
  alpha = NULL,
  contrasts = NULL,
  method = c("analytic", "bootstrap"),
  B = 1000L,
  seed = NULL,
  multiplier = c("rademacher", "mammen", "webb"),
  lambda_c = 1,
  riesz_max_iter = 5000L,
  riesz_tol = 1e-09,
  cv_time_budget = Inf
)
```

## Arguments

- result:

  A fitted object of class `"fetwfe"`, `"etwfe"`, `"betwfe"`, or
  `"twfeCovs"`.

- family:

  Character; one of `"event_study"`, `"cohort"`, `"all_post_treatment"`,
  or `"custom"`. See Details for each family's resolution.

- alpha:

  Numeric in `(0, 1)`; significance level. Default `NULL`, which
  inherits the fit's own `alpha` (the level its `catt_df` / displayed
  simultaneous bands were built at), so `simultaneousCIs(fit)`
  reproduces the fit's shown simultaneous bounds. Pass an explicit value
  to override.

- contrasts:

  For `family = "custom"`, a `K x num_treats` matrix whose rows give the
  `K` linear combinations of the underlying per-`(g, t)`
  treatment-effect vector (the `multcomp::glht()` convention;
  `num_treats` is the number of estimated effects, equal to `G` for
  `twfeCovs`). Ignored for the other families. Note that the `"custom"`
  family omits the cohort-probability variance term (`Sigma_2 = 0`), so
  a custom contrast that pools across cohorts in a probability-weighted
  way is anti-conservative (its band can under-cover); use
  `family = "cohort"` for cohort-pooled effects.

- method:

  Character; how the simultaneous critical value is computed.
  `"analytic"` (default) uses the exact multivariate-normal sup-t
  quantile via
  [`mvtnorm::qmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/qmvnorm.html).
  `"bootstrap"` uses the multiplier bootstrap of Chernozhukov,
  Chetverikov & Kato (2013): it perturbs the per-unit influence
  functions and reads the sup-t critical value off the bootstrap
  distribution. The two are asymptotically equivalent; the bootstrap
  scales better to large effect families and is
  heteroskedasticity/cluster-robust by construction (it always reports
  cluster-robust per-unit standard errors, regardless of the fit's
  `se_type`). `method = "bootstrap"` supports all four families. For
  `family = "event_study"` the bootstrap additionally perturbs the
  per-unit cohort-probability (propensity) influence function with its
  own independent multiplier stream – a two-channel bootstrap matching
  the analytic variance `Sigma = Sigma_1 + Sigma_2` (the event-study
  family is the one whose `Sigma_2` is non-zero, because each event-time
  effect pools across cohorts weighted by the estimated cohort
  probabilities). When the (full) design is **high-dimensional
  (`p >= NT`)** – where the analytic Gram inverse need not exist – the
  bootstrap uses the full-design **desparsified** construction of
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  (per-effect nodewise directions) generalized to the family. This
  desparsified `p >= NT` band path
  ([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  fits only) generalizes the
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  construction — whose overall-ATT coverage is validated near-nominally
  in simulation (Theorem `debiased.highdim.thm`, Faletto 2025) — to a
  family of functionals. Its family-wise coverage is likewise validated
  near-nominally (Theorem `debiased.highdim.joint.thm`): ~0.92 over the
  **event-study** family (`K = 7`) at the `p >= NT` anchor of Faletto
  (2025), the band sibling of the scalar result (~0.94); the cohort and
  all-post-treatment families are covered by the same theorem but were
  not separately simulated. As with the scalar, validity rests on the
  sparsity / restricted-eigenvalue primitives (assumed, not proved), so
  inspect the returned `feasibility` / `converged` diagnostics. A
  non-[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  `p >= NT` fit (e.g.
  [`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md))
  has no desparsified band, so it instead falls back to the fixed-`p`
  selected-support band and emits a
  [`warning()`](https://rdrr.io/r/base/warning.html): a \#308 coverage
  study shows that band substantially **under-covers** in the `p >= NT`
  regime – even when the selected support is low-dimensional – because
  the bridge shrinks the band center toward zero and the post-selection
  standard error understates variability, so the band is returned but
  should be treated as unreliable (use a
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
  fit for valid high-dimensional bands). The desparsified path covers
  all four families (a high-dimensional `family = "event_study"` fit
  additionally carries the propensity channel `F_pi`). In the
  high-dimensional regime the band is centered on the **debiased**
  estimate (the high-dimensional FETWFE theory correction, equal to
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)'s
  point estimate for the matching contrast), not the post-selection
  bridge estimate; fixed-p bands center on the (unbiased) bridge
  estimate as before. Because this desparsified bootstrap band is built
  from the empirical per-unit influence function – needing neither GLS
  whitening nor `sig_eps_sq` – it also accepts a high-dimensional
  **`gls = FALSE`** fetwfe fit (`calc_ses = FALSE`), the band analog of
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)'s
  Omega-free SE (#307, \#313). `method = "analytic"` still requires
  valid analytic SEs (a `q < 1`, `gls = TRUE`, rank-satisfied fit). **At
  `p >= NT`, `method = "bootstrap"` is the valid simultaneous band** —
  the desparsified Theorem `debiased.highdim.joint.thm` band above. A
  `method = "analytic"` call there instead returns the analytic band on
  the *selected support*: a post-selection band that under-covers (it is
  not the uniformly-valid desparsified band), and a `gls = FALSE` fit
  has no analytic SEs at all — so pass `method = "bootstrap"` for
  high-dimensional bands.

- B:

  Integer; number of multiplier-bootstrap replicates
  (`method = "bootstrap"` only). Default `1000`.

- seed:

  Optional integer; if supplied, the bootstrap draws are reproducible
  (the ambient random-number stream is saved and restored around them).
  If `NULL` (default), the draws come from the ambient generator, so
  results vary run to run. Ignored when `method = "analytic"`.

- multiplier:

  Character; the multiplier-weight distribution for the bootstrap:
  `"rademacher"` (default, `+/-1`), `"mammen"` (the Mammen 1993
  two-point distribution), or `"webb"` (the Webb 2013 six-point
  distribution, preferred with very few clusters). Ignored when
  `method = "analytic"`.

- lambda_c, riesz_max_iter, riesz_tol:

  Controls for the high-dimensional (`p >= NT`) bootstrap, where each
  effect's debiasing direction is a nodewise (desparsified-lasso)
  `riesz_lasso()` solve. `lambda_c` is the leading constant of the
  penalty `lambda_node = lambda_c * max(|a|) * sqrt(log p / N)`: either
  a single positive number (the default `1.0` = theory scale) or the
  string `"cv"`, which selects the constant **per fit** by
  cross-validation (the same CV
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  uses, on the overall-ATT direction, so one constant serves the point
  estimate and every band effect; \#295). **Smaller fixed values can
  leave directions infeasible (a warning fires and those bands are
  unreliable).** Ignored when `p < NT` or `method = "analytic"`. The
  default stays the fixed `1.0`; the family-wise coverage validated in
  Faletto (2025) used the `"cv"` selector (`lambda_c` ~ 0.54,
  feasibility 1.0), the same constant that serves the scalar
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  and every band effect.

- cv_time_budget:

  Numeric; a wall-clock backstop (in seconds) for the `lambda_c = "cv"`
  cross-validation (the same \#384 backstop
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  offers). `Inf` (the default) leaves selection fully deterministic /
  reproducible; a finite value stops the CV early and falls back to the
  theory scale (`lambda_c = 1.0`) with a warning, so an adversarial
  high-dimensional draw cannot spin indefinitely (#384). Whether a
  finite budget fires depends on machine speed, so the result can too.
  Ignored unless `lambda_c = "cv"` and `p >= NT`.

## Value

An object of S3 class `"simultaneous_cis"`: a list with

- ci:

  A data frame with columns `effect`, `estimate`, `simultaneous_ci_low`,
  `simultaneous_ci_high`, `pointwise_ci_low`, `pointwise_ci_high` (one
  row per effect in the family).

- adjusted_p_values:

  Numeric vector of length `K`: the single-step max-T
  multiplicity-adjusted (family-wise) p-value for each effect, the exact
  dual of the simultaneous band (a coefficient lies outside the
  `(1 - alpha)` band iff its adjusted p-value is `< alpha`). Computed
  via
  [`mvtnorm::pmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/pmvnorm.html)
  over the same correlation matrix the band uses (or, under
  `se_type = "conservative"`, the Bonferroni adjustment
  `min(1, K * pointwise_p)`). `NA` for degenerate (zero-variance)
  effects. (#200)

- critical_value:

  The simultaneous critical value `c_{1 - alpha}` (or, when the fit used
  `se_type = "conservative"`, the Bonferroni critical value
  `qnorm(1 - alpha/(2K))` – see Details).

- pointwise_critical_value:

  `qnorm(1 - alpha/2)`, for reference.

- bonferroni_critical_value:

  `qnorm(1 - alpha/(2K))`, for reference.

- family:

  The requested family (character).

- alpha:

  The significance level used.

- K:

  The number of effects in the family (integer).

For `method = "bootstrap"` the list additionally carries `method`, `B`,
`seed`, `multiplier`, and `regime` (`"fixed-p"` or
`"high-dimensional"`); the `critical_value` is the bootstrap sup-t
quantile and the standard errors backing `ci` are the cluster-robust
per-unit ones (so for a non-`cluster` fit the bootstrap band may differ
from the analytic homoskedastic band). A `"high-dimensional"` fit
additionally carries per-effect `feasibility`, `converged`, and
`lambda_node` (the nodewise-direction diagnostics), plus `lambda_c` (the
leading constant used) and `lambda_c_selection` (`"fixed"` or `"cv"`).
When `lambda_c = "cv"` it also carries `lambda_cv`, the cross-validation
diagnostics (the grid, the per-grid feasibility flags and CV losses, and
whether the theory-scale fallback fired), matching
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)'s
`lambda_cv`.

## Details

**Family resolution and `K`.** `"event_study"` resolves to one effect
per post-treatment event time `e = 0, ..., T - 2` (`K = T - 1`);
`"cohort"` to one effect per treated cohort (`K = G`);
`"all_post_treatment"` to one effect per `(g, t)` cell
(`K = num_treats`); `"custom"` to the `K = nrow(contrasts)`
user-supplied contrasts.

**Joint covariance.** The `K x K` covariance `Sigma = Sigma_1 + Sigma_2`
is reconstructed at call time from the fit's stored slots (design
matrix, selected support, `theta_hat` / `beta_hat`,
`cohort_probs_overall`, `sig_eps_sq`). `Sigma_1` is the
regression-coefficient piece and `Sigma_2` the cohort-probability piece,
generalizing the package's per-point variance machinery (the same
machinery
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
uses). By construction `sqrt(diag(Sigma))` equals the package's existing
per-point standard errors for the corresponding effects. The `Sigma`
blocks are not persisted on the fit; re-derivation is sub-second.

**Degenerate (zero-variance) effects.** An effect whose entire
contribution to the selected support is zeroed by the bridge penalty –
or, in scattered-cohort panels, an event time with an empty valid-cohort
set – has a standard error of exactly 0 by construction, so its
simultaneous and pointwise CIs collapse to a point at the estimate and
it is excluded from the joint correlation matrix (it adds no family-wise
risk; the critical value is computed over the non-degenerate
sub-family). This `se = 0` convention is the simultaneous-CI analog of
the `NA` standard error
[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
reports for the same structurally-degenerate event times; both assign
the effect an estimate of 0.

**Paper grounding.** Theorem (c') tight Gaussianity (Faletto 2025,
`paper_arxiv.tex:1233`) guarantees the joint asymptotic normality;
Assumption (Psi-IF) (assumption equation `paper_arxiv.tex:2013`;
in-prose discussion at paper line 1268) is the influence-function
condition the package's default cohort-sample-proportions estimator
satisfies; the fixed-dim framing follows the paper's AE point 1(d).

**Conservative fallback.** When the fit was made with
`se_type = "conservative"`, the function falls back to
Bonferroni-corrected pointwise CIs (the Cauchy-Schwarz upper bound used
for the conservative scalar SE does not generalize to a `K x K`
covariance matrix) and emits a brief
[`message()`](https://rdrr.io/r/base/message.html). The
`$critical_value` field is set to the Bonferroni value in this branch.

**Numerical integration.** The critical value is computed via
`mvtnorm::qmvnorm(..., algorithm = mvtnorm::GenzBretz())` (mvtnorm's
default quasi-Monte Carlo integrator; sub-second through `K` up to about
100, so no `K` cap is of practical concern for FETWFE families).
`mvtnorm` is an `Imports` dependency (as of version 1.16.0, when
simultaneous bands became the default reported confidence interval; see
the `ci_type` argument of
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)).
The function uses it only when `K > 1` and `se_type != "conservative"`
(the `K = 1` and conservative paths bypass the dependency), and retains
a defensive [`stop()`](https://rdrr.io/r/base/stop.html) with an
actionable message if it is somehow unavailable (e.g., a corrupted
install).

**Determinism contract.** The function is deterministic in its inputs:
the same fit plus the same `family`, `alpha`, and `contrasts` always
produces the same critical value across calls. This is achieved by
wrapping the internal
[`mvtnorm::qmvnorm()`](https://rdrr.io/pkg/mvtnorm/man/qmvnorm.html)
call with a save/restore of the caller's `.Random.seed` and a fixed
internal `set.seed(1L)` immediately before the call. The function does
NOT mutate the caller's `.Random.seed` (the save/restore via
[`on.exit()`](https://rdrr.io/r/base/on.exit.html) leaves the caller's
RNG state identical pre- and post-call), matching the convention adopted
by `R/fetwfe_core.R::getBetaCV()` in PR \#181 / v1.13.5. Users do not
need to call [`set.seed()`](https://rdrr.io/r/base/Random.html) before
`simultaneousCIs()` to get reproducible results, and downstream
RNG-using code observes no perturbation.

## References

Faletto, G. (2025). Fused Extended Two-Way Fixed Effects for
Difference-in- Differences with Staggered Adoptions. arXiv:2312.05985.

Hothorn, T., Bretz, F., & Westfall, P. (2008). Simultaneous Inference in
General Parametric Models. *Biometrical Journal* 50(3), 346-363.

## See also

[`eventStudy()`](https://gregfaletto.github.io/fetwfePackage/reference/eventStudy.md)
for the per-point event-study estimates and Wald intervals that
`family = "event_study"` provides simultaneous bands over.

## Examples

``` r
# \donttest{
  coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
  sim <- simulateData(coefs, N = 120, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 123)
  fit <- fetwfeWithSimulatedData(sim)
  sci <- simultaneousCIs(fit, family = "event_study", alpha = 0.05)
  print(sci)
#> Parametric simultaneous (1 - alpha) confidence intervals
#> Family: event_study, K = 5, alpha = 0.05
#> Simultaneous critical value: 2.4454 (pointwise 1.9600, Bonferroni 2.5758)
#> 
#>  effect estimate simultaneous_ci_low simultaneous_ci_high pointwise_ci_low
#>      e0 0.000000           0.0000000             0.000000         0.000000
#>      e1 1.354096           0.9933308             1.714862         1.064944
#>      e2 3.242404           2.8548785             3.629930         2.931804
#>      e3 3.700280           2.9969740             4.403586         3.136583
#>      e4 5.751528           5.3111327             6.191922         5.398553
#>  pointwise_ci_high
#>           0.000000
#>           1.643248
#>           3.553005
#>           4.263977
#>           6.104502
# }
```
