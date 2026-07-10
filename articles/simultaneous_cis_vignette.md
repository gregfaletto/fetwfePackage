# Simultaneous confidence bands and the zero-effect test

``` r

library(fetwfe)
```

This vignette covers the package’s **family-wise inference** tools for
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
[`betwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/betwfe.md),
[`etwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/etwfe.md),
and
[`twfeCovs()`](https://gregfaletto.github.io/fetwfePackage/reference/twfeCovs.md):
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md),
which returns confidence bands that hold *jointly* across a family of
effects (not merely one effect at a time), and the
**selection-consistency zero-effect test** — the sense in which
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)’s
`selected` flag is an implicit test of whether a cohort’s effect is
truly zero.

It is a companion to
[`vignette("inference_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/inference_vignette.md)
(how the package’s standard errors are computed) and
[`vignette("high_dimensional_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/high_dimensional_vignette.md)
(debiased inference when $`p \geq NT`$). We reuse that first vignette’s
small running example:

``` r

set.seed(2026)
sim_coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
sim_data <- simulateData(
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
```

## 1. Simultaneous confidence intervals via `simultaneousCIs()`

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
setup above:

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
regime. It is validated in simulation
([`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
fits only): the scalar overall-ATT coverage (Theorem
`debiased.highdim.thm`) and the family-wise *band* coverage (Theorem
`debiased.highdim.joint.thm`, ≈0.92 over the event-study family at the
$`p > NT`$ anchor of Faletto 2025) are both near-nominal, so the
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
default stays the fixed `1.0`; the overall-ATT coverage validated in
Faletto (2025) is robust across the penalty values studied (a CV
selector below the theory scale and a fixed constant above it,
bracketing the default `1.0`). The simulator can generate the `p > NT`
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
fit for valid high-dimensional simultaneous bands. See
[`vignette("high_dimensional_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/high_dimensional_vignette.md)
for
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
and the accompanying high-dimensional diagnostics.

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

## 2. Selection consistency as an implicit zero-effect test

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
of Theorem 6.4 — the same machinery documented in
[`vignette("inference_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/inference_vignette.md).
So `selected = TRUE` cohorts come with the usual valid CIs and finite
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

## References

Faletto, G. (2025). “Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions.” *arXiv preprint*
[arXiv:2312.05985](https://arxiv.org/abs/2312.05985).
