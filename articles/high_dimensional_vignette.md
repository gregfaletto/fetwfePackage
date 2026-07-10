# High-dimensional inference: debiasedATT() when p \>= NT

``` r

library(fetwfe)
```

This vignette covers **high-dimensional inference** in `fetwfe`: what
changes when the design has as many or more columns than observations
($`p \geq NT`$), and how
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
restores a uniformly valid confidence interval for the overall ATT.

It is a companion to
[`vignette("inference_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/inference_vignette.md)
(how the package’s standard errors are computed in the low-dimensional
regime) and
[`vignette("simultaneous_cis_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/simultaneous_cis_vignette.md)
(family-wise simultaneous bands).

## 1. High-dimensional inference ($`p \geq NT`$): `debiasedATT()`

The standard errors and simultaneous bands documented in the companion
vignettes assume a **low-dimensional** design ($`p < NT`$), where the
analytic Gram inverse exists and the fit’s standard errors are the
asymptotically valid ones of
[`vignette("inference_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/inference_vignette.md).
When the design is **high-dimensional** — the number of design columns
$`p`$ (cohort and time fixed effects, covariates, and all their
treatment interactions) meets or exceeds the number of observations
$`NT`$ — two things change:

1.  The fused fit is a *regularized* (bridge) estimate, and its
    pointwise standard errors (and the simultaneous bands of
    [`vignette("simultaneous_cis_vignette")`](https://gregfaletto.github.io/fetwfePackage/articles/simultaneous_cis_vignette.md))
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
which implements the high-dimensional debiased construction of **Theorem
`debiased.highdim.thm`** (Faletto 2025). Its coverage is robust across
the penalty values studied — a cross-validated selector below the theory
scale and a fixed constant above it bracket the default `1.0` — so the
default is fine; below we let `lambda_c = "cv"` select the constant per
fit as an alternative.

``` r

db <- debiasedATT(fit_hd, lambda_c = "cv")
db
#> Debiased overall ATT (high-dimensional regime)
#>   Estimate:   -2.3684
#>   Std. Error: 0.2507
#>   95% CI: [-2.8597, -1.8771]
#> High-dimensional (p >= NT) desparsified interval
#>   [overall-ATT coverage validated near-nominally in simulation
#>   (Theorem debiased.highdim.thm, Faletto 2025); diagnostics below]
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
| Debiased ([`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md) with `lambda_c = "cv"`) | **≈ 94%** |
| Fused plug-in (pointwise) | ≈ 77% |

The fused interval under-covers because of the post-selection
non-uniformity; the debiased interval is near-nominal — the direct
simulation confirmation of **Theorem `debiased.highdim.thm`**, robust
across the penalty values studied (a CV selector below the theory scale
and a fixed constant above it, bracketing the default `1.0`). The
per-direction `converged` / `feasibility` diagnostics above remain worth
a glance to confirm the nodewise directions are well-behaved on your own
fit.

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

One scope caveat: **Theorem `debiased.highdim.thm`** (validated in
simulation) covers the **scalar overall-ATT**
([`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)),
but the **family-wise band** (Theorem `debiased.highdim.joint.thm`)
coverage at $`p \geq NT`$ is *not* itself simulation-validated — the
paper states it experimental-grade — so treat these high-dimensional
bands as the experimental part of this path.

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

Faletto, G. (2025). “Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions.” *arXiv preprint*
[arXiv:2312.05985](https://arxiv.org/abs/2312.05985).
