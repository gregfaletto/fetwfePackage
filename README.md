# Fused Extended Two-Way Fixed Effects

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/fetwfe)](https://CRAN.R-project.org/package=fetwfe)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/fetwfe)](https://CRAN.R-project.org/package=fetwfe)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The `{fetwfe}` package implements *fused extended two-way fixed effects* (FETWFE), a methodology for estimating treatment effects in difference-in-differences with staggered adoptions.

* For a brief introduction to the methodology, as well as background on difference-in-differences with staggered adoptions and motivation for FETWFE, see this [blog post](https://gregoryfaletto.com/2023/12/13/new-paper-fused-extended-two-way-fixed-effects-for-difference-in-differences-with-staggered-adoptions/).
* For more detailed slides on the methodology (but less detailed than the paper), see [here](https://gregoryfaletto.com/2024/02/11/presentation-on-fused-extended-two-way-fixed-effects/).
* Check out the most recent draft of the full paper [here](https://arxiv.org/abs/2312.05985).
* This [blog post](https://gregoryfaletto.com/2025/01/03/new-r-fetwfe-package-implementing-fused-extended-two-way-fixed-effects/) explains a little more about what the package does under the hood, if you're interested.

# fetwfePackage
 
To install the `{fetwfe}` package, simply use

```R
install.packages("fetwfe")
```

You can also install the latest development version by using
```R
# install.packages("remotes")  # if needed
remotes::install_github("gregfaletto/fetwfePackage")
```

The primary function in the `{fetwfe}` is `fetwfe()`, which implements fused extended two-way fixed effects. Here's some example code applying `fetwfe()` to the `divorce` data set from the `bacondecomp` package -- the unilateral ("no-fault") divorce / female-suicide panel of Stevenson and Wolfers (2006):

```R
library(fetwfe)
library(bacondecomp)

data(divorce)

# Restrict to the female subset (`sex == 2`). `changed` is already the absorbing
# 0/1 divorce-reform indicator, and the response is the elasticity-scaled female
# suicide rate. This reproduces the empirical application in Faletto (2025,
# Sec. 8.2): the three covariates are passed as controls (one, `murderrate`, is
# auto-dropped because it is missing in 1964, leaving two), and the noise
# variances are supplied (precomputed by REML) to keep the call fast.
divorce_f <- divorce[divorce$sex == 2, ]

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

# Overall ATT is about -6% on the elasticity-scaled female suicide rate, with a
# 95% CI that excludes zero; FETWFE retains heterogeneous cohort effects here.
summary(res)
```

For vignettes and full documentation, check out the [page for the `{fetwfe}` package on CRAN](https://CRAN.R-project.org/package=fetwfe).

## Penalty and fusion structures

A core differentiator of FETWFE is the *geometry* of its fusion penalty, chosen via the `fusion_structure` argument of `fetwfe()`:

- **`"cohort"`** (the default) fuses treatment effects within and between cohorts — the two-way fusion structure of Faletto (2025).
- **`"event_study"`** instead fuses effects at the same time since treatment (event time `e = t - g`) across cohorts — the package realization of the paper's event-study-penalty theory.

For full control, the `fusion_matrix` argument accepts any invertible `num_treats x num_treats` differencing matrix `D_N`, encoding an arbitrary fusion structure beyond the two built-ins. The same `fusion_structure` option is available in `genCoefs()` for simulation studies.

For guidance on which to use, see the vignette *"Choosing a fusion structure: cohort vs. event-study penalties"* — `vignette("fusion_structure_vignette", package = "fetwfe")`, also on the [CRAN page](https://CRAN.R-project.org/package=fetwfe).

## References
- Faletto, G (2025). *Fused Extended Two-Way Fixed Effects for Difference-in-Differences with Staggered Adoptions*. [arXiv preprint arXiv:2312.05985](https://arxiv.org/abs/2312.05985).

