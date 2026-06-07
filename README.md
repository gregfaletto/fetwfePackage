# Fused Extended Two-Way Fixed Effects

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

## References
- Faletto, G (2025). *Fused Extended Two-Way Fixed Effects for Difference-in-Differences with Staggered Adoptions*. [arXiv preprint arXiv:2312.05985](https://arxiv.org/abs/2312.05985).

