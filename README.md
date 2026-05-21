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

The primary function in the `{fetwfe}` is `fetwfe()`, which implements fused extended two-way fixed effects. Here's some example code applying `fetwfe()` to the `castle` data set from the `bacondecomp` package (the same example used in the package vignette):

```R
library(fetwfe)
library(bacondecomp)

data(castle)

# Response: the log homicide rate. Treatment: `cdl` records the share of
# the year the castle-doctrine law was in effect, so `cdl > 0` gives the
# absorbing 0/1 treatment indicator `fetwfe()` requires.
castle$l_homicide <- log(castle$homicide)
castle$treated <- as.integer(castle$cdl > 0)

res <- fetwfe(
    pdata = castle,
    time_var = "year",
    unit_var = "state",
    treatment = "treated",
    response = "l_homicide",
    verbose = TRUE)

summary(res)
```

For vignettes and full documentation, check out the [page for the `{fetwfe}` package on CRAN](https://CRAN.R-project.org/package=fetwfe).

## References
- Faletto, G (2025). *Fused Extended Two-Way Fixed Effects for Difference-in-Differences with Staggered Adoptions*. [arXiv preprint arXiv:2312.05985](https://arxiv.org/abs/2312.05985).

