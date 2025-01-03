# fetwfePackage
 
To install the `fetwfe` package, simply use

```
# install.packages("remotes")  # if needed
remotes::install_github("gregfaletto/yourPackageName")
library(fetwfe)
```

The `fetwfe` package for now contains one function, `fetwfe()`, which implements fused extended two-way fixed effects. See [Faletto (2024)](https://arxiv.org/abs/2312.05985) for details on the methodology. Here's some example code that implements the data application from the paper:

```
library(fetwfe)
library(bacondecomp)

set.seed(23451)

data(divorce)

res <- fetwfe(
    pdata=divorce[divorce$sex == 2, ],
    time_var="year",
    unit_var="st",
    treatment="changed",
    covs=c("murderrate", "lnpersinc", "afdcrolls"),
    response="suiciderate_elast_jag",
    q=0.5,
    verbose=TRUE)
```
