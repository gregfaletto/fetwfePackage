# Fused Extended Two-Way Fixed Effects

The `{fetwfe}` package implements *fused extended two-way fixed effects* (FETWFE), a methodology for estimating treatment effects in difference-in-differences with staggered adoptions.

* For a brief introduction to the methodology, as well as background on difference-in-differences with staggered adoptions and motivation for FETWFE, see this [blog post](https://gregoryfaletto.com/2023/12/13/new-paper-fused-extended-two-way-fixed-effects-for-difference-in-differences-with-staggered-adoptions/).
* For more detailed slides on the methodology (but less detailed than the paper), see [here](https://gregoryfaletto.com/2024/02/11/presentation-on-fused-extended-two-way-fixed-effects/).
* Check out the most recent draft of the full paper [here](https://arxiv.org/abs/2312.05985).
* This [blog post](https://gregoryfaletto.com/2025/01/03/new-r-fetwfe-package-implementing-fused-extended-two-way-fixed-effects/) explains a little more about what the package does under the hood, if you're interested.

# fetwfePackage
 
To install the `{fetwfe}` package, simply use

```R
# install.packages("remotes")  # if needed
remotes::install_github("gregfaletto/fetwfePackage")
library(fetwfe)
```

The `{fetwfe}` package contains one function, `fetwfe()`, which implements fused extended two-way fixed effects. Here's some example code that implements the data application from the paper:

```R
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

# Average treatment effect on the treated units (in percentage point
# units)
100 * res$att_hat

# Conservative 95% confidence interval for ATT (in percentage point units)

low_att <- 100 * (res$att_hat - qnorm(1 - 0.05 / 2) * res$att_se)
high_att <- 100 * (res$att_hat + qnorm(1 - 0.05 / 2) * res$att_se)

c(low_att, high_att)

# Cohort average treatment effects and confidence intervals (in percentage
# point units)

catt_df_pct <- res$catt_df
catt_df_pct[["Estimated TE"]] <- 100 * catt_df_pct[["Estimated TE"]]
catt_df_pct[["SE"]] <- 100 * catt_df_pct[["SE"]]
catt_df_pct[["ConfIntLow"]] <- 100 * catt_df_pct[["ConfIntLow"]]
catt_df_pct[["ConfIntHigh"]] <- 100 * catt_df_pct[["ConfIntHigh"]]

catt_df_pct
```

# Documentation

Some documentation for the `fetwfe()` function follows.

## Description
This function implements fused extended two-way fixed effects without sample splitting. It estimates the overall average treatment effect on the treated (ATT) as well as the cohort average treatment effects on the treated units (CATT).

## Usage
```R
fetwfe(
  pdata,
  time_var,
  unit_var,
  treatment,
  covs,
  response,
  indep_counts = NA,
  sig_eps_sq = NA,
  sig_eps_c_sq = NA,
  lambda.max = NA,
  lambda.min = NA,
  nlambda = 100,
  q = 0.5,
  verbose = FALSE,
  alpha = 0.05
)
```

## Arguments
- `pdata`: A data frame containing the panel dataset. Each row in this data frame should represent an observation of a unit at a particular time. The data frame must include columns that meet the specifications described below.
- `time_var`: A character string specifying the name of the column in the data frame that contains the variable for the time period. The values in this column are expected to be integers, such as years. Recommended formats for representing dates include `YYYY`, `YYYYMM`, or `YYYYMMDD`, depending on what is appropriate for your dataset.
- `unit_var`: A character string specifying the name of the column in the data frame that contains a variable identifying each unit. The values in this column are expected to be characters, representing the "name" of each unit.
- `treatment`: A character string specifying the name of the column in the data frame that contains the treatment indicator variable. The values in this column are expected to be integers. Specifically, the value should be `0` if the unit was untreated at the given time period and `1` if the unit was treated. The treatment is expected to be an absorbing state, meaning that once a unit is treated, it remains treated in all subsequent time periods. Units that are treated in the first time period will be automatically removed from the analysis. Ensure that there are untreated units in the final time period to allow proper estimation.
- `covs`: A character vector specifying the names of the columns in the data frame that contain covariates. All of these columns must contain numeric or integer values. If categorical variables are included, they should be encoded, for example, as binary indicators, before being used in this function. At least one covariate must be provided.
- `response`: A character string specifying the name of the column in the data frame that contains the response variable for each unit at each time. This variable must be numeric or integer.
- `indep_counts`: An optional argument. If a large number of units are available, the dataset can be split into two subsets, with one half containing the data provided to this function in `pdata`, and the other half summarized in `indep_counts`. This should be an integer vector representing the number of units in the untreated cohort and each of the treated cohorts. This split allows for exact standard error estimation. The length of `indep_counts` must match the number of cohorts in the data. All values in `indep_counts` must be strictly positive, and their sum must match the total number of units in `pdata`. The default is `NA`, which calculates conservative standard errors.
- `sig_eps_sq`: An optional numeric argument specifying the variance of the row-level independent and identically distributed (IID) noise assumed for each observation. If the variance is known, it is recommended to provide this value. If unknown, it will be estimated using ridge regression as described in Section 26.5.1 of Pesaran (2015). The default is `NA`.
- `sig_eps_c_sq`: An optional numeric argument specifying the variance of the unit-level IID noise. If the variance is known, it is recommended to provide this value. If unknown, it will be estimated using ridge regression as described in Section 26.5.1 of Pesaran (2015). The default is `NA`.
- `lambda.max`: An optional numeric argument specifying the maximum value for the penalty parameter lambda in the grid search used to select the model by the Bayesian Information Criterion (BIC). If not provided, it will be automatically selected. This parameter should be chosen such that the smallest model (corresponding to `lambda.max`) selects very few features, while the largest model (corresponding to `lambda.min`) selects many features. The default is `NA`.
- `lambda.min`: An optional numeric argument specifying the minimum value for the penalty parameter lambda in the grid search. The default is `NA`.
- `nlambda`: An optional integer specifying the total number of penalty parameters to be considered in the grid search. The default is `100`.
- `q`: An optional numeric argument that determines the type of penalty used for fusion regularization. When `q = 1`, the method applies the lasso penalty, while for values of `q` between `0` and `1`, standard errors and confidence intervals can be computed. When `q = 2`, ridge regression is applied. See Faletto (2024) for details. The default is `0.5`.
- `verbose`: A logical argument. If set to `TRUE`, the function will print additional details on its progress during execution. The default is `FALSE`.
- `alpha`: A numeric argument specifying the significance level for calculating confidence intervals. The function will return confidence intervals with coverage `1 - alpha`. The default is `0.05`.

## Value
The function returns a named list with the following components:
- `att_hat`: The estimated overall average treatment effect on the treated.
- `att_se`: The standard error for the overall average treatment effect. If `indep_counts` is provided, this standard error is asymptotically exact. Otherwise, it is asymptotically conservative.
- `catt_hats`: A named vector containing the estimated average treatment effects for each cohort.
- `catt_ses`: A named vector of standard errors for the cohort average treatment effects, if `q` is less than `1`.
- `cohort_probs`: A vector of the estimated probabilities of being in each cohort, conditional on treatment. These probabilities are calculated from the counts of units in each cohort.
- `catt_df`: A data frame displaying cohort names, cohort average treatment effects, standard errors, and confidence intervals.
- `beta_hat`: The estimated vector of coefficients.
- `treat_inds`: The indices of the coefficients in `beta_hat` corresponding to the treatment effects for each cohort at each time.
- `treat_int_inds`: The indices of the coefficients in `beta_hat` corresponding to the interactions between the treatment effects and covariates.
- `sig_eps_sq`: The provided or estimated variance of row-level IID noise.
- `sig_eps_c_sq`: The provided or estimated variance of unit-level IID noise.
- `lambda.max`, `lambda.min`: The selected or provided values for the penalty parameter lambda.
- `lambda.max_model_size`, `lambda.min_model_size`, `lambda_star_model_size`: The sizes of the models corresponding to these penalty parameters.
- `lambda_star`: The penalty parameter lambda selected by BIC.
- `X_ints`: The design matrix containing all interactions, time and cohort dummies, and other variables.
- `y`: The vector of response variables.
- `X_final`, `y_final`: The transformed design matrix and response vector after applying adjustments.
- `N`, `T`, `R`, `d`, `p`: The number of units, time periods, treated cohorts, covariates, and features in the final dataset.


## References
- Faletto, G (2024). *Fused Extended Two-Way Fixed Effects for Difference-in-Differences with Staggered Adoptions*. [arXiv preprint arXiv:2312.05985](https://arxiv.org/abs/2312.05985).
- Pesaran, M. H. (2015). *Time Series and Panel Data Econometrics*. Oxford University Press. [URL](https://ideas.repec.org/b/oxp/obooks/9780198759980.html).

