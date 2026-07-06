# Simulation Vignette for FETWFE: From Coefficients to True Treatment Effects

``` r

# Load necessary libraries
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(fetwfe)
```

## Introduction

This vignette demonstrates how to conduct simulation studies with the
[fetwfe](https://gregfaletto.github.io/fetwfePackage/) package. In
particular, we will:

- **Generate a vector of coefficients** with
  [`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md).
  These coefficients produce unit‐and‐time specific responses that
  respect difference‐in‐differences assumptions (e.g., conditional
  parallel trends) and the sparsity assumptions behind FETWFE.
- **Simulate a panel data set** with
  [`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md).
- **Estimate treatment effects** via the FETWFE estimator with
  [`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md).
- **Extract the true treatment effects** using
  [`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md).

The workflow here follows the simulation‐study design outlined in [the
paper](https://arxiv.org/abs/2312.05985), so you may wish to skim its
setup section for additional context.

## Simulation Workflow Using Piping

Below is a complete simulation pipeline, step by step.

### Step 1: Generate Simulation Coefficients

The
[`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
function returns an object of class `"FETWFE_coefs"`, containing both
the coefficient vector and its simulation parameters. The
`fusion_structure` argument (mirroring
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md))
additionally controls the *basis* in which the true treatment effects
are sparse — `"cohort"` (the default) or `"event_study"`; see the
vignette
[`vignette("fusion_structure_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/fusion_structure_vignette.md).
In this example we set:

- **G** (number of treated cohorts) = 3
- **T** (number of time periods) = 6
- **d** (number of covariates) = 2  
- **density** (sparsity level) = 0.1  
- **eff_size** (effect‐size multiplier) = 2

``` r

# Generate the coefficient object for simulation
sim_coefs <- genCoefs(
  G         = 3, 
  T         = 4, 
  d         = 2, 
  density   = 0.1, 
  eff_size  = 2, 
  seed      = 101
)
```

(Again, for more details on the meaning of these parameters, see the
simulation study section of [the
paper](https://arxiv.org/abs/2312.05985).)

### Step 2: Simulate Panel Data

Next, we simulate a panel data set using the generated coefficient
object with the
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)
function. With
[`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md),
we generate:

- `N` units, each assigned to one of the cohorts
- Time‐invariant covariates drawn from a specified distribution
- Outcomes at times 1 through T, using our simulated coefficients

Here we choose:

- `N` (number of units) as 60,
- `sig_eps_sq` (observation-level noise variance) as 1,
- `sig_eps_c_sq` (unit-level noise variance) as 1, and
- use the default `"gaussian"` distribution for the covariates.

``` r

# Simulate panel data based on the coefficients
sim_data <- simulateData(
  sim_coefs,
  N = 60,
  sig_eps_sq = 1,
  sig_eps_c_sq = 1,
  distribution = "gaussian",
  seed = 101
  )
```

The dataframe is stored in `sim_data$pdata`, so we can take a quick look
at the results:

``` r

head(sim_data$pdata)
#>   time   unit treatment            y       cov1      cov2
#> 1    1 unit01         0  0.001990395  0.4061679 0.1262211
#> 2    2 unit01         0 -0.019280260  0.4061679 0.1262211
#> 3    3 unit01         0 -1.424784613  0.4061679 0.1262211
#> 4    4 unit01         0 -0.737871496  0.4061679 0.1262211
#> 5    1 unit02         0  0.203581047 -0.5242428 0.8761651
#> 6    2 unit02         0  0.213280262 -0.5242428 0.8761651
```

### Step 3: Run the FETWFE Estimator on Simulated Data

We then run the estimator on the simulated data using
[`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md).
(We could get the same results by manually unpacking `sim_data` and
passing the arguments appropriately to `fewtfe()`.
[`fetwfeWithSimulatedData()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfeWithSimulatedData.md)
is just a wrapper function that takes care of this for us.)

``` r

result <- fetwfeWithSimulatedData(sim_data)
```

We can now extract the results from `result` in the same way that we can
with the standard
[`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md)
function.

``` r

summary(result)
#> Summary of Fused Extended Two-Way Fixed Effects
#> ================================================
#> 
#> Overall ATT: 0.0843  (SE = 0.1372, p = 0.5389, 95% CI = [-0.1845, 0.3531])
#> Selected: TRUE
#> 
#> CATT (preview) [simultaneous 95% CI]:
#>  cohort   estimate        se     ci_low    ci_high      p_value selected
#>       2 -1.3369236 0.1972001 -1.7779399 -0.8959073 2.411404e-11     TRUE
#>       3  0.9165759 0.1264333  0.6338217  1.1993300 8.366641e-13     TRUE
#>       4  0.0000000 0.0000000  0.0000000  0.0000000           NA    FALSE
#> 
#> Event Study (preview) [simultaneous 95% CI]:
#>  event_time n_cohorts   estimate        se     ci_low   ci_high      p_value
#>           0         3  0.0000000 0.0000000  0.0000000  0.000000           NA
#>           1         2  0.5095182 0.3915374 -0.3636595  1.382696 3.417505e-01
#>           2         1 -2.0053854 0.2958002 -2.6650570 -1.345714 2.411404e-11
#> 
#> Model Details:
#>   Units (N)           : 60
#>   Time periods (T)    : 4
#>   Treated cohorts (G) : 3
#>   Covariates (d)      : 2
#>   Features (p)        : 38
#>   Selected size       : 4
#>   Lambda*             : 0.0984
```

### Step 4: Extract True Treatment Effects

To evaluate the estimated ATT, we can compute the true treatment effects
using the original coefficient object. The
[`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md)
function extracts both the overall average treatment effect and the
cohort-specific effects.

``` r

# Extract the true treatment effects
true_tes <- getTes(sim_coefs)

# Print the true overall treatment effect
cat("True Overall ATT:", true_tes$att_true, "\n")
#> True Overall ATT: -0.1111111

# Print the cohort-specific treatment effects
print(true_tes$actual_cohort_tes)
#> [1] -1.333333  1.000000  0.000000
```

We can use this to calculate metrics to evaluate our estimated treatment
effect, like squared error:

``` r

squared_error <- (result$att_hat - true_tes$att_true)^2

cat("Squared error of ATT estimate:", squared_error, "\n")
#> Squared error of ATT estimate: 0.03817985
```

### Combining the Workflow in One Pipeline

You can also chain the simulation functions together with the pipe
operator. The following code generates the coefficients, simulates the
data, and runs the estimator all in one pipeline:

``` r

coefs <- genCoefs(G = 3, T = 4, d = 2, density = 0.1, eff_size = 2, seed = 2025)

result_piped <- coefs |>
  simulateData(N = 60, sig_eps_sq = 1, sig_eps_c_sq = 1, seed = 2025) |>
  fetwfeWithSimulatedData()

cat("Estimated Overall ATT from piped workflow:", result_piped$att_hat, "\n")
#> Estimated Overall ATT from piped workflow: -0.7183584

true_tes_piped <- coefs |> getTes()

# Print the true overall treatment effect
cat("True Overall ATT:", true_tes_piped$att_true, "\n")
#> True Overall ATT: -0.6666667

# Print the squared estimation error
squared_error_piped = (result_piped$att_hat - true_tes_piped$att_true)^2

cat("Squared estimation error:", squared_error_piped, "\n")
#> Squared estimation error: 0.00267203
```

## Covariate-dependent cohort assignment

The simulator workflow above generates panels in which cohort assignment
is independent of the covariates. Each unit is assigned to one of the
$`G + 1`$ cohorts with uniform probability $`1/(G + 1)`$, regardless of
its covariate values. This is the simplest data-generating process (DGP)
for cohort assignment, but in many realistic settings, cohort membership
is correlated with observable characteristics — adoption decisions
depend on firm size, industry, region, and so on.

Version 1.14.0 of [fetwfe](https://gregfaletto.github.io/fetwfePackage/)
introduces two covariate-dependent DGPs to
[`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
via the new `assignment_type` and `assignment_strength` arguments:

- `assignment_type = "marginal"` (default): the original uniform
  marginal DGP. Preserves pre-1.14.0 behavior exactly.
- `assignment_type = "multinomial"`: a multinomial-logit
  propensity-score model. Each unit’s cohort probability is
  ``` math
  \pi_g(x) = \frac{\exp(\gamma_g^\top x)}{\sum_{g' = 0}^{G} \exp(\gamma_{g'}^\top x)}, \quad g = 0, 1, \ldots, G,
  ```
  with $`\gamma_0 \equiv 0`$ for the never-treated reference cohort and
  $`\gamma_g`$ drawn from a Gaussian (scaled by `assignment_strength`)
  for each treated cohort $`g`$.
- `assignment_type = "ordered"`: an ordered-logit (proportional-odds)
  model that exploits the natural ordering of the adoption times. The
  ordinal scale runs through treated cohorts in temporal-adoption order
  with never-treated at the top: slot 1 is cohort 1 (the
  earliest-adopting cohort), slot 2 is cohort 2, …, slot $`G`$ is cohort
  $`G`$ (the latest-adopting cohort), and slot $`G + 1`$ is
  never-treated. Cumulative probabilities follow the standard McCullagh
  proportional-odds parameterization with the subtraction convention,
  ``` math
  P(W' \le g \mid x) = \mathrm{plogis}(\alpha_g - \gamma^\top x), \quad g = 1, \ldots, G,
  ```
  with a single shared $`\gamma`$ vector and cutpoints $`\alpha_g`$
  chosen so the marginal cohort probabilities are approximately uniform
  $`1/(G + 1)`$. As $`\gamma^\top x`$ increases monotonically, the modal
  cohort traverses the ordinal scale in order: cohort 1 $`\to`$ cohort 2
  $`\to \ldots \to`$ cohort $`G`$$`\to`$ never-treated. Low
  $`\gamma^\top x`$ concentrates mass at the *bottom* of the scale
  (earliest-adopting cohorts); high $`\gamma^\top x`$ concentrates mass
  at the *top* (never-treated, with later cohorts at intermediate
  values). $`\gamma`$ is the propensity for being *later* on the
  treatment-timing scale.

Both DGPs are the canonical reference models named in the FETWFE paper
(Faletto 2025, line 1016).

### Multinomial-logit DGP

We generate a panel where adoption probabilities depend on the
covariates via a multinomial-logit model. The `assignment_strength`
parameter scales the logit coefficients: at `assignment_strength = 0`,
the DGP collapses to the marginal case (uniform $`1/(G + 1)`$
probabilities, regardless of $`X`$); larger values produce stronger
covariate-cohort coupling.

``` r

sim_coefs_mn <- genCoefs(
  G                   = 3,
  T                   = 4,
  d                   = 2,
  density             = 0.1,
  eff_size            = 2,
  assignment_type     = "multinomial",
  assignment_strength = 1.0,
  seed                = 101
)
```

The returned `FETWFE_coefs` object now carries three new slots:

``` r

sim_coefs_mn$assignment_type
#> [1] "multinomial"
sim_coefs_mn$assignment_strength
#> [1] 1
str(sim_coefs_mn$assignment_coefs)
#> List of 6
#>  $ type                : chr "multinomial"
#>  $ strength            : num 1
#>  $ coefs               : num [1:2, 1:3] 0.181 0.785 -1.353 1.983 1.238 ...
#>  $ interactions        : NULL
#>  $ delta               : NULL
#>  $ interaction_strength: NULL
```

`assignment_coefs$coefs` is the $`d \times G`$ matrix of $`\gamma_g`$
vectors; column $`g`$ holds the propensity coefficients for treated
cohort $`g`$. We can simulate a panel and verify that the empirical
cohort proportions reflect covariate-cohort coupling rather than uniform
marginal assignment:

``` r

sim_mn <- simulateData(
  sim_coefs_mn,
  N            = 200,
  sig_eps_sq   = 1,
  sig_eps_c_sq = 1,
  seed         = 101
)
cat("Empirical cohort proportions:\n")
#> Empirical cohort proportions:
print(round(sim_mn$assignments / sum(sim_mn$assignments), 3))
#> [1] 0.285 0.180 0.285 0.250
```

We can fit FETWFE on this panel exactly as before:

``` r

fit_mn <- fetwfeWithSimulatedData(sim_mn, verbose = FALSE)
cat("Estimated overall ATT:", fit_mn$att_hat, "\n")
#> Estimated overall ATT: 0.1820036
```

### Ordered-logit DGP

The proportional-odds variant uses a single shared $`\gamma`$ vector and
$`G`$ cutpoints. It is appropriate when the cohorts have a natural
ordering — earlier adopters being “more like” later adopters than
never-treated units, for instance. The ordinal scale runs in temporal
order with never-treated at the top: slot 1 = cohort 1 (first treated),
…, slot $`G`$ = cohort $`G`$ (last treated), slot $`G + 1`$ =
never-treated. As $`\gamma^\top x`$ increases, the modal cohort moves up
the scale from cohort 1 toward never-treated.

``` r

sim_coefs_ord <- genCoefs(
  G                   = 3,
  T                   = 4,
  d                   = 2,
  density             = 0.1,
  eff_size            = 2,
  assignment_type     = "ordered",
  assignment_strength = 1.0,
  seed                = 101
)
```

``` r

str(sim_coefs_ord$assignment_coefs)
#> List of 7
#>  $ type                : chr "ordered"
#>  $ strength            : num 1
#>  $ coefs               : num [1:2] 0.181 0.785
#>  $ cutpoints           : num [1:3] -1.2532 -0.00765 1.23988
#>  $ interactions        : NULL
#>  $ delta               : NULL
#>  $ interaction_strength: NULL
```

Here `assignment_coefs$coefs` is a length-$`d`$ vector (the shared
$`\gamma`$) and `cutpoints` is the length-$`G`$ vector of cutpoints
$`\alpha_g`$. The cutpoints are chosen by root-finding on the
marginal-uniform condition: for each $`g \in 1, \ldots, G`$, the package
solves for $`\alpha_g`$ such that
$`E_X[\mathrm{plogis}(\alpha_g - \gamma^\top X)] = g / (G + 1)`$. This
guarantees that the marginal cohort proportions are approximately
uniform $`1/(G + 1)`$ for any choice of `assignment_strength` — the
strength controls the *conditional* coupling between $`X`$ and $`W`$,
not the marginal cohort proportions.

``` r

sim_ord <- simulateData(
  sim_coefs_ord,
  N            = 200,
  sig_eps_sq   = 1,
  sig_eps_c_sq = 1,
  seed         = 101
)
cat("Empirical cohort proportions (ordered):\n")
#> Empirical cohort proportions (ordered):
print(round(sim_ord$assignments / sum(sim_ord$assignments), 3))
#> [1] 0.315 0.220 0.250 0.215
```

### Effect of assignment strength

The `assignment_strength` parameter controls how strongly the covariates
predict cohort membership. Below is a side-by-side comparison of the
cohort distribution at three strength levels:

``` r

for (s in c(0.0, 1.0, 5.0)) {
  coefs_s <- genCoefs(
    G                   = 3,
    T                   = 4,
    d                   = 2,
    density             = 0.1,
    eff_size            = 2,
    assignment_type     = "multinomial",
    assignment_strength = s,
    seed                = 101
  )
  sim_s <- simulateData(
    coefs_s,
    N            = 200,
    sig_eps_sq   = 1,
    sig_eps_c_sq = 1,
    seed         = 101
  )
  cat(sprintf(
    "strength = %.1f -> proportions: %s\n",
    s,
    paste(
      round(sim_s$assignments / sum(sim_s$assignments), 3),
      collapse = ", "
    )
  ))
}
#> strength = 0.0 -> proportions: 0.32, 0.235, 0.215, 0.23
#> strength = 1.0 -> proportions: 0.285, 0.18, 0.285, 0.25
#> strength = 5.0 -> proportions: 0.3, 0.04, 0.36, 0.3
```

At `assignment_strength = 0` the cohort proportions are approximately
uniform $`1/(G + 1) = 0.25`$. At higher strength, the proportions
deviate from uniform because the covariate distribution drives the
assignment.

### Nonlinear propensity via interactions

Sometimes the propensity model needs to be nonlinear in the covariates —
for example, to stress-test estimators that assume a linear propensity.
The `assignment_interactions` argument injects multiplicative
interaction terms into the propensity model only; the outcome model
continues to use the original covariates. The argument accepts a list of
integer pairs naming the columns of $`X`$ whose product enters the
propensity:

``` r

sim_coefs_int <- genCoefs(
  G                               = 3,
  T                               = 4,
  d                               = 2,
  density                         = 0.1,
  eff_size                        = 2,
  assignment_type                 = "multinomial",
  assignment_strength             = 1.0,
  assignment_interactions         = list(c(1, 2)),
  assignment_interaction_strength = 1.0,
  seed                            = 101
)
str(sim_coefs_int$assignment_coefs$interactions)
#> List of 1
#>  $ : int [1:2] 1 2
str(sim_coefs_int$assignment_coefs$delta)
#>  num [1, 1:3] 0.896 0.254 0.55
```

Self-interactions are allowed and produce quadratic terms $`x_j^2`$:

``` r

sim_coefs_quad <- genCoefs(
  G                       = 3,
  T                       = 4,
  d                       = 2,
  density                 = 0.1,
  eff_size                = 2,
  assignment_type         = "multinomial",
  assignment_interactions = list(c(1, 1)),
  seed                    = 101
)
```

The interaction columns are constructed internally inside the propensity
model each time it is evaluated; the simulated `pdata` keeps just the
original $`d`$ covariate columns (the outcome design matrix is
unchanged). Pairs are canonicalized to `c(min, max)` and duplicates are
silently deduplicated.

### Truth derivation under non-marginal DGPs

Under the default marginal DGP, the overall true ATT is the
equal-weighted average of the cohort-specific effects:
``` math
\tau_\mathrm{ATT} = \frac{1}{G} \sum_{g = 1}^{G} \tau_\mathrm{ATT}(g).
```

Under covariate-dependent DGPs this aggregation must use cohort-share
weights instead, matching Faletto (2025) Eq. `att.estimator.weighted`
(line 837).
[`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md)
automatically computes the propensity-weighted aggregation when the
input `FETWFE_coefs` object was built with a non-marginal
`assignment_type`:

``` r

te_mn <- getTes(sim_coefs_mn)
cat("Multinomial-DGP cohort weights:\n")
#> Multinomial-DGP cohort weights:
print(round(te_mn$cohort_weights, 3))
#> [1] 0.219 0.426 0.355
cat("Propensity-weighted overall ATT:", round(te_mn$att_true, 4), "\n")
#> Propensity-weighted overall ATT: 0.1343
cat("(Uniform-weighted ATT for comparison:",
    round(mean(te_mn$actual_cohort_tes), 4), ")\n")
#> (Uniform-weighted ATT for comparison: -0.1111 )
```

The `cohort_weights` slot is new in 1.14.0. Under the marginal DGP it is
exactly uniform $`1/G`$; under the multinomial or ordered DGPs it
reflects $`E[\pi_g(X)] / \sum_{g' \text{ treated}} E[\pi_{g'}(X)]`$,
estimated by Monte Carlo integration over the covariate distribution.

The per-cohort CATT vector `actual_cohort_tes` is unchanged across DGPs
— it is intrinsic to the coefficient vector $`\beta`$ and does not
depend on how cohorts are assigned.

## Targeted-sparsity coefficients for high-dimensional DGPs

The `density` argument spreads the true signal uniformly across *all*
coefficient coordinates. That is fine for low-dimensional designs, but
when the design is **high-dimensional** ($`p \geq NT`$) — many
covariates and/or short panels — the treatment-effect coordinates are a
tiny fraction of the $`p`$ entries, so a density-based draw almost never
lands signal there: the overall ATT comes out (near) zero and the
data-generating process is *degenerate* for studying treatment-effect
inference.

Version 1.48.0 added a **targeted-sparsity** mode to
[`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
(the mutually-exclusive arguments `n_signal_cohorts` and
`treat_base_levels`) that builds the coefficient vector
*deterministically*, placing the treatment signal directly on a few
cohorts’ fused base levels and leaving every other coordinate at zero.
This builds a sparse, **non-degenerate, heterogeneous** truth with
positive cohort-weight variance $`V_2 > 0`$ (and, for the
`n_signal_cohorts` form, a guaranteed non-zero overall ATT — a free-form
`treat_base_levels` with mixed-sign levels could cancel to a zero ATT) —
the data-generating process behind the high-dimensional debiased-ATT
coverage study (see the high-dimensional section of
[`vignette("inference_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/inference_vignette.md)
and
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)).

`n_signal_cohorts = k` is the quick form: it places a heterogeneous
signal on `k` cohorts (cohort 1 stays at the baseline; cohorts
`2..(k + 1)` receive the fused base levels `eff_size * (1, 2, ..., k)`),
giving exactly `k` non-zero coordinates.

``` r

# A high-dimensional design (d = 12 covariates) with a controlled, sparse,
# heterogeneous truth -- signal on 2 cohorts' fused base levels:
hd_coefs <- genCoefs(
  G = 3, T = 5, d = 12,
  density = 0.2, # ignored in targeted-sparsity mode, but still required
  eff_size = 2,
  n_signal_cohorts = 2,
  seed = 2
)
te <- getTes(hd_coefs)
# A guaranteed non-zero, heterogeneous truth (which a density-based draw at this
# dimension almost never produces):
c(overall_ATT = te$att_true)
#> overall_ATT 
#>    2.666667
te$actual_cohort_tes # per-cohort: baseline 0, then two distinct effects
#> [1] 0 2 6
```

For full control, `treat_base_levels` takes an explicit length-`G`
vector of per-cohort *fused base levels* — not the per-cohort ATTs:
under `fusion_structure = "cohort"` a vector `(b_1, ..., b_G)` maps to
the cumulative cohort effects `cumsum(b)`, so `c(0, m, 0)` yields cohort
ATTs `(0, m, m)`. See
[`?genCoefs`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md).

## Conclusion

In this vignette, we walked through how to use the simulation functions
in the [fetwfe](https://gregfaletto.github.io/fetwfePackage/) package to
simulate data and run simulations similar to the ones in the simulation
studies section of [the FETWFE paper](https://arxiv.org/abs/2312.05985).
Version 1.14.0 adds covariate-dependent cohort-assignment DGPs
(multinomial-logit and ordered-logit) so Monte Carlo studies can
stress-test FETWFE against the realistic regime where cohort selection
depends on observable characteristics, and version 1.48.0 adds a
targeted-sparsity coefficient mode (`n_signal_cohorts` /
`treat_base_levels`) for controlled, non-degenerate high-dimensional
data-generating processes.

This pipeline streamlines simulation experiments so you can rapidly
evaluate FETWFE’s performance under varying scenarios. For more details,
consult the package documentation or reach out to the author.
