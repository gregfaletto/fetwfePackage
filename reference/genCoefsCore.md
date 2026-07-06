# Generate Coefficient Vector for Data Generation (core)

This function generates a coefficient vector `beta` along with a sparse
auxiliary vector `theta` for simulation studies of the fused extended
two-way fixed effects estimator. The returned `beta` is formatted to
align with the design matrix created by
[`simulateDataCore()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateDataCore.md),
and is a valid input for the `beta` argument of that function. The
vector `theta` is sparse, with nonzero entries occurring with
probability `density` and scaled by `eff_size`. See the simulation
studies section of Faletto (2025) for details.

## Usage

``` r
genCoefsCore(
  G = NULL,
  T,
  d,
  density,
  eff_size,
  fusion_structure = c("cohort", "event_study"),
  n_signal_cohorts = NULL,
  treat_base_levels = NULL,
  seed = NULL,
  R = NULL
)
```

## Arguments

- G:

  Integer. The number of treated cohorts (treatment is assumed to start
  in periods 2 to `G + 1`). Defaults to `NULL`; supply either `G` or the
  deprecated alias `R` (described below).

- T:

  Integer. The total number of time periods.

- d:

  Integer. The number of time-invariant covariates. If `d > 0`,
  additional terms corresponding to covariate main effects and
  interactions are included in `beta`.

- density:

  Numeric in (0,1\]. The probability that any given entry in the initial
  coefficient vector `theta` is nonzero. `density = 1` gives a fully
  dense (non-sparse) coefficient vector.

- eff_size:

  Numeric. The magnitude used to scale nonzero entries in `theta`. Each
  nonzero entry is set to `eff_size` or `-eff_size` (with a 60 percent
  chance for a positive value).

- fusion_structure:

  Character. One of `"cohort"` (default) or `"event_study"`. Selects the
  inverse fusion transform applied to the treatment-effect block,
  controlling the basis in which the true treatment effects are sparse.
  `"cohort"` is byte-identical to previous behavior; `"event_study"`
  fuses effects at the same time since treatment across cohorts. See
  [`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
  for details.

- n_signal_cohorts, treat_base_levels:

  (Optional) Targeted-sparsity mode (#332). Both `NULL` (the default) is
  the unchanged, byte-identical uniform-`density` path. Otherwise the
  treatment signal is placed deterministically on the per-cohort fused
  base levels for a sparse, non-degenerate, heterogeneous
  high-dimensional DGP. See
  [`genCoefs()`](https://gregfaletto.github.io/fetwfePackage/reference/genCoefs.md)
  for the full description; the two are mutually exclusive.

- seed:

  (Optional) Integer. Seed for reproducibility. `NA` (or `NULL`) means
  "draw from the ambient random-number generator" — no
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called.

- R:

  Deprecated. The former name for `G`; still accepted with a deprecation
  warning, and will be removed in a future release. Use `G`.

## Value

A list with two elements:

- `beta`:

  A numeric vector representing the full coefficient vector after the
  inverse fusion transform.

- theta:

  A numeric vector representing the coefficient vector in the
  transformed feature space. `theta` is a sparse vector, which aligns
  with an assumption that deviations from the restrictions encoded in
  the FETWFE model are sparse. `beta` is derived from `theta`.

## Details

The length of `beta` is given by \$\$p = G + (T - 1) + d + dG +
d(T - 1) + \mathit{num\\treats} + (\mathit{num\\treats} \times d)\$\$,
where the number of treatment parameters is defined as
\$\$\mathit{num\\treats} = T \times G - \frac{G(G+1)}{2}\$\$.

The function operates in two steps:

1.  It first creates a sparse vector `theta` of length \\p\\, with
    nonzero entries occurring with probability `density`. Nonzero
    entries are set to `eff_size` or `-eff_size` (with a 60\\

2.  The full coefficient vector `beta` is then computed by applying an
    inverse fusion transform to `theta` using internal routines:
    `genBackwardsInvFusionTransformMat()` for the fixed-effect blocks
    and, for the treatment-effect block,
    `genInvTwoWayFusionTransformMat()` when
    `fusion_structure = "cohort"` or
    `genInvEventStudyFusionTransformMat()` when
    `fusion_structure = "event_study"`.

## References

Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions. *arXiv preprint
arXiv:2312.05985*. <https://arxiv.org/abs/2312.05985>.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Set parameters for the coefficient generation
  G <- 3         # Number of treated cohorts
  T <- 6         # Total number of time periods
  d <- 2         # Number of covariates
  density <- 0.1 # Probability that an entry in the initial vector is nonzero
  eff_size <- 1.5  # Scaling factor for nonzero coefficients
  seed <- 789    # Seed for reproducibility

  # Generate coefficients using genCoefsCore()
  coefs_core <- genCoefsCore(G = G, T = T, d = d, density = density,
  eff_size = eff_size, seed = seed)
  beta <- coefs_core$beta
  theta <- coefs_core$theta

  # For diagnostic purposes, compute the expected length of beta.
  # The length p is defined internally as:
  #   p = G + (T - 1) + d + d*G + d*(T - 1) + num_treats + num_treats*d,
  # where num_treats = T * G - (G*(G+1))/2.
  num_treats <- T * G - (G * (G + 1)) / 2
  p_expected <- G + (T - 1) + d + d * G + d * (T - 1) + num_treats + num_treats * d

  cat("Length of beta:", length(beta), "\nExpected length:", p_expected, "\n")
} # }
```
