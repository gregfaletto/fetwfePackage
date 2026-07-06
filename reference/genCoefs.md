# Generate Coefficient Vector for Data Generation

This function generates a coefficient vector `beta` for simulation
studies of the fused extended two-way fixed effects estimator. It
returns an S3 object of class `"FETWFE_coefs"` containing `beta` along
with simulation parameters `G`, `T`, and `d`. See the simulation studies
section of Faletto (2025) for details.

## Usage

``` r
genCoefs(
  G = NULL,
  T,
  d,
  density,
  eff_size,
  fusion_structure = c("cohort", "event_study"),
  assignment_type = c("marginal", "multinomial", "ordered"),
  assignment_strength = 1,
  assignment_interactions = NULL,
  assignment_interaction_strength = NULL,
  n_signal_cohorts = NULL,
  treat_base_levels = NULL,
  seed = NULL,
  verbose = FALSE,
  R = NULL
)
```

## Arguments

- G:

  Optional integer. The number of treated cohorts (treatment is assumed
  to start in periods 2 to `G + 1`). Defaults to `NULL`; supply either
  `G` or the deprecated alias `R` (described below).

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

  Character. One of `"cohort"` (default) or `"event_study"`. Controls
  the basis in which the true treatment-effect coefficients are sparse.
  Under `"cohort"` the sparse transformed coefficients `theta` are
  mapped back through the default two-way (cohort) inverse fusion
  transform (byte-identical to previous behavior); under `"event_study"`
  they are mapped through the event-study inverse fusion transform, so
  the true treatment effects are sparse in the event-study basis
  (effects sharing the same time since treatment, \\e = t - g\\, tend to
  be equal across cohorts). This is the simulation-side companion to the
  `fusion_structure` argument of
  [`fetwfe()`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md).

- assignment_type:

  Character. One of `"marginal"` (default), `"multinomial"`, or
  `"ordered"`. Selects the data-generating process for cohort assignment
  in
  [`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md).

  - `"marginal"`: each unit's cohort is drawn uniformly, independent of
    its covariates. Original pre-1.14.0 behavior, preserved
    byte-identically.

  - `"multinomial"`: cohort assignment follows a multinomial-logit
    propensity-score model \\\pi_g(x) = \exp(\gamma_g^\top x) /
    \sum\_{g'} \exp(\gamma\_{g'}^\top x)\\ with \\\gamma_0 \equiv 0\\
    (never-treated reference). Requires \\d \ge 1\\.

  - `"ordered"`: cohort assignment follows a proportional-odds
    (ordered-logit) model with cumulative probabilities \\P(W \le g
    \| x) = \mathrm{plogis}(\alpha_g - \gamma^\top x)\\. Cutpoints
    \\\alpha_g\\ are chosen so the marginal cohort probabilities are
    approximately uniform 1/(G + 1). Requires \\d \ge 1\\.

- assignment_strength:

  Non-negative numeric scalar. Scales the logit coefficients in the
  propensity-score model. `0` reduces both non-marginal types to the
  uniform marginal distribution by construction. Larger values produce
  stronger covariate-cohort coupling. Defaults to `1.0`. Ignored when
  `assignment_type = "marginal"`.

- assignment_interactions:

  Optional. A list of length-2 integer vectors, each naming a pair of
  covariate indices \\(j, k)\\ in \\\[1, d\]\\ whose elementwise product
  \\x\_{i,j} \cdot x\_{i,k}\\ enters the propensity model as an
  additional column. Self-interactions `c(j, j)` are allowed and yield a
  quadratic term \\x_j^2\\. Unordered pairs are canonicalized to
  `c(min(j, k), max(j, k))` and duplicates are silently deduplicated
  (inspect `coefs$assignment_coefs$interactions` on the returned object
  to verify the retained list). The interaction columns enter the
  propensity model only; the outcome model continues to use the original
  covariates. Defaults to `NULL` (no interactions — v1.14.0 behavior is
  preserved byte-identically). Passing
  [`list()`](https://rdrr.io/r/base/list.html) (empty list) is treated
  as equivalent to `NULL` — no interactions specified, behavior is
  identical to the v1.14.0 marginal-cohort path within `multinomial` /
  `ordered` DGPs. Errors when `assignment_type = "marginal"` (the
  marginal DGP has no propensity model to augment). New in 1.14.1.

- assignment_interaction_strength:

  Optional non-negative numeric scalar. Scales the Gaussian draws of the
  interaction coefficients independently of `assignment_strength`.
  Defaults to `NULL`, which means "fall through to
  `assignment_strength`". Useful when stress-testing whether the
  nonlinear-propensity angle alone drives downstream differences
  (without simultaneously cranking the linear angle). Ignored when
  `assignment_interactions = NULL`. New in 1.14.1.

- n_signal_cohorts:

  (Optional) Targeted-sparsity mode (#332). A single positive integer
  `k` in `1..(G - 1)`. When supplied, the coefficient vector is built
  *deterministically* with the treatment signal placed on exactly `k`
  cohorts' fused base levels (and the covariate / interaction blocks
  left at zero), giving `||theta||_0 = k` non-zeros with a guaranteed
  non-zero, *heterogeneous* cohort effect (overall ATT \\\neq 0\\ and
  positive cohort-weight variance \\V_2 \> 0\\). This produces the
  simultaneously sparse and non-degenerate high-dimensional (\\p \ge
  NT\\) data-generating process the desparsified / nodewise debiased-ATT
  coverage study needs, which the default uniform `density` placement
  cannot (it almost never lands signal on the tiny treatment-effect
  coordinate block). Cohort 1 is held at the baseline and cohorts
  `2..(k + 1)` receive the distinct magnitudes
  `eff_size * (1, 2, ..., k)` (so `eff_size` must be non-zero). Defaults
  to `NULL`; `NULL` (with `treat_base_levels = NULL`) is the unchanged,
  byte-identical uniform-`density` behavior. Mutually exclusive with
  `treat_base_levels`.

- treat_base_levels:

  (Optional) Targeted-sparsity mode (#332), explicit form. A finite
  numeric vector of length `G` giving each cohort's fused base level
  directly (placed on `getTreatInds()[getFirstInds()]`); every other
  coordinate is left at zero, so `||theta||_0` is the number of non-zero
  entries. Note these are *fused* base levels, not the per-cohort ATTs:
  under `fusion_structure = "cohort"` a vector `(b_1, ..., b_G)` maps to
  the cumulative cohort effects `cumsum(b)` (e.g. `c(0, m, 0)` gives
  cohort ATTs `(0, m, m)`), and `"event_study"` uses a different but
  still heterogeneity-preserving map. The result must be non-degenerate
  and heterogeneous (\\V_2 \> 0\\); an all-zero or constant-effect
  specification is rejected. Defaults to `NULL`. Mutually exclusive with
  `n_signal_cohorts`.

- seed:

  (Optional) Integer. Seed for reproducibility. Three deterministic
  offsets share this seed: the main coefficient draw uses `seed`; the
  assignment coefficients use `seed + 1L`; the Monte Carlo integration
  in
  [`getTes()`](https://gregfaletto.github.io/fetwfePackage/reference/getTes.md)
  uses `seed + 2L`. `NA` (or `NULL`) means "draw from the ambient
  random-number generator" — no
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called, so any
  preceding [`set.seed()`](https://rdrr.io/r/base/Random.html) is
  respected and the result varies across calls.

- verbose:

  Logical. If `TRUE`, emit a
  [`message()`](https://rdrr.io/r/base/message.html) when
  `assignment_interactions` canonicalization removes duplicate or
  unordered pairs (e.g., when the user passes both `c(1, 2)` and
  `c(2, 1)`). Default `FALSE` (silent — users can verify the final
  canonical list via `coefs$assignment_coefs$interactions`).

- R:

  Deprecated. The former name for `G`; still accepted with a deprecation
  warning, and will be removed in a future release. Use `G`.

## Value

An object of class `"FETWFE_coefs"`, which is a list containing:

- beta:

  A numeric vector representing the full coefficient vector after the
  inverse fusion transform.

- theta:

  A numeric vector representing the coefficient vector in the
  transformed feature space. `theta` is a sparse vector, which aligns
  with an assumption that deviations from the restrictions encoded in
  the FETWFE model are sparse. `beta` is derived from `theta`.

- fusion_structure:

  The fusion structure (`"cohort"` or `"event_study"`) used to build the
  treatment-effect coefficients.

- G:

  The provided number of treated cohorts.

- R:

  Deprecated alias for `G`, retained for backward compatibility;
  populated with the same value. Use `G`. Will be removed in a future
  release.

- T:

  The provided number of time periods.

- d:

  The provided number of covariates.

- seed:

  The provided seed.

- assignment_type:

  The selected cohort-assignment DGP (`"marginal"` / `"multinomial"` /
  `"ordered"`). New in 1.14.0.

- assignment_strength:

  The scaling factor applied to the assignment coefficients. New in
  1.14.0.

- assignment_interaction_strength:

  The scaling factor applied to the interaction coefficients. `NULL`
  when no interactions were specified or when the user passed `NULL`
  (the fall-through default). New in 1.14.1.

- assignment_coefs:

  `NULL` when `assignment_type = "marginal"`; otherwise a list with
  elements `type`, `strength`, `coefs` (the gamma matrix or vector), and
  (for ordered) `cutpoints`. Starting in 1.14.1, `assignment_coefs` also
  carries the sub-slots `interactions` (the canonicalized + deduplicated
  list of pairs, or `NULL`), `delta` (the interaction coefficient matrix
  for multinomial or vector for ordered, or `NULL`), and
  `interaction_strength` (the effective scaling factor used for the
  `delta` draws, or `NULL` when no interactions). New in 1.14.0;
  `interactions`, `delta`, and `interaction_strength` sub-slots new in
  1.14.1.

## Details

Optional arguments `assignment_type` and `assignment_strength` control
whether cohort membership in the simulated panel is drawn marginally
(independent of the covariates, the original behavior) or from a
covariate-dependent propensity-score model — either a multinomial-logit
or an ordered-logit (proportional-odds) model. The default
`assignment_type = "marginal"` preserves the pre-1.14.0 behavior
byte-identically. See
[`vignette("simulation_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/simulation_vignette.md)
for worked examples.

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

The multinomial-logit and proportional-odds reference DGPs are the
canonical parametric propensity-score models named in Faletto (2025)
line 1016; the propensity-weighted population-truth aggregation matches
Eq. `att.estimator.weighted` (line 837).

Note on the random-number generator: passing an explicit numeric `seed`
calls `set.seed(seed)` internally and **leaves the global RNG advanced**
after the call returns. This is deliberate — it makes a simulation both
reproducible (the same `seed` always yields the same coefficients) and
varying (drawing data afterwards consumes the advanced stream). To draw
from / preserve the ambient random stream instead — without calling
[`set.seed()`](https://rdrr.io/r/base/Random.html) — pass `seed = NA`
(or `seed = NULL`).

## References

Faletto, G (2025). Fused Extended Two-Way Fixed Effects for
Difference-in-Differences with Staggered Adoptions. *arXiv preprint
arXiv:2312.05985*. <https://arxiv.org/abs/2312.05985>.

## See also

[`fetwfe`](https://gregfaletto.github.io/fetwfePackage/reference/fetwfe.md),
whose `fusion_structure` argument this mirrors on the estimation side;
[`vignette("fusion_structure_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/fusion_structure_vignette.md)
for the cohort-vs-event-study distinction, and
[`vignette("simulation_vignette", package = "fetwfe")`](https://gregfaletto.github.io/fetwfePackage/articles/simulation_vignette.md)
for the full simulation pipeline.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate coefficients
  coefs <- genCoefs(G = 5, T = 30, d = 12, density = 0.1, eff_size = 2, seed = 123)

  # Simulate data using the coefficients
  sim_data <- simulateData(coefs, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)

  # Event-study-sparse truth: treatment effects that share the same time
  # since treatment are fused across cohorts (the simulation-side companion
  # to fetwfe()'s fusion_structure = "event_study").
  coefs_es <- genCoefs(
    G = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
    fusion_structure = "event_study", seed = 123
  )

  # Covariate-dependent cohort assignment: multinomial-logit DGP
  coefs_mn <- genCoefs(
    G = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
    assignment_type = "multinomial", assignment_strength = 1.0,
    seed = 123
  )
  sim_mn <- simulateData(coefs_mn, N = 120, sig_eps_sq = 5, sig_eps_c_sq = 5, seed = 123)

  # Covariate-dependent cohort assignment with nonlinear propensity
  # (multinomial-logit + a single x1*x2 interaction term in the propensity
  # model only; outcome model continues to use plain X):
  coefs_int <- genCoefs(
    G = 5, T = 30, d = 12, density = 0.1, eff_size = 2,
    assignment_type = "multinomial",
    assignment_interactions = list(c(1, 2)),
    assignment_interaction_strength = 1.5,
    seed = 123
  )
} # }
```
