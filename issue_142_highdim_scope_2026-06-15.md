# \#142 design + the high-dimensional scope decision — 2026-06-15

Builds on `issue_142_assessment_2026-06-15.md`. That doc recommended
**defer** for the *fixed-p* case (the analytic
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
already ships exact sup-t bands). The **high-dimensional angle flips
that**: at `p >= NT` there is no analytic band at all, so the multiplier
bootstrap is the *principled* path, not a redundant one. This doc
records the design and isolates the one decision that sets the scope.

------------------------------------------------------------------------

## THE DECISION (yours, as paper author)

**How should the high-dimensional (`p >= NT`) simultaneous bands be
scoped and labeled?**

### Option A — General theory → ship experimental *(my recommendation)*

Validity rests on the high-dim simultaneous-inference literature
(Belloni–Chernozhukov–Hansen 2014; Chernozhukov–Chetverikov–Kato 2013)
that Theorem 6.6 (`debiased.highdim.thm`) already builds on, but there
is **no FETWFE theorem for the joint/family case yet**. Build it; label
the `p >= NT` path **experimental** (exactly like high-dim
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md));
track coverage validation in a follow-up issue (the \#31 → \#295
pattern). - Pro: ships the actual motivation; honest; reuses the
established `#31` caveat machinery. - Con: a second experimental path
whose finite-sample coverage isn’t yet validated.

### Option B — Paper theorem covers it → ship supported

There is (or will be) a FETWFE theorem extending the high-dim debiased
CLT to joint families / sup-t calibration. Implement and ship `p >= NT`
as a **supported** feature, docs citing the theorem. - Needs: the
theorem label (so docs can cite it like `debiased.highdim.thm`). - Pro:
fully supported. Con: blocks on the theorem existing.

### Option C — Fixed-p only for now

Implement just the fixed-p bootstrap (literal \#142: `did`-parity,
large-K scaling, validated against the analytic `Sigma`); **defer**
high-dim until a paper-side theorem lands. Build `F` so it’s reusable
later. - Pro: lowest risk, fully grounded. Con: doesn’t deliver the
high-dim motivation now.

**Recommendation: A** — it ships the high-dim value (the reason to do
\#142 at all) while being honest about validation, mirroring how we
shipped high-dim
[`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md).
Pick **B** only if a joint-family high-dim theorem already exists (tell
me its label); pick **C** if you’d rather not add a second experimental
surface right now.

### Sub-questions (any option)

1.  Is the `event_study` family (which needs the new per-unit
    *propensity* IF — see below) in the first cut, or do we start with
    `cohort` / `all_post_treatment` / `custom` (where `Sigma_2 = 0`,
    regression IF only) and add `event_study` second?
2.  Defaults: `B = 1000`? multiplier weights `rademacher` vs `mammen`
    (two-point)?

------------------------------------------------------------------------

## Why high-dim is the real value (confirmed by grounding)

- [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
  **errors at `p >= NT`**: `getGramInv()` (`R/gls_machinery.R:316–375`)
  detects the singular Gram, returns `calc_ses = FALSE`, and
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
  stops with *“the Gram matrix on the selected support is not
  invertible”* (`R/simultaneous_cis.R:456`). So there is **no analytic
  sup-t band** in high-dim — the bootstrap is the only path.
- The high-dim literature calibrates max-statistics over debiased-lasso
  estimators via the **multiplier bootstrap on per-unit IFs**, precisely
  because the joint covariance of the nodewise/desparsified estimators
  is hard to pin down analytically. This is the same theory family
  Theorem 6.6 cites for the single-effect debiased ATT.
- Bonus: **\#293 (just merged) is exactly the data generator needed to
  test this** — it can now produce `p > NT` panels with a known
  coefficient vector.

------------------------------------------------------------------------

## Architecture (one design, both regimes)

The bootstrap perturbs a per-unit influence-function matrix `F`
(`N_units x K`, one column per effect in the family).
`F = F_reg + F_pi`.

**The bridge (reused, already exists).**
[`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md)
builds, for any family, `psi_tes_mat` (`K x num_treats`, via
`.build_psi_tes_for_family()`, `R/simultaneous_cis.R:846–924`) and then
`Psi <- t(psi_tes_mat %*% d_inv_treat_sel)` (`p_sel x K`,
`R/simultaneous_cis.R:486`) — column `k` of `Psi` **is** effect `k`’s
direction in selected coefficient space. This is the input the IF
construction needs.

**Regression channel `F_reg` (generalizes `debiasedATT`).** For effect
`k`: - debiasing direction `v_k`: - fixed-p:
`v_k <- solve(Sig + ridge, Psi[, k])` (`R/debiased_att.R:356`
pattern); - high-dim: `v_k <- riesz_lasso(Sig, Psi[, k], lambda_node_k)`
(`R/riesz_lasso.R:31`, reuse \#31). - per-row score
`score_k <- (X %*% v_k) * resid` and per-unit IF
`F_reg[, k] <- tapply(score_k, unit, sum)` — exactly `debiasedATT`’s
`score`/`unit_scores` lines (`R/debiased_att.R:408–413`), now `K`
columns. - Shared across `k`: `X`, `resid`, `Sig`, the unit index.
Per-effect: `Psi[,k]`, `v_k`, `lambda_node_k`. - Cost: `K` linear solves
(fixed-p) or `K` `riesz_lasso` coordinate-descent calls (high-dim,
`O(K · max_iter · p^2)`), `Sig` formed once.

**Propensity channel `F_pi` (the one NEW derivation).** The package has
only the *aggregate* `Sigma_2` (`theta' J_k' Sigma_pi_hat J_l theta`,
`R/variance_machinery.R:1259–1284`), no per-unit summand. But the
cohort-probability estimator is `pi_hat_g = (1/N) sum_i 1{c_i = g}`, so
its per-unit IF is `(e_{c_i} - pi_bar)` (`e` a unit indicator, `pi_bar`
the cohort-prob vector), whose empirical covariance is exactly
`Sigma_pi_hat = .multinomial_cov(...)`. Effect `k`’s propensity IF is
then `F_pi[i, k] = (J_k %*% theta_sel)' (e_{c_i} - pi_bar)`, with `J_k`
the per-effect Jacobian already built by `.build_j_list_for_family()`
(`R/simultaneous_cis.R:942–980`). **Only `event_study` has
`Sigma_2 != 0`**; `cohort` / `all_post_treatment` / `custom` give
`J_k = 0`, so `F_pi = 0` there.

**The bootstrap.** Draw `xi_i^{(b)}` (Rademacher / Mammen two-point),
`b = 1..B`;
`T^{(b)} = max_k |sum_i xi_i^{(b)} F[i,k]| / sqrt(N * v_N^{(k)})`;
critical value `c_hat = quantile_{1-alpha}({T^{(b)}})`; band
`estimate_k +/- c_hat * sqrt(v_N^{(k)} / (NT))`, where
`v_N^{(k)} = N * diag(Sigma)_k` (or `N * crossprod(F)[k,k]/N`).

------------------------------------------------------------------------

## Validation anchor (why fixed-p comes first)

In fixed-p there IS analytic ground truth, so the `F`-construction can
be *proven* correct: - `crossprod(F) / N_units ≈ Sigma` (the analytic
`Sigma_1 + Sigma_2`), and - bootstrap `c_hat ≈ qmvnorm` analytic `c`.

In high-dim neither exists. So: **build + validate `F` in fixed-p, then
reuse the identical construction in high-dim with `riesz_lasso`
directions.** This de-risks the part that can’t be checked directly.

------------------------------------------------------------------------

## Phasing

- **Phase 1 (fixed-p).** Build `F` (regression + propensity), the
  multiplier bootstrap, and `method = "bootstrap"` on
  [`simultaneousCIs()`](https://gregfaletto.github.io/fetwfePackage/reference/simultaneousCIs.md).
  Validate `crossprod(F)/N ≈ Sigma` and bootstrap `c ≈ qmvnorm c`.
  Delivers the literal \#142 (did-parity, large-K) AND the validated `F`
  machinery. *Supported.*
- **Phase 2 (high-dim).** Per-effect `riesz_lasso` `v_k` + the same
  bootstrap; this is where the analytic path errors, so the bootstrap is
  the only band. *Scope/label per the decision above.*

------------------------------------------------------------------------

## API (proposed)

    simultaneousCIs(result, family = c("event_study","cohort","all_post_treatment","custom"),
                    alpha = 0.05, contrasts = NULL,
                    method = c("analytic","bootstrap"), B = 1000, seed = NULL,
                    weights = c("rademacher","mammen"))

- Default `method = "analytic"` → byte-identical to today
  (back-compatible).
- `method = "bootstrap"` works in both regimes. At `p >= NT` the
  analytic path errors today, so the error message should point users to
  `method = "bootstrap"` (or we auto-switch with a note).
- Return: same S3 `"simultaneous_cis"`; bootstrap path fills
  `critical_value` with `c_hat` and adds `B` / `seed` / `method` fields.

------------------------------------------------------------------------

## Tests (planned)

- `crossprod(F)/N ≈ Sigma` (fixed-p; the F-correctness anchor) —
  regression-only and, for `event_study`, including `F_pi` vs `Sigma_2`.
- bootstrap `c ≈ qmvnorm c` at large `B` (fixed-p).
- `c_hat ∈ [qnorm(1-alpha/2) ≈ 1.96, Bonferroni]`; strictly tighter than
  Bonferroni under positive correlation.
- simultaneous coverage on a simulated DGP (~95% across the family).
- high-dim end-to-end on a `#293` `p > NT` fixture: bands produced,
  finite, deterministic given `seed`.
- determinism; mutation checks (e.g. wrong unit aggregation / missing
  `F_pi`).

------------------------------------------------------------------------

## Open methodological note

The per-unit propensity IF (`F_pi`) is the only derivation not already
in the package. It’s standard (the multinomial score), and it’s
checkable against the existing `Sigma_2`
(`crossprod(F_pi)/N ≈ Sigma_2`), so it’s low-risk — but it is new code,
and it’s the reason the `event_study` family is called out separately in
sub-question 1 above.
