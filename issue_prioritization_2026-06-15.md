# Issue prioritization (paper-anchored) — 2026-06-15

A suggested prioritization of the open `fetwfePackage` issues,
**anchored on the FETWFE paper resubmission**. This is a
snapshot/suggestion as of 2026-06-15, not an authoritative roadmap —
reweight as priorities shift.

> **Status update (end of 2026-06-15 session).** Tier 1 **\#31 shipped**
> (PR \#294, v1.29.0) with follow-up **\#295** for coverage/`lambda_c`
> tuning. Tier 2 **\#146 shipped** (verified already done in v1.12.0,
> closed). **\#293 shipped** (PR \#296, v1.30.0 —
> [`simulateData()`](https://gregfaletto.github.io/fetwfePackage/reference/simulateData.md)
> can now make `p > NT` panels, unblocking high-dim ground-truth tests).
> **\#142 is SHIPPED** (updated 2026-06-17): Phase 1 fixed-`p` bootstrap
> (PR \#297) + Phase 2 high-dim desparsified bootstrap (PR \#298) +
> Phase 3 `event_study` (PR \#300) — all three merged, \#142 closed. The
> live thread has moved to the **high-dim debiased-band cluster** (#299
> / PR \#305, \#295, \#303, \#304). After that, the paper-anchored
> Tier-3 items **\#5 / \#30 Layer 2 / \#26 / \#33** below. Canonical
> sequencing doc: `.workflow/PRIORITIZATION.md`.

**Note on the top item:** the paper now contains a second,
high-dimensional regime (Theorem 6.6, the `p > NT` desparsified
debiased-ATT). \#31 is the package realization of it. A *validated
reference* implementation of the nodewise debiasing direction
(`riesz_lasso`, the ℓ1-penalized Riesz representer) already exists on
the paper side (`fetwfe` repo, `simulations/highdim_functions.R`), so
\#31 is a port of a tested algorithm rather than a from-scratch
derivation.

## Tier 1 — the live thread (paper-driven)

- **\#31** — Extend
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md)
  to `p ≥ NT` via a nodewise (desparsified-lasso) debiasing direction.
  The package realization of paper Theorem 6.6. Same estimator as the
  shipped fixed-`p`
  [`debiasedATT()`](https://gregfaletto.github.io/fetwfePackage/reference/debiasedATT.md);
  only the construction of the debiasing direction `v` changes
  (exact/ridged inverse → nodewise sparse inverse). Lift the `p < NT`
  guard.

## Tier 2 — paper ↔︎ package sync, small & high-leverage

- **\#146** — Surface unit-scaled `tilde_v_N` from the paper’s Item-4
  rate rewrite.
- **\#142** — Multiplier-bootstrap simultaneous CIs over
  treatment-effect families (the package side of the paper’s
  simultaneous event-study bands).

## Tier 3 — credibility / what reviewers & DiD users expect

- **\#5** — Support unbalanced panels (the most common real-world
  usability gap).
- **\#30** — Pre-trends diagnostic (Wald test and/or `pre_only` fit
  mode).
- **\#26** — Support anticipation periods.
- **\#33** — `predict.fetwfe()` for CATT at user-specified `x` (API
  completion).

## Tier 4 — broader features (larger, less urgent)

- **\#35** — Heterogeneous trends (relaxing parallel trends).
- **\#54** — OLS-after-selection FETWFE variant (lasso-select +
  OLS-refit).
- **\#34** — Generalized smooth-functional aggregation weights `f(π)`.
- **\#32** — User-supplied propensity-score estimator for CATT-x
  weighting.
- **\#161** — `fitCohortLogit()` for covariate-dependent cohort
  assignment.
- **\#7** — ATT-like estimators without cohort-membership probabilities.
- **\#4** — Accommodate unequal numbers of labeled and unlabeled units.

## Tier 5 — infrastructure / tuning / cleanup

- **\#279** — Tune `add_ridge` by CV over a grid.
- **\#283** — pkgdown site for the vignettes (discoverability).
- **\#201** — Remove the deprecated R cohort-count argument/field
  (follow-up to \#41).
- **\#155** — Deprecate/remove top-level duplicates of `$internal`
  slots.

## If the goal shifts from the resubmission to package adoption

Promote **\#5 (unbalanced panels)** to the top — it is the single most
common blocker for applied users. The rest of the ordering is otherwise
stable.
