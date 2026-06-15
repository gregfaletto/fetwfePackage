# Assessment: #142 (multiplier-bootstrap simultaneous CIs) — 2026-06-15

A working assessment, written before committing to an implementation. **Not** a decision —
a record of the concerns and the design space so the call can be made deliberately.

## TL;DR

- The package **already ships exact, principled simultaneous (sup-t) bands** via the
  analytic `simultaneousCIs()` (`R/simultaneous_cis.R`, ~1088 lines): the critical value
  is the `(1 - alpha)` quantile of `max_k |Z_k|`, `Z ~ N(0, cov2cor(Sigma))`, computed
  deterministically with `mvtnorm::qmvnorm()`. Families (`event_study`, `cohort`,
  `all_post_treatment`, `custom`), single-step max-T adjusted p-values, and a Bonferroni
  comparison are all already there.
- The multiplier bootstrap that #142 specifies is **asymptotically equivalent** to those
  bands. Its genuine, non-redundant value is narrower than the issue implies (see
  "Where the bootstrap actually differs").
- The **foundation #142 names is not in the package**: it assumes "unit-level influence
  functions (IFs) as first-class objects, or the bootstrap has nothing to perturb." But
  `simultaneousCIs()` builds its joint covariance `Sigma` **analytically** (Gram-inverse
  sandwich + multinomial delta-method), **not** from per-unit IFs. The per-unit IF matrix
  the bootstrap resamples would have to be built first.
- Net: #142 is a **moderate build with mostly-redundant output**. It's worth doing only
  for specific reasons (below); if those don't weigh heavily, a higher-leverage issue is a
  better use of the slot.

## What #142 asks for

A Chernozhukov–Chetverikov–Kato (2013) multiplier bootstrap over a family of
treatment-effect estimates:

1. Per-effect per-unit IF summand `f_i^{(k)}` for each effect `k` in the family.
2. Multiplier weights `xi_i^{(b)}` (Mammen / Rademacher), `b = 1..B`.
3. Sup-t statistic `T^{(b)} = max_k |(1/sqrt N) sum_i xi_i^{(b)} f_i^{(k)}| / sqrt(v_N^{(k)})`.
4. Critical value `c_hat = quantile_{1-alpha}({T^{(b)}})`.
5. Band: `theta_hat_k +/- c_hat * sqrt(v_N^{(k)} / (NT))`.

Plus tests (coverage; `c_hat` between pointwise `~1.96` and Bonferroni) and an
`inference_vignette.Rmd` section.

## Key finding: the analytic machinery is not IF-based

`simultaneousCIs()` assembles `Sigma = Sigma_1 + Sigma_2`:

- `Sigma_1` (regression channel) — `.assemble_joint_cov_var1()` from `Psi` (`p_sel x K`)
  and the **Gram inverse** `getGramInv()` (`R/simultaneous_cis.R:445-498`). This is the
  closed-form sandwich `Psi' G_inv Psi * sig_eps_sq`, not an average of per-unit terms.
- `Sigma_2` (cohort-probability channel) — `.assemble_joint_cov_var2()` from the
  **multinomial covariance** of the cohort proportions and per-effect Jacobians
  (`R/simultaneous_cis.R:499-525`).

There is **no `N_units x K` per-unit IF matrix anywhere** in this path. So the issue's
stated prerequisite ("without unit-level IFs as first-class objects, the bootstrap has
nothing to perturb") is literally unmet. The "Gaussianity upgrade" the issue depends on
(`#141`, tight Gaussian variance, landed in 1.12.0) delivered `hat_v_N` / `tilde_v_N` and
the variance *components* — but not the per-unit IFs.

## Why the bootstrap is (mostly) redundant here

The multiplier bootstrap with Rademacher/Mammen weights, as `B -> infinity`, reproduces
`max_k |Z_k|` for `Z ~ N(0, Sigma_emp)`, where `Sigma_emp = (1/N) sum_i f_i f_i'` is the
**empirical IF covariance**. That is the same object `qmvnorm()` already integrates
analytically. Under the paper's fixed-dim tight-Gaussianity (Theorem (c')), the MVN law is
the *exact* asymptotic law, so:

- The analytic bands are the principled answer and are **deterministic** (no Monte-Carlo
  noise, no seed).
- A bootstrap band is, at best, a noisier estimate of the same quantity.

## Where the bootstrap actually differs (the honest case *for* it)

Three real, non-redundant reasons — worth weighing:

1. **Heteroskedasticity / cluster robustness for free.** The default analytic `Sigma_1`
   uses the homoskedastic `sig_eps_sq * Psi' G_inv Psi`. The IF bootstrap's `Sigma_emp`
   uses *actual per-unit residual products*, so it is naturally cluster/heteroskedasticity-
   robust without the user choosing `se_type = "cluster"`. (Caveat: the package already
   offers a cluster-robust sandwich, so this overlaps with existing functionality.)
2. **Scaling in large families.** `qmvnorm()` integration gets slow / less accurate as `K`
   grows (dozens of event-study × cohort effects). A bootstrap scales linearly in `B*K` and
   stays stable. This is the most defensible advantage.
3. **`did`-parity / reviewer expectation.** `did::aggte(..., cband = TRUE)` is a bootstrap;
   a JoE reviewer asked for simultaneous bands. Shipping a bootstrap labelled like the
   competition has communication value even if the analytic band is as good or better.

None of these is a *statistical* necessity given the analytic bands; they are
robustness / scaling / familiarity arguments.

## Cost to build (moderate, reusable)

The per-unit IF matrix is assemblable from pieces that already exist:

- **Regression channel.** The IF of `psi_k' tau_hat` is
  `f_i^{reg,k} = psi_k' A_sel G_inv X_i' eps_i` (unit `i`'s rows on the selected support,
  times residuals). `debiasedATT()` already computes exactly this in **one** direction:
  `score <- (X %*% v) * resid; unit_scores <- tapply(score, unit, sum)`
  (`R/debiased_att.R`). Generalizing `v` (a vector) to a `K`-column direction matrix
  `D = (Psi' A_sel G_inv)'` yields the `N_units x K` regression-IF matrix directly.
- **Propensity channel.** `Sigma_2` already uses the per-cohort multinomial Jacobians; the
  per-unit propensity IF is `J_k' (e_{cohort_i} - pi)`, assemblable from the same Jacobians
  and each unit's cohort membership.

So: build `F = F_reg + F_pi` (`N_units x K`), then the bootstrap loop is ~20 lines. Estimate:
a few focused days incl. tests + vignette — not a multi-week effort, but real, and it adds a
second, Monte-Carlo, seed-dependent inference path to maintain alongside the analytic one.

## Options

1. **Full #142 — build the IF matrix + faithful CCK bootstrap.** Most work; exactly what the
   issue specifies; delivers the heteroskedasticity-robust + large-K-scalable bootstrap and
   `did`-parity. Adds a maintained Monte-Carlo path.
2. **Lighter — Gaussian multiplier on the existing analytic `Sigma`.** Draw
   `Z^{(b)} ~ N(0, Sigma)` and take the sup-t quantile. Cheap (reuses `Sigma`), gives a
   `cband`-style B-replication API matching the analytic bands in expectation — but it is
   **not** the per-unit-IF bootstrap the issue specifies (no heteroskedasticity robustness,
   no large-K speed win; it's just Monte-Carlo of what `qmvnorm` does exactly). Mostly a
   cosmetic / API-parity move. Hard to justify over the existing deterministic bands.
3. **Defer #142; take a higher-leverage issue.** The analytic `simultaneousCIs()` already
   covers the core need with an exact, deterministic band. Spend the slot on e.g. **#5**
   (unbalanced panels — the top applied-user blocker), **#30** (pre-trends diagnostic), or
   **#33** (`predict.fetwfe`, which already has a pushed branch). Optionally annotate #142
   with this assessment and leave it open as "nice-to-have, after higher-leverage work."

## Recommendation

Lean **Option 3 (defer)** unless one of the bootstrap's real advantages is actually wanted:
- if **large-K event-study families** are a near-term use case → Option 1 is justified
  (the `qmvnorm` scaling limit is the one place the analytic path genuinely strains);
- if the goal is purely **reviewer optics / `did`-parity** → Option 1 (faithful) is worth
  more than Option 2 (the lighter version invites the question "why is this random when the
  analytic band is exact?").

Option 2 is hard to recommend: it spends effort to add Monte-Carlo noise to an exact result.

## Open questions (would change the recommendation)

- Does the JoE reviewer request specifically demand a **bootstrap**, or just **simultaneous
  bands** (which already ship)? If the latter, #142 may be substantially satisfied already.
- Is there a concrete near-term analysis with a **large** effect family where `qmvnorm` is
  too slow/inaccurate? That's the strongest single reason to build it.
- Paper-side (`paper-same-data-gaussianity/SPEC.md`): is the multiplier bootstrap going into
  the paper as a method, or only mentioned as an "optional follow-up"? If it's load-bearing
  in the paper, parity argues for Option 1.

## Appendix — if we do Option 1, the build order

1. `.build_unit_if_matrix(result, psi_tes_mat, ...)` → `N_units x K` matrix `F`
   (regression IFs via the generalized `debiasedATT` score; propensity IFs via the `Sigma_2`
   Jacobians). Validate `crossprod(F)/N_units ≈ Sigma` (the analytic covariance) as a
   correctness check — they must agree up to the homoskedastic-vs-empirical difference.
2. `.multiplier_bootstrap(F, v_N, alpha, B, weights = c("rademacher","mammen"), seed)` →
   `c_hat`.
3. Thread into `simultaneousCIs(result, ..., method = c("analytic","bootstrap"), B, seed)`;
   default stays `"analytic"` (deterministic, back-compatible).
4. Tests: `crossprod(F)/N ≈ Sigma`; coverage on a simulated DGP; `c_hat` in
   `[1.96, Bonferroni]`; determinism given `seed`; agreement with the analytic `c` up to
   Monte-Carlo error at large `B`.
5. `inference_vignette.Rmd`: a section comparing pointwise / Bonferroni / analytic-sup-t /
   bootstrap-sup-t on one simulated event-study family.
