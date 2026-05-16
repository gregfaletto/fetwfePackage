# Tests for the BETWFE partial-selection SE fix.
#
# Background: prior to v1.9.2, getPsiRUnfused() normalized psi_r by k_sel
# (count of *selected* coefficients in cohort r's treatment block) rather
# than k_full (the full block size). The cohort point estimate
# cohort_tes[r] = mean(tes[first_ind_r:last_ind_r]) averages over the full
# block (unselected entries are exact zeros post-bridge), so the variance
# formula targeted the wrong estimand whenever the bridge solver zeroed
# out some-but-not-all of a cohort's coefficients. Per-cohort SEs were
# inflated by factor k_full / k_sel; overall att_se shifted (direction
# depends on cohort-dispersion structure). See
# .plans/follow-ups/issue-draft-betwfe-se-partial-selection.md.

test_that("getPsiRUnfused invariant: t(psi_r) %*% theta_sel == cohort_tes[r]", {
	# Constructed scenario: cohort 1 occupies indices 1..4 (k_full = 4); the
	# bridge selected positions 1, 2, 5 across all cohorts, so within cohort
	# 1 only positions 1 and 2 are selected (k_sel = 2). theta_hat_treat_sel
	# has 3 entries (length of sel_treat_inds_shifted).
	psi_r <- fetwfe:::getPsiRUnfused(
		first_ind_r = 1,
		last_ind_r = 4,
		sel_treat_inds_shifted = c(1L, 2L, 5L),
		gram_inv = diag(3)
	)
	# Post-fix: psi_r = c(1/4, 1/4, 0).
	expect_equal(psi_r, c(0.25, 0.25, 0))

	# Invariant: psi_r %*% theta_sel matches cohort_tes[r] for a full
	# cohort-block tes vector (zeros at unselected positions).
	theta_sel <- c(1.0, 2.0, 3.0) # the bridge-selected coefficient values
	# cohort_tes[1] = mean(tes[1:4]) where tes = c(1, 2, 0, 0)
	# (positions 3-4 unselected → 0).
	cohort_te_expected <- (1.0 + 2.0 + 0 + 0) / 4
	expect_equal(as.numeric(psi_r %*% theta_sel), cohort_te_expected)
})

test_that("getPsiRUnfused is bit-identical to pre-fix on full-selection (ETWFE/twfeCovs path)", {
	# When sel_treat_inds_shifted = 1:num_treats (the ETWFE / twfeCovs call
	# shape), every index is "selected", so k_sel = k_full. The post-fix
	# psi_r[inds_r] <- 1 / k_full equals the pre-fix
	# psi_r[inds_r] <- 1; psi_r <- psi_r / sum(psi_r) = 1 / k_full.
	num_treats <- 7L
	psi_r <- fetwfe:::getPsiRUnfused(
		first_ind_r = 1,
		last_ind_r = 4,
		sel_treat_inds_shifted = seq_len(num_treats),
		gram_inv = diag(num_treats)
	)
	# Cohort 1 occupies indices 1..4 (k_full=4), all selected.
	expect_equal(psi_r, c(0.25, 0.25, 0.25, 0.25, 0, 0, 0))
})

test_that("BETWFE SE/SD ratio is ~1 under partial selection (MC)", {
	skip_on_cran()

	# A DGP where cohort 1 reliably has partial selection: cohort_betas
	# of c(4, 4, 0, 0) so bridge selects positions 1-2 but zeros 3-4 at
	# this signal-to-noise level. Two cohorts (R=2), short panel (T=5).
	# `density = 0.5` controls the sparsity in genCoefs internally; the
	# `eff_size = 4` and `sig_eps_sq = 1` combination gives a clear SNR.
	#
	# Critical: simulateData(coefs_obj) reads coefs_obj$seed and calls
	# set.seed() internally. Wrapping the loop in set.seed(i) does NOT
	# vary noise across reps. Must mutate coefs_obj$seed.
	coefs_template <- genCoefs(
		R = 2,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 4,
		seed = 42
	)

	n_reps <- 30
	catt_hat_1 <- numeric(n_reps)
	catt_se_1 <- numeric(n_reps)
	att_hat <- numeric(n_reps)
	att_se <- numeric(n_reps)

	for (i in seq_len(n_reps)) {
		coefs_i <- coefs_template
		coefs_i$seed <- 8000L + i
		sim <- simulateData(
			coefs_i,
			N = 300,
			sig_eps_sq = 1,
			sig_eps_c_sq = 1
		)
		res <- betwfeWithSimulatedData(sim)
		catt_hat_1[i] <- res$catt_df[["Estimated TE"]][1]
		catt_se_1[i] <- res$catt_df$SE[1]
		att_hat[i] <- res$att_hat
		att_se[i] <- res$att_se
	}

	# Per-cohort SE/SD ratio should be ~1 post-fix (pre-fix it was ~2).
	# Allow generous slack against MC noise on 30 reps.
	cohort1_ratio <- mean(catt_se_1, na.rm = TRUE) / sd(catt_hat_1)
	expect_gt(cohort1_ratio, 0.7)
	expect_lt(cohort1_ratio, 1.5)

	# Overall att_se/sd ratio is also informative but the direction effect
	# is more complex (att_var_2 can grow under the fix). Sanity bounds.
	att_ratio <- mean(att_se, na.rm = TRUE) / sd(att_hat)
	expect_gt(att_ratio, 0.5)
	expect_lt(att_ratio, 3.0)
})
