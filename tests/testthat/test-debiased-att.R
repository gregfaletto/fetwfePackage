# Tests for debiasedATT() / debiasedATTWithSimulatedData() (#291).
#
# The pivotal test reproduces the validated paper-repo reference algorithm
# `simulations/method_functions.R::debiased_fetwfe` (inlined here, since the
# paper repo is not a package dependency) and asserts the accessor matches it to
# < 1e-10 on a cohort fit. (Verified bite-y: perturbing the ridge constant, the
# N_T definition, or the per-unit cluster sums in R/debiased_att.R blows this
# comparison far past 1e-10.)

# Inlined reference: the exact arithmetic of `debiased_fetwfe`, hard-coding
# fusion_structure = "cohort" and taking N_T from the simulated pdata (as the
# reference does). Valid only for cohort fits (the reference's hard-coded
# "cohort" transform); the accessor generalizes it by threading fit$fusion_structure.
.ref_debiased_fetwfe <- function(fit, dat) {
	G <- fit$G
	Tt <- fit$T
	d <- fit$d
	ti <- fit$treat_inds
	cp <- fit$cohort_probs
	num_treats <- length(ti)
	p <- length(fit$beta_hat)
	fi <- getFirstInds(G = G, T = Tt)
	sizes <- (Tt - 1):(Tt - G)
	cohort_of_treat <- rep(seq_len(G), times = sizes)
	a_beta <- numeric(p)
	for (g in seq_len(G)) {
		idx <- ti[cohort_of_treat == g]
		a_beta[idx] <- cp[g] / length(idx)
	}
	A <- genFullInvFusionTransformMat(
		first_inds = fi,
		T = Tt,
		G = G,
		d = d,
		num_treats = num_treats,
		fusion_structure = "cohort",
		d_inv_treat = NULL
	)
	a_theta <- c(0, as.numeric(crossprod(A, a_beta)))
	X <- fit$internal$X_final
	y <- fit$internal$y_final
	th <- fit$internal$theta_hat
	n <- nrow(X)
	Sig <- crossprod(X) / n
	v <- solve(Sig + (1e-6 * mean(diag(Sig))) * diag(ncol(Sig)), a_theta[-1])
	r <- as.numeric(y - th[1] - X %*% th[-1])
	g <- as.numeric((X %*% v) * r)
	tau <- sum(a_theta * th) + mean(g)
	unit <- rep(seq_len(n / Tt), each = Tt)
	Gi <- tapply(g, unit, sum)
	se_theta_sq <- sum(Gi^2) / n^2
	pd <- dat$pdata
	N_T <- length(unique(pd[[dat$unit_var]][pd[[dat$treatment]] == 1]))
	var_weight <- sum(cp * (fit$catt_hats - fit$att_hat)^2) / N_T
	c(att = tau, se = sqrt(se_theta_sq + var_weight))
}

.make_fit <- function(fusion_structure = "cohort", q = 0.5, seed = 7) {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.6,
		eff_size = 1.5,
		seed = seed,
		fusion_structure = fusion_structure
	)
	dat <- simulateData(
		coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = seed
	)
	# Single-sample fit (drop the auto-populated indep_counts): the debiased SE is
	# validated for the single-sample case, and this avoids the two-sample warning.
	dat$indep_counts <- NA
	list(
		fit = fetwfeWithSimulatedData(
			dat,
			q = q,
			fusion_structure = fusion_structure
		),
		dat = dat
	)
}

test_that("debiasedATT matches the validated reference algorithm to < 1e-10", {
	f <- .make_fit("cohort")
	db <- debiasedATT(f$fit)
	ref <- .ref_debiased_fetwfe(f$fit, f$dat)
	expect_lt(abs(db$att - ref[["att"]]), 1e-10)
	expect_lt(abs(db$se - ref[["se"]]), 1e-10)
})

test_that("debiasedATT return shape and Wald interval", {
	f <- .make_fit("cohort")
	db <- debiasedATT(f$fit, alpha = 0.05)
	expect_identical(
		names(db),
		c("att", "se", "ci_low", "ci_high", "var_reg", "var_weight")
	)
	expect_true(all(vapply(
		db,
		function(x) is.numeric(x) && length(x) == 1L,
		NA
	)))
	# SE is the sum of the two channels; CI is the Wald interval.
	expect_equal(db$se, sqrt(db$var_reg + db$var_weight))
	z <- stats::qnorm(0.975)
	expect_equal(db$ci_low, db$att - z * db$se)
	expect_equal(db$ci_high, db$att + z * db$se)
	expect_gt(db$var_reg, 0)
	expect_gte(db$var_weight, 0)
})

test_that("debiased estimate differs from the fused att_hat", {
	f <- .make_fit("cohort")
	db <- debiasedATT(f$fit)
	# By construction the debiased (ETWFE-direction) estimate is not the fused one.
	expect_false(isTRUE(all.equal(db$att, f$fit$att_hat)))
})

test_that("debiasedATT generalizes to event-study fits (fusion_structure threaded)", {
	# The reference hard-codes the 'cohort' transform; an event-study fit needs the
	# fit's own fusion_structure or the ATT-direction identity fails. The accessor
	# threads fit$fusion_structure, so its internal identity guard passes and it
	# returns a finite, well-formed result.
	f <- .make_fit("event_study")
	expect_identical(f$fit$fusion_structure, "event_study")
	db <- expect_no_error(debiasedATT(f$fit))
	expect_true(is.finite(db$att) && is.finite(db$se) && db$se > 0)
	expect_equal(db$se, sqrt(db$var_reg + db$var_weight))
	# The reference's hard-coded 'cohort' transform would violate the identity on
	# this fit (so the accessor's threading is load-bearing, not cosmetic).
	G <- f$fit$G
	Tt <- f$fit$T
	d <- f$fit$d
	ti <- f$fit$treat_inds
	cp <- f$fit$cohort_probs
	p <- length(f$fit$beta_hat)
	fi <- getFirstInds(G = G, T = Tt)
	cot <- rep(seq_len(G), times = (Tt - 1):(Tt - G))
	a_beta <- numeric(p)
	for (g in seq_len(G)) {
		idx <- ti[cot == g]
		a_beta[idx] <- cp[g] / length(idx)
	}
	A_wrong <- genFullInvFusionTransformMat(
		first_inds = fi,
		T = Tt,
		G = G,
		d = d,
		num_treats = length(ti),
		fusion_structure = "cohort",
		d_inv_treat = NULL
	)
	a_theta_wrong <- c(0, as.numeric(crossprod(A_wrong, a_beta)))
	expect_gt(
		abs(sum(a_theta_wrong * f$fit$internal$theta_hat) - f$fit$att_hat),
		1e-6
	)
})

test_that("debiasedATT handles a custom fusion_matrix fit (#236)", {
	# A custom fusion_matrix fit has fusion_structure == "cohort" but a non-NULL
	# internal$d_inv_treat; the accessor must thread that (not hard-code the
	# built-in cohort transform), or the ATT-direction identity fails.
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.6,
		eff_size = 1.5,
		seed = 7
	)
	dat <- simulateData(
		coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 7
	)
	dat$indep_counts <- NA # single-sample (no two-sample warning)
	nt <- getNumTreats(G = 3, T = 5)
	D <- diag(nt)
	D[lower.tri(D)] <- 1 # finite, invertible num_treats x num_treats
	fit_custom <- fetwfeWithSimulatedData(dat, q = 0.5, fusion_matrix = D)
	expect_false(is.null(fit_custom$internal$d_inv_treat))
	db <- expect_no_error(debiasedATT(fit_custom))
	expect_true(is.finite(db$att) && is.finite(db$se) && db$se > 0)
	expect_equal(db$se, sqrt(db$var_reg + db$var_weight))
})

test_that("debiasedATTWithSimulatedData matches debiasedATT on the same fit", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.6,
		eff_size = 1.5,
		seed = 11
	)
	dat <- simulateData(
		coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 11
	)
	dat$indep_counts <- NA # single-sample (no two-sample warning)
	w <- debiasedATTWithSimulatedData(dat, q = 0.5)
	db <- debiasedATT(fetwfeWithSimulatedData(dat, q = 0.5))
	expect_equal(w$att, db$att)
	expect_equal(w$se, db$se)
})

test_that("debiasedATT warns on indep_counts fits, not on single-sample fits", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.6,
		eff_size = 1.5,
		seed = 7
	)
	dat <- simulateData(
		coefs,
		N = 150,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 7
	)
	# indep_counts fit (sim default) -> the cohort-weight channel is not validated
	# for the two-sample case, so the accessor warns.
	fit_indep <- fetwfeWithSimulatedData(dat, q = 0.5)
	expect_true(isTRUE(fit_indep$indep_counts_used))
	expect_warning(debiasedATT(fit_indep), "indep_counts")

	# Single-sample fit -> no warning.
	dat_single <- dat
	dat_single$indep_counts <- NA
	fit_single <- fetwfeWithSimulatedData(dat_single, q = 0.5)
	expect_false(isTRUE(fit_single$indep_counts_used))
	expect_no_warning(debiasedATT(fit_single))
})

test_that("debiasedATT input guards", {
	f <- .make_fit("cohort")

	# q >= 1 (no valid SE / calc_ses FALSE)
	fit_q2 <- fetwfeWithSimulatedData(f$dat, q = 2)
	expect_error(debiasedATT(fit_q2), "q < 1")

	# add_ridge = TRUE (ridge-augmented design) is unsupported
	fit_ridge <- fetwfeWithSimulatedData(f$dat, q = 0.5, add_ridge = TRUE)
	expect_error(debiasedATT(fit_ridge), "add_ridge")

	# non-fetwfe estimators are unsupported
	expect_error(
		debiasedATT(etwfeWithSimulatedData(f$dat)),
		"requires a fitted"
	)
	expect_error(debiasedATT(list(a = 1)), "requires a fitted")

	# alpha out of (0, 1)
	expect_error(debiasedATT(f$fit, alpha = 1.5), "alpha")
	expect_error(debiasedATT(f$fit, alpha = 0), "alpha")

	# p >= NT high-dimensional regime (synthetic: shrink BOTH X_final and y_final
	# to fewer rows than columns so the p >= NT guard fires -- keeping their row
	# counts equal so the add_ridge guard does not pre-empt it)
	bad <- f$fit
	pcols <- ncol(bad$internal$X_final)
	bad$internal$X_final <- bad$internal$X_final[
		seq_len(pcols - 1L),
		,
		drop = FALSE
	]
	bad$internal$y_final <- bad$internal$y_final[seq_len(pcols - 1L)]
	expect_error(debiasedATT(bad), "high-dimensional")
})
