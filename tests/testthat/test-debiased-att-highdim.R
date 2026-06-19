# Tests for the high-dimensional (p >= NT) regime of debiasedATT() (#31).
#
# The fixed-p path builds the debiasing direction by the exact/ridged inverse
# v = (Sig + tiny*I)^{-1} a, which is invalid once p >= NT (Sig = X'X/n is
# singular). The high-dim path instead uses the nodewise (desparsified-lasso)
# Riesz representer `riesz_lasso()` of paper Theorem `debiased.highdim.thm`.
# These tests pin (a) the nodewise solver's fidelity and KKT/boundary behavior
# and (b) the accessor end-to-end on a genuine p >= NT fit.

# ---- An independent, naive reference for riesz_lasso ------------------------
# Same coordinate-descent minimizer of (1/2) v'Sig v - a'v + lambda||v||_1, but
# WITHOUT the package's O(p^2) rank-1 maintenance of `Sig %*% v`: it recomputes
# (Sig v)_j by a full inner product each coordinate. The two must agree to
# machine precision (the problem is strictly convex when Sig is nonsingular, so
# the minimizer is unique) -- this validates exactly the rank-1 update, the one
# optimization that could harbor a bug.
.ref_riesz_naive <- function(Sig, a, lambda, max_iter = 5000L, tol = 1e-12) {
	p <- length(a)
	v <- numeric(p)
	dg <- diag(Sig)
	conv <- FALSE
	for (it in seq_len(max_iter)) {
		mc <- 0
		for (j in seq_len(p)) {
			if (dg[j] <= 0) {
				next
			}
			Sv_j <- sum(Sig[j, ] * v) # full recompute of (Sig v)_j
			cj <- Sv_j - dg[j] * v[j]
			rho <- a[j] - cj
			vj <- sign(rho) * max(0, abs(rho) - lambda) / dg[j]
			dv <- vj - v[j]
			if (dv != 0) {
				v[j] <- vj
				ad <- abs(dv)
				if (ad > mc) {
					mc <- ad
				}
			}
		}
		if (mc < tol) {
			conv <- TRUE
			break
		}
	}
	v
}

# A reasonably-conditioned nonsingular Gram for the solver-level tests.
.psd_gram <- function(p = 40L, n = 200L, seed = 11L) {
	set.seed(seed)
	M <- matrix(stats::rnorm(n * p), n, p)
	crossprod(M) / n
}

# ---- A genuine p >= NT fetwfe fit (the simulator cannot make one: it enforces
# N >= (G+1)(d+1), so build the raw panel directly and supply the noise
# variances to bypass REML). N=20, T=8, d=12 => p = 376 >> NT = 160. -----------
.make_highdim_panel <- function() {
	set.seed(3)
	N <- 20L
	T <- 8L
	d <- 12L
	cohort_of_unit <- c(rep(0L, 8), rep(2L, 4), rep(3L, 4), rep(4L, 4))
	eff <- c(`2` = 0.5, `3` = 1.75, `4` = 3.0)
	covs <- matrix(stats::rnorm(N * d), N, d)
	do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cohort_of_unit[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:T,
				treat = as.integer(g > 0 & (1:T) >= g)
			)
			for (j in 1:d) {
				df[[paste0("x", j)]] <- covs[i, j]
			}
			te <- if (g > 0) eff[[as.character(g)]] else 0
			df$y <- 0.3 *
				(1:T) /
				T +
				te * df$treat +
				0.2 * covs[i, 1] +
				stats::rnorm(T, 0, 0.4)
			df
		})
	)
}

.make_highdim_fit <- function() {
	fetwfe(
		pdata = .make_highdim_panel(),
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:12),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 0.16,
		sig_eps_c_sq = 0.1
	)
}

# gls = FALSE sibling (#307): the SAME panel, but NO supplied variances -- the fit
# that errors today (the REML guard at p >= N(T-1)) and is the whole point of this
# path. calc_ses = FALSE; the debiasedATT() cluster-robust sandwich SE needs no Omega.
.make_highdim_fit_nogls <- function() {
	fetwfe(
		pdata = .make_highdim_panel(),
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:12),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
}

# Build the high-dim fit once (the fit is the slow part); reuse across tests.
hd_fix <- .make_highdim_fit()
hd_fix_nogls <- .make_highdim_fit_nogls()

# ============================ solver-level tests =============================

test_that("riesz_lasso matches an independent naive implementation (machine precision)", {
	Sig <- .psd_gram(p = 40L)
	set.seed(101)
	a <- stats::rnorm(40) * 0.5
	lam <- 0.05
	v_pkg <- riesz_lasso(Sig, a, lam, max_iter = 5000L, tol = 1e-12)
	v_nai <- .ref_riesz_naive(Sig, a, lam, max_iter = 5000L, tol = 1e-12)
	expect_lt(max(abs(v_pkg - v_nai)), 1e-12)
	expect_true(attr(v_pkg, "converged"))
})

test_that("riesz_lasso reduces to the exact inverse as lambda -> 0 (fixed-p boundary)", {
	Sig <- .psd_gram(p = 40L)
	set.seed(102)
	a <- stats::rnorm(40) * 0.5
	v <- riesz_lasso(Sig, a, 1e-10, max_iter = 2e5L, tol = 1e-12)
	expect_lt(max(abs(v - solve(Sig, a))), 1e-7)
	expect_true(attr(v, "converged"))
})

test_that("riesz_lasso satisfies its KKT feasibility certificate when converged", {
	Sig <- .psd_gram(p = 30L, seed = 7L)
	set.seed(103)
	a <- stats::rnorm(30)
	lam <- 0.1
	v <- riesz_lasso(Sig, a, lam, max_iter = 5000L, tol = 1e-12)
	# KKT for the l1-penalized objective: ||Sig v - a||_inf <= lambda.
	expect_true(attr(v, "converged"))
	expect_lte(attr(v, "feasibility"), lam * (1 + 1e-9))
	# attr matches a direct recomputation of the residual sup-norm.
	expect_equal(attr(v, "feasibility"), max(abs(Sig %*% v - a)))
})

test_that("riesz_lasso leaves all-zero (dead) columns at zero without NaNs", {
	Sig <- .psd_gram(p = 20L, seed = 5L)
	Sig[5, ] <- 0
	Sig[, 5] <- 0
	set.seed(104)
	a <- stats::rnorm(20)
	v <- riesz_lasso(Sig, a, 0.05)
	expect_identical(v[5], 0)
	expect_true(all(is.finite(v)))
})

test_that("riesz_lasso is deterministic", {
	Sig <- .psd_gram(p = 25L, seed = 9L)
	set.seed(105)
	a <- stats::rnorm(25)
	v1 <- riesz_lasso(Sig, a, 0.07)
	v2 <- riesz_lasso(Sig, a, 0.07)
	expect_identical(as.numeric(v1), as.numeric(v2))
	expect_identical(attr(v1, "feasibility"), attr(v2, "feasibility"))
})

test_that("lambda_node_default is the theory-scaled rate c * scale * sqrt(log p / N)", {
	expect_equal(lambda_node_default(p = 100, N = 25), sqrt(log(100) / 25))
	expect_equal(
		lambda_node_default(p = 376, N = 20, c = 2, scale = 3),
		2 * 3 * sqrt(log(376) / 20)
	)
})

# ====================== accessor end-to-end (p >= NT) ========================

test_that("the high-dim fixture really is p >= NT and fits with q < 1", {
	X <- hd_fix$internal$X_final
	expect_gte(ncol(X), nrow(X)) # p >= NT
	expect_true(isTRUE(hd_fix$internal$calc_ses))
})

test_that("debiasedATT runs in the p >= NT regime and returns nodewise diagnostics", {
	db <- debiasedATT(hd_fix)
	# the 6 fixed-p elements PLUS the high-dim diagnostics (feasibility / converged
	# / lambda_node + the #295 lambda_c / lambda_c_selection; lambda_cv only under
	# lambda_c = "cv").
	expect_identical(
		names(db),
		c(
			"att",
			"se",
			"ci_low",
			"ci_high",
			"var_reg",
			"var_weight",
			"feasibility",
			"converged",
			"lambda_node",
			"lambda_c",
			"lambda_c_selection"
		)
	)
	expect_identical(db$lambda_c_selection, "fixed") # default is the fixed 1.0
	expect_true(is.finite(db$att))
	expect_gt(db$se, 0)
	expect_true(is.finite(db$se))
	expect_gt(db$var_reg, 0)
	expect_gt(db$var_weight, 0)
	expect_equal(db$se, sqrt(db$var_reg + db$var_weight))
	# Wald interval ordering.
	expect_lt(db$ci_low, db$att)
	expect_lt(db$att, db$ci_high)
	# Debiased estimate differs from the fused plug-in (it is a different
	# estimator, not a re-report of att_hat).
	expect_false(isTRUE(all.equal(db$att, hd_fix$att_hat)))
})

test_that("the nodewise direction meets its feasibility certificate end-to-end", {
	db <- debiasedATT(hd_fix)
	expect_true(db$converged)
	# KKT: ||Sig v - a||_inf <= lambda_node (the relaxed-inverse bound Theorem
	# `debiased.highdim.thm` uses). With converged == TRUE this must hold.
	expect_lte(db$feasibility, db$lambda_node * (1 + 1e-9))
	expect_gt(db$lambda_node, 0)
})

test_that("the high-dim branch is load-bearing: the fixed-p inverse explodes on the singular Gram", {
	# Reconstruct Sig and the theta-space target a exactly as the accessor does,
	# then compare the nodewise direction against the WRONG fixed-p formula
	# solve(Sig + tiny*I, a). On the singular high-dim Gram the latter blows up;
	# the nodewise direction stays O(1). This is why the regime branch exists.
	X <- hd_fix$internal$X_final
	n <- nrow(X)
	p <- ncol(X)
	Sig <- crossprod(X) / n
	G <- hd_fix$G
	Tt <- hd_fix$T
	dd <- hd_fix$d
	ti <- hd_fix$treat_inds
	cp <- hd_fix$cohort_probs
	num_treats <- length(ti)
	fi <- getFirstInds(G = G, T = Tt)
	cot <- rep(seq_len(G), times = (Tt - 1):(Tt - G))
	a_beta <- numeric(p)
	for (g in seq_len(G)) {
		idx <- ti[cot == g]
		a_beta[idx] <- cp[g] / length(idx)
	}
	A <- genFullInvFusionTransformMat(
		first_inds = fi,
		T = Tt,
		G = G,
		d = dd,
		num_treats = num_treats,
		fusion_structure = hd_fix$fusion_structure,
		d_inv_treat = hd_fix$internal$d_inv_treat
	)
	a_th <- as.numeric(crossprod(A, a_beta))

	lam_node <- lambda_node_default(
		p = p,
		N = n / Tt,
		c = 1.0,
		scale = max(abs(a_th))
	)
	v_node <- riesz_lasso(Sig, a_th, lam_node)
	v_fixedp <- solve(Sig + (1e-6 * mean(diag(Sig))) * diag(p), a_th)

	expect_lt(max(abs(v_node)), 50) # nodewise: O(1)
	expect_gt(max(abs(v_fixedp)), 1e4) # fixed-p inverse: blows up
	expect_gt(max(abs(v_fixedp)) / max(abs(v_node)), 100)

	# End-to-end pin: the accessor's lambda_node and debiased att must equal an
	# independent reconstruction from the same a_th / v_node above (the high-dim
	# analog of the fixed-p reference test). This pins the scale = max(|a|), the
	# N = clusters scaling, and the score/correction formula through the accessor.
	# The high-dim nuisance is the internal q=1 fused lasso (#303), NOT the bridge
	# theta_hat -- reconstruct it the same way the accessor does (n / Tt units), so
	# this stays an INDEPENDENT reconstruction (the deterministic seed gives the
	# identical nuisance => exact match), not a tautology.
	th_q1 <- fetwfe:::.fit_q1_nuisance(
		X,
		as.numeric(hd_fix$internal$y_final),
		n / Tt,
		Tt
	)
	resid <- as.numeric(hd_fix$internal$y_final - th_q1[1] - X %*% th_q1[-1])
	att_recompute <- sum(c(0, a_th) * th_q1) + mean((X %*% v_node) * resid)
	db <- debiasedATT(hd_fix)
	expect_equal(db$lambda_node, lam_node)
	expect_equal(db$att, att_recompute)
})

test_that("lambda_c scales the nodewise penalty and shrinks the debiasing direction", {
	db1 <- debiasedATT(hd_fix, lambda_c = 1)
	db2 <- debiasedATT(hd_fix, lambda_c = 2)
	# lambda_node is exactly linear in lambda_c (everything else held fixed).
	expect_equal(db2$lambda_node, 2 * db1$lambda_node)
	# At lambda_c = 2 the penalty lambda_node exceeds ||a||_inf, so v = 0 is optimal
	# and the debiased estimate collapses to the *q = 1 plug-in* (#303) -- NOT
	# att_hat (the bridge plug-in). A much larger lambda_c gives the same v = 0, so
	# both collapse to the identical q = 1 plug-in; and that plug-in differs from
	# the bridge att_hat -- the mutation signal that the nuisance actually moved.
	expect_equal(db2$att, debiasedATT(hd_fix, lambda_c = 100)$att)
	expect_false(isTRUE(all.equal(db2$att, hd_fix$att_hat)))
	# At the default lambda_c = 1 the correction is active (att != the plug-in).
	expect_false(isTRUE(all.equal(db1$att, db2$att)))
	expect_false(isTRUE(all.equal(db1$att, hd_fix$att_hat)))
})

test_that("debiasedATT warns when the nodewise solver does not fully resolve", {
	# Tiny penalty + capped iterations => the KKT feasibility constraint cannot be
	# met within the budget => the feasibility warning fires.
	expect_warning(
		debiasedATT(hd_fix, lambda_c = 1e-6, riesz_max_iter = 200L),
		"feasibility"
	)
	# A single sweep is feasible (loose penalty) but not converged => the softer
	# "more iterations" warning fires instead.
	expect_warning(
		debiasedATT(hd_fix, riesz_max_iter = 1L),
		"more iterations"
	)
})

test_that("debiasedATT is deterministic in the high-dim regime", {
	a <- debiasedATT(hd_fix)
	b <- debiasedATT(hd_fix)
	expect_identical(a$att, b$att)
	expect_identical(a$se, b$se)
	expect_identical(a$lambda_node, b$lambda_node)
})

test_that("debiasedATT validates the high-dim solver arguments", {
	expect_error(debiasedATT(hd_fix, lambda_c = 0), "lambda_c")
	expect_error(debiasedATT(hd_fix, lambda_c = -1), "lambda_c")
	expect_error(debiasedATT(hd_fix, lambda_c = c(1, 2)), "lambda_c")
	expect_error(debiasedATT(hd_fix, riesz_max_iter = 0), "riesz_max_iter")
	expect_error(debiasedATT(hd_fix, riesz_tol = 0), "riesz_tol")
	expect_error(debiasedATT(hd_fix, riesz_tol = -1e-9), "riesz_tol")
})

# ==============================================================================
# #307: Omega-free high-dim cluster-robust SE (gls = FALSE). The headline -- a
# p >= NT fit WITHOUT supplied variances (impossible today: the REML guard) --
# plus the load-bearing plug-in V2 == att_var_2 cross-check.
# ==============================================================================

test_that("gls = FALSE fits a p >= NT panel without supplied variances (#307)", {
	expect_s3_class(hd_fix_nogls, "fetwfe")
	# genuinely high-dim, and no whitening / variance components were estimated.
	expect_gte(
		ncol(hd_fix_nogls$internal$X_final),
		nrow(hd_fix_nogls$internal$X_final)
	)
	expect_false(isTRUE(hd_fix_nogls$internal$calc_ses))
	expect_true(is.na(hd_fix_nogls$sig_eps_sq))
	expect_true(is.na(hd_fix_nogls$sig_eps_c_sq))
	# the oracle SE machinery did not run, so att_var_2 is absent.
	av2 <- hd_fix_nogls$internal$variance_components$att_var_2
	expect_true(is.null(av2) || is.na(av2))
	# the headline ATT is still produced.
	expect_true(is.finite(hd_fix_nogls$att_hat))
})

test_that("debiasedATT() returns a finite cluster-robust SE on a gls = FALSE fit (#307)", {
	db <- debiasedATT(hd_fix_nogls)
	expect_true(all(
		c(
			"att",
			"se",
			"ci_low",
			"ci_high",
			"var_reg",
			"var_weight",
			"feasibility",
			"converged",
			"lambda_node"
		) %in%
			names(db)
	))
	expect_true(is.finite(db$att))
	expect_true(is.finite(db$se) && db$se > 0)
	expect_true(is.finite(db$var_reg) && db$var_reg > 0) # V1 unit-clustered sandwich
	expect_true(is.finite(db$var_weight) && db$var_weight > 0) # V2 plug-in
	expect_equal(db$se, sqrt(db$var_reg + db$var_weight))
	expect_equal(db$ci_low, db$att - stats::qnorm(0.975) * db$se)
})

test_that("debiasedATT() is deterministic on a gls = FALSE fit (#307)", {
	a <- debiasedATT(hd_fix_nogls)
	b <- debiasedATT(hd_fix_nogls)
	expect_identical(a$att, b$att)
	expect_identical(a$se, b$se)
	expect_identical(a$lambda_node, b$lambda_node)
})

# The load-bearing validation: the plug-in V2 (`.plugin_v2`, used when att_var_2 is
# absent) EQUALS the oracle `att_var_2` the whitened path uses. hd_fix has equal
# treated cohorts; the inline fit below has UNEQUAL cohorts (4 / 6 / 8 treated
# units), so the match is not a balanced-design coincidence. Mutation-checkable:
# a wrong N_tau (e.g. N instead of N_treated) or dropping the cohort_probs
# weighting in `.plugin_v2` breaks this past 1e-10.
test_that(".plugin_v2 equals att_var_2 (equal, unequal, two-sample cohorts) (#307)", {
	expect_equal(
		fetwfe:::.plugin_v2(hd_fix),
		hd_fix$internal$variance_components$att_var_2,
		tolerance = 1e-10
	)
	set.seed(7)
	Nq <- 24L
	Tq <- 6L
	cohort_of_unit <- c(rep(0L, 6), rep(3L, 4), rep(4L, 6), rep(5L, 8))
	eff <- c(`3` = 1.0, `4` = 2.0, `5` = -1.0)
	cv <- stats::rnorm(Nq)
	rows <- do.call(
		rbind,
		lapply(seq_len(Nq), function(i) {
			g <- cohort_of_unit[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:Tq,
				treat = as.integer(g > 0 & (1:Tq) >= g),
				x1 = cv[i]
			)
			te <- if (g > 0) eff[[as.character(g)]] else 0
			df$y <- 0.1 *
				(1:Tq) +
				te * df$treat +
				0.3 * cv[i] +
				stats::rnorm(Tq, 0, 0.5)
			df
		})
	)
	fit_uneq <- fetwfe(
		pdata = rows,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = "x1",
		response = "y",
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 0.25,
		sig_eps_c_sq = 0.1
	)
	# confirm the treated cohorts really are unequal (so the test bites weighting).
	expect_gt(max(fit_uneq$cohort_probs) - min(fit_uneq$cohort_probs), 0.05)
	expect_equal(
		fetwfe:::.plugin_v2(fit_uneq),
		fit_uneq$internal$variance_components$att_var_2,
		tolerance = 1e-10
	)
	# two-sample (indep_counts): the plug-in must match the TWO-SAMPLE att_var_2
	# (the case the naive single-sample (1/N_T) recompute gets wrong). simulateData()
	# auto-populates indep_counts; fetwfeWithSimulatedData() keeps it.
	coefs_ic <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.6,
		eff_size = 1.5,
		seed = 7
	)
	dat_ic <- simulateData(
		coefs_ic,
		N = 80,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 7
	)
	fit_indep <- fetwfeWithSimulatedData(dat_ic, q = 0.5)
	expect_true(isTRUE(fit_indep$indep_counts_used))
	expect_equal(
		fetwfe:::.plugin_v2(fit_indep),
		fit_indep$internal$variance_components$att_var_2,
		tolerance = 1e-10
	)
})

test_that("the att_var_2 path is byte-unchanged: supplied-var var_weight == att_var_2 (#307)", {
	# the plug-in fallback must NOT perturb the whitened / supplied-variance path.
	db <- debiasedATT(hd_fix)
	expect_equal(db$var_weight, hd_fix$internal$variance_components$att_var_2)
})

# A gls = FALSE FIXED-p (p < NT) fit must error naming the REAL cause (the fixed-p
# cluster-robust path is the #312 follow-up), not the misleading "requires q < 1"
# (the fit IS q < 1). This is the PR #311 review's boundary-gap clarity ask.
test_that("a gls = FALSE fixed-p fit errors with the real cause, not 'requires q < 1' (#307/#312)", {
	set.seed(5)
	N <- 30L
	Tt <- 5L
	cohort_of_unit <- c(rep(0L, 12), rep(3L, 9), rep(4L, 9))
	eff <- c(`3` = 1, `4` = 2)
	cv <- stats::rnorm(N)
	rows <- do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cohort_of_unit[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:Tt,
				treat = as.integer(g > 0 & (1:Tt) >= g),
				x1 = cv[i]
			)
			te <- if (g > 0) eff[[as.character(g)]] else 0
			df$y <- 0.1 *
				(1:Tt) +
				te * df$treat +
				0.3 * cv[i] +
				stats::rnorm(Tt, 0, 0.5)
			df
		})
	)
	fit <- fetwfe(
		pdata = rows,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = "x1",
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
	expect_lt(ncol(fit$internal$X_final), nrow(fit$internal$X_final)) # fixed-p
	# the message names gls = FALSE + the #312 follow-up (NOT "requires q < 1"; the
	# "bridge selection" q-attribution lives in the other branch only).
	expect_error(debiasedATT(fit), "gls = FALSE")
	expect_error(debiasedATT(fit), "follow-up")
})
