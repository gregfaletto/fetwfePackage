# Tests for #312: Omega-free fixed-p (p < NT) cluster-robust debiasedATT() SE on
# gls = FALSE fits. #307 made the high-dim (p >= NT) SE Omega-free; #312 extends
# the SAME to the fixed-p regime. The fixed-p var_reg is the exact-inverse OLS
# direction + unit-clustered residual sandwich (sig_eps_sq-free), and var_weight
# falls back to the propensity-only plug-in (.plugin_v2). The gate now accepts any
# un-whitened fit (`is.na(sig_eps_sq)`) regardless of q, rejecting only a WHITENED
# q >= 1 fit (`!is.na(sig_eps_sq)`).

.mk_fixedp_fit_312 <- function(gls, q = 0.5) {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 60,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	args <- list(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = q,
		verbose = FALSE,
		gls = gls
	)
	if (gls) {
		args$sig_eps_sq <- 1
		args$sig_eps_c_sq <- 0.5
	}
	do.call(fetwfe, args)
}

test_that("the fixed-p var_reg equals an independent unit-clustered sandwich reconstruction (#312)", {
	# Load-bearing correctness pin: rebuild the cluster-robust sandwich from scratch
	# (OLS theta_hat, exact-inverse direction, OLS residuals) -- NOT via the accessor.
	f <- .mk_fixedp_fit_312(gls = FALSE, q = 0.5)
	expect_lt(ncol(f$internal$X_final), nrow(f$internal$X_final))
	db <- debiasedATT(f)
	X <- f$internal$X_final
	y <- as.numeric(f$internal$y_final)
	n <- nrow(X)
	p <- ncol(X)
	Tt <- f$T
	N <- f$N
	G <- f$G
	ti <- f$treat_inds
	A <- genFullInvFusionTransformMat(
		first_inds = getFirstInds(G = G, T = Tt),
		T = Tt,
		G = G,
		d = f$d,
		num_treats = length(ti),
		fusion_structure = f$fusion_structure,
		d_inv_treat = f$internal$d_inv_treat
	)
	sizes <- (Tt - 1):(Tt - G)
	cot <- rep(seq_len(G), times = sizes)
	a_beta <- numeric(p)
	for (g in seq_len(G)) {
		idx <- ti[cot == g]
		a_beta[idx] <- f$cohort_probs[g] / length(idx)
	}
	a_theta <- c(0, as.numeric(crossprod(A, a_beta)))
	Sig <- crossprod(X) / n
	v <- solve(Sig + 1e-6 * mean(diag(Sig)) * diag(p), a_theta[-1])
	th <- f$internal$theta_hat
	resid <- as.numeric(y - th[1] - X %*% th[-1])
	score <- as.numeric((X %*% v) * resid)
	unit <- rep(seq_len(N), each = Tt)
	var_reg_ref <- sum(tapply(score, unit, sum)^2) / n^2
	expect_equal(db$var_reg, var_reg_ref, tolerance = 1e-10)
	# var_weight is the Omega-free plug-in; both channels positive on this fixture.
	expect_identical(db$var_weight, fetwfe:::.plugin_v2(f))
	expect_gt(db$var_weight, 0)
	expect_equal(db$se, sqrt(db$var_reg + db$var_weight))
})

test_that("gls = FALSE and gls = TRUE fixed-p agree on the point estimate but differ in SE (whitening buys efficiency)", {
	# Structure check (NOT byte-equality): same target ATT, different var_reg
	# (un-whitened sandwich is the less efficient, still-valid estimate).
	f0 <- .mk_fixedp_fit_312(gls = FALSE, q = 0.5)
	f1 <- .mk_fixedp_fit_312(gls = TRUE, q = 0.5)
	d0 <- debiasedATT(f0)
	d1 <- debiasedATT(f1)
	expect_lt(abs(d0$att - d1$att), 0.25) # same overall-ATT target
	expect_false(isTRUE(all.equal(d0$var_reg, d1$var_reg))) # whitening changes var_reg
})

test_that("a gls = FALSE fixed-p fit with q >= 1 is accepted (the SE is q-independent; #312)", {
	# The cluster-robust sandwich + exact-inverse OLS identity are valid for ANY q;
	# is.na(sig_eps_sq) (un-whitened) is what the gate keys on, not q. Pins the
	# "accept regardless of q" decision so it can't be silently re-narrowed.
	f <- .mk_fixedp_fit_312(gls = FALSE, q = 1)
	expect_true(is.na(f$sig_eps_sq))
	db <- debiasedATT(f)
	expect_true(is.finite(db$att))
	expect_gt(db$se, 0)
})

test_that("a WHITENED (gls = TRUE) q >= 1 fixed-p fit still errors (the kept boundary)", {
	f <- .mk_fixedp_fit_312(gls = TRUE, q = 2)
	expect_false(isTRUE(f$internal$calc_ses))
	expect_false(is.na(f$sig_eps_sq)) # REML estimated it even at q >= 1
	expect_error(debiasedATT(f), "q < 1")
})

# ---- the boundary band N(T-1) <= p < NT (REML infeasible, gls = FALSE only) ----
.mk_band_panel_312 <- function() {
	set.seed(4)
	N <- 12L
	T <- 5L
	d <- 2L
	cohort_of_unit <- c(rep(0L, 3), rep(2L, 3), rep(3L, 3), rep(4L, 3))
	eff <- c(`2` = 1, `3` = 2.5, `4` = 4)
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
			df$y <- 0.2 *
				(1:T) +
				te * df$treat +
				0.3 * covs[i, 1] +
				stats::rnorm(T, 0, 0.5)
			df
		})
	)
}

test_that("the boundary band N(T-1) <= p < NT gets a finite SE (gls = TRUE infeasible there)", {
	panel <- .mk_band_panel_312()
	fb <- fetwfe(
		pdata = panel,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = c("x1", "x2"),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
	p <- ncol(fb$internal$X_final)
	n <- nrow(fb$internal$X_final)
	N <- fb$N
	Tt <- fb$T
	expect_gte(p, N * (Tt - 1)) # in the band ...
	expect_lt(p, n) # ... but still fixed-p
	# gls = TRUE is infeasible in the band (REML guard p < N(T-1)).
	expect_error(
		fetwfe(
			pdata = panel,
			time_var = "year",
			unit_var = "unit",
			treatment = "treat",
			covs = c("x1", "x2"),
			response = "y",
			q = 0.5,
			verbose = FALSE,
			gls = TRUE
		)
	)
	db <- debiasedATT(fb)
	expect_true(is.finite(db$att))
	expect_true(is.finite(db$se) && db$se > 0)
})
