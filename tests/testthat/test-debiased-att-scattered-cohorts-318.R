# Tests for #318: debiasedATT() supports non-consecutive (scattered-cohort)
# adoption. Its overall-ATT weight vector `a_theta` previously used
# `getFirstInds(G, T)` + `sizes <- (T-1):(T-G)`, both of which assume consecutive
# adoption (offsets 2..G+1); on scattered adoption (gaps, e.g. years 2,4,7 --
# bacondecomp::divorce) the treatment cells were mis-assigned and the ATT-identity
# guard correctly fired. The fix resolves the ACTUAL offsets via the #174
# dispatcher (`.resolve_cohort_offsets_and_first_inds`), so the identity
# `a_theta' theta_hat == att_hat` holds for scattered cohorts too, in both regimes.

# Raw panel with `adopt` adoption years (integer-coercible cohort labels, so the
# #174 dispatcher derives the real offsets rather than the consecutive fallback).
.mk_panel_318 <- function(adopt, d = 3L, N = 24L, T = 8L, seed = 5) {
	set.seed(seed)
	covs <- matrix(stats::rnorm(N * d), N, d)
	cohort_of_unit <- c(
		rep(0L, N - 3L * 6L),
		rep(adopt[1], 6),
		rep(adopt[2], 6),
		rep(adopt[3], 6)
	)
	eff <- stats::setNames(c(0.5, 2.0, 3.5), as.character(adopt))
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

.fit_318 <- function(panel, gls, ...) {
	fetwfe(
		pdata = panel,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = grep("^x", names(panel), value = TRUE),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = gls,
		...
	)
}

# Independent reconstruction of the overall-ATT identity check the guard performs.
.att_identity_gap <- function(fit) {
	G <- fit$G
	T <- fit$T
	nt <- length(fit$treat_inds)
	p <- ncol(fit$internal$X_final)
	fi <- fetwfe:::.resolve_cohort_offsets_and_first_inds(
		fit,
		G = G,
		T = T
	)$first_inds
	sizes <- diff(c(fi, nt + 1L))
	cot <- rep(seq_len(G), times = sizes)
	a_beta <- numeric(p)
	for (g in seq_len(G)) {
		idx <- fit$treat_inds[cot == g]
		a_beta[idx] <- fit$cohort_probs[g] / length(idx)
	}
	A <- genFullInvFusionTransformMat(
		first_inds = fi,
		T = T,
		G = G,
		d = fit$d,
		num_treats = nt,
		fusion_structure = fit$fusion_structure,
		d_inv_treat = fit$internal$d_inv_treat
	)
	a_theta <- c(0, as.numeric(crossprod(A, a_beta)))
	abs(sum(a_theta * fit$internal$theta_hat) - fit$att_hat)
}

test_that("debiasedATT() runs on a scattered-cohort FIXED-p fit (#318)", {
	f <- .fit_318(
		.mk_panel_318(c(2L, 4L, 7L), d = 3L),
		gls = TRUE,
		sig_eps_sq = 0.16,
		sig_eps_c_sq = 0.1
	)
	expect_lt(ncol(f$internal$X_final), nrow(f$internal$X_final)) # fixed-p
	# the resolved offsets are the SCATTERED ones, not consecutive.
	offs <- fetwfe:::.resolve_cohort_offsets_and_first_inds(
		f,
		G = f$G,
		T = f$T
	)
	expect_identical(as.integer(offs$cohort_offsets_int), c(2L, 4L, 7L))
	expect_false(identical(
		as.integer(offs$first_inds),
		as.integer(getFirstInds(G = f$G, T = f$T))
	))
	db <- debiasedATT(f)
	expect_true(is.finite(db$att) && is.finite(db$se) && db$se > 0)
	# the ATT-direction identity now holds (acceptance criterion, 1e-6 tol).
	expect_lt(.att_identity_gap(f), 1e-6)
})

test_that("debiasedATT() runs on a scattered-cohort HIGH-dim fit (#318)", {
	f <- .fit_318(.mk_panel_318(c(2L, 4L, 7L), d = 10L), gls = FALSE)
	expect_gte(ncol(f$internal$X_final), nrow(f$internal$X_final)) # high-dim
	db <- debiasedATT(f)
	expect_true(is.finite(db$att) && is.finite(db$se) && db$se > 0)
	expect_lt(.att_identity_gap(f), 1e-6)
})

test_that("consecutive-cohort fits are byte-identical (the resolver matches getFirstInds) (#318)", {
	f <- .fit_318(
		.mk_panel_318(c(2L, 3L, 4L), d = 3L),
		gls = TRUE,
		sig_eps_sq = 0.16,
		sig_eps_c_sq = 0.1
	)
	# consecutive adoption: the resolver returns exactly getFirstInds() and the
	# block sizes reduce to (T-1):(T-G), so the construction is unchanged.
	fi_res <- fetwfe:::.resolve_cohort_offsets_and_first_inds(
		f,
		G = f$G,
		T = f$T
	)$first_inds
	expect_identical(
		as.integer(fi_res),
		as.integer(getFirstInds(G = f$G, T = f$T))
	)
	nt <- length(f$treat_inds)
	expect_identical(
		as.integer(diff(c(fi_res, nt + 1L))),
		as.integer((f$T - 1L):(f$T - f$G))
	)
	db <- debiasedATT(f)
	expect_true(is.finite(db$att) && db$se > 0)
})
