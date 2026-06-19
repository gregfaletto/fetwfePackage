# Regression test for #323: the high-dim `lambda_c = "cv"` band must CV its penalty
# constant on the SAME overall-ATT direction debiasedATT() uses, INCLUDING under
# scattered (non-consecutive) cohort adoption.
#
# `.simultaneous_cis_impl()` built the CV-selection direction `a_att` with hard-coded
# consecutive-adoption block sizes `(T-1):(T-G)`. Under scattered adoption (offsets
# 2,4,7 here) the real per-cohort sizes differ -- (7,5,2) vs the consecutive (7,6,5)
# -- and the consecutive sizes both mis-weight AND OVERRUN `treat_inds` (they sum to
# 18 > num_treats = 14), so the band CV'd a wrong/fallback constant and diverged from
# debiasedATT(). #318 had made debiasedATT() support scattered cohorts, invalidating
# the stale comment that claimed "debiasedATT() errors on scattered, so no divergence
# can ship". The fix uses `diff(c(first_inds, num_treats + 1L))` (the resolver's
# sizes, matching debiased_att.R) -- byte-identical on consecutive adoption.
#
# Mutation-checkable: reverting the construction at R/simultaneous_cis.R to
# `(T-1):(T-G)` makes the band CV a different constant on this scattered fit, so the
# overall-ATT band center no longer equals debiasedATT(lambda_c = "cv")$att (the #295
# "one constant for point + band" identity) -> the final expect_equal FAILs.

.mk_panel_323 <- function(
	adopt = c(2L, 4L, 7L),
	d = 10L,
	N = 24L,
	T = 8L,
	seed = 5
) {
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

test_that("high-dim lambda_c='cv' band CVs the same constant as debiasedATT() on scattered cohorts (#323)", {
	fit <- fetwfe(
		pdata = .mk_panel_323(),
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:10),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
	# fixture invariants: genuinely high-dim full design + scattered adoption.
	expect_gte(ncol(fit$internal$X_final), nrow(fit$internal$X_final))
	G <- fit$G
	T_ <- fit$T
	num_treats <- length(fit$treat_inds)
	first_inds <- fetwfe:::.resolve_cohort_offsets_and_first_inds(
		fit,
		G = G,
		T = T_
	)$first_inds
	sizes <- diff(c(first_inds, num_treats + 1L))
	# the bug's root: the resolved sizes differ from the consecutive assumption AND
	# the consecutive sizes overrun treat_inds under scattered adoption.
	expect_identical(as.integer(sizes), c(7L, 5L, 2L))
	expect_false(identical(as.integer(sizes), as.integer((T_ - 1L):(T_ - G))))
	expect_gt(sum((T_ - 1L):(T_ - G)), num_treats)

	# overall-ATT custom contrast, built with the CORRECT resolved sizes.
	cohort_of_treat <- rep(seq_len(G), times = sizes)
	contrast <- numeric(num_treats)
	for (g in seq_len(G)) {
		idx <- which(cohort_of_treat == g)
		contrast[idx] <- fit$cohort_probs[g] / length(idx)
	}
	expect_equal(sum(contrast), sum(fit$cohort_probs[seq_len(G)]))

	sc <- suppressWarnings(simultaneousCIs(
		fit,
		family = "custom",
		contrasts = matrix(contrast, nrow = 1L),
		method = "bootstrap",
		B = 200,
		seed = 1,
		lambda_c = "cv"
	))
	expect_identical(sc$regime, "high-dimensional")
	db <- debiasedATT(fit, lambda_c = "cv")
	# #295 identity: one CV'd lambda_c serves point + band, so the overall-ATT band
	# center equals debiasedATT()$att. The buggy consecutive-size a_att CVs a
	# different constant on scattered data and breaks this.
	expect_equal(sc$ci$estimate, db$att, tolerance = 1e-7)
})
