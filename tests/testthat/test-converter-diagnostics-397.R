# Tests for #397: converter/validator input diagnostics.
# (1) attgtToFetwfeDf() no longer emits the always-false "treatment reversal"
#     warning (it measured input row order, not a real reversal, which is
#     impossible by construction); output is row-order-invariant.
# (2) Inf-coded never-treated (fect/staggered/some did) converts as never-treated
#     (mapped to 0), instead of the misleading "Missing values" error.
# (3)/(4) The estimator input validator names NA response / NA time / mixed-NA
#     indep_counts violations instead of a bare internal assertion or a crash.

.attgt_397 <- local({
	set.seed(397)
	data.frame(
		unit = rep(1:6, each = 4),
		period = rep(1:4, times = 6),
		G = rep(c(0, 0, 2, 2, 3, 3), each = 4),
		Y = stats::rnorm(24),
		x1 = rep(stats::rnorm(6), each = 4)
	)
})
.conv_397 <- function(d) {
	attgtToFetwfeDf(
		d,
		yname = "Y",
		tname = "period",
		idname = "unit",
		gname = "G",
		covars = "x1"
	)
}

test_that("attgtToFetwfeDf: no false 'reversal' warning; output is row-order-invariant; Inf accepted (#397)", {
	asc <- .conv_397(.attgt_397[order(.attgt_397$period), ])

	# Time-descending input previously tripped the (always-false) reversal
	# warning; it must not now, and the converted output is identical regardless
	# of input row order (the converter sorts internally).
	w <- testthat::capture_warnings(
		desc <- .conv_397(.attgt_397[order(-.attgt_397$period), ])
	)
	expect_false(any(grepl("switch from treated back", w)))
	expect_identical(asc, desc)

	# Inf-coded never-treated converts identically to the 0-coded panel (was:
	# "Missing values in `G`." because as.integer(Inf) is NA).
	attgt_inf <- .attgt_397
	attgt_inf$G[attgt_inf$G == 0] <- Inf
	expect_identical(.conv_397(.attgt_397), .conv_397(attgt_inf))
})

test_that("the estimator input validator names NA response / time / indep_counts violations (#397)", {
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 1)
	sim <- simulateData(
		cf,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	base <- list(
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		verbose = FALSE
	)

	# NA response: named violation, NOT the bare `all(!is.na(data)) is not TRUE`.
	pd_resp <- sim$pdata
	pd_resp[[sim$response]][1] <- NA
	expect_error(
		do.call(etwfe, c(list(pdata = pd_resp), base)),
		"response column.*must not contain missing"
	)

	# NA time: named violation, not a confusing idCohorts() balance error.
	pd_time <- sim$pdata
	pd_time[[sim$time_var]][1] <- NA_integer_
	expect_error(
		do.call(etwfe, c(list(pdata = pd_time), base)),
		"time_var column.*must not contain missing"
	)

	# Mixed-NA indep_counts: named violation, not "missing value where TRUE/FALSE
	# needed" (the validator formerly crashed on `any(c(18, NA, 6) <= 0)` -> NA).
	expect_error(
		do.call(
			etwfe,
			c(list(pdata = sim$pdata, indep_counts = c(18L, NA, 6L)), base)
		),
		"indep_counts must not contain missing"
	)

	# A normal fit (no NA anywhere) is unaffected.
	expect_s3_class(do.call(etwfe, c(list(pdata = sim$pdata), base)), "etwfe")
})
