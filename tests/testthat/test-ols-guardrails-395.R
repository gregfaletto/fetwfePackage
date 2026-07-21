# Tests for #395: OLS-pair (etwfe/twfeCovs) guard rails.
# (1) A rank-deficient design (e.g. perfectly collinear covariates) now errors
#     informatively for BOTH OLS estimators, instead of the bare
#     `all(!is.na(beta_hat_slopes)) is not TRUE` assertion.
# (2) The d+1-per-cohort gate stays for etwfe() (its full design needs it) but is
#     skipped for twfeCovs() (collapsed design). The twfeCovs-succeeds case is in
#     test-twfeCovs.R Test 30; here we check the etwfe() gate still fires.

.mk_panel_395 <- function(cohort_of_unit, T = 4L, collinear = FALSE, seed = 1) {
	set.seed(seed)
	rows <- list()
	for (i in seq_along(cohort_of_unit)) {
		g <- cohort_of_unit[i]
		x1 <- stats::rnorm(1)
		for (t in 1:T) {
			treat <- as.integer(g > 0 & t >= g)
			row <- data.frame(
				unit = sprintf("u%02d", i),
				year = t,
				treat = treat,
				x1 = x1
			)
			if (collinear) {
				row$x2 <- 2 * x1
			}
			row$y <- 0.3 * t + 1.5 * treat + 0.5 * x1 + stats::rnorm(1, 0, 0.3)
			rows[[length(rows) + 1L]] <- row
		}
	}
	do.call(rbind, rows)
}

test_that("collinear covariates yield an informative rank-deficiency error, not a bare assertion (#395)", {
	# Cohorts of 4 units (>= d + 1 = 3), so the d+1 gate does NOT pre-empt: the
	# design reaches lm(), which aliases the collinear column -> the informative
	# rank-deficiency error (defect 1), for both OLS estimators.
	df <- .mk_panel_395(c(rep(0L, 4), rep(2L, 4), rep(3L, 4)), collinear = TRUE)
	args <- list(
		pdata = df,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		response = "y",
		covs = c("x1", "x2"),
		verbose = FALSE
	)
	# New message names the cause (+ the add_ridge escape hatch); it is NOT the
	# bare `all(!is.na(beta_hat_slopes)) is not TRUE` assertion.
	expect_error(do.call(etwfe, args), "rank-deficient")
	expect_error(do.call(twfeCovs, args), "rank-deficient")
	# The penalized estimator fits the identical collinear input fine (contract).
	expect_s3_class(do.call(fetwfe, c(args, list(q = 0.5))), "fetwfe")
})

test_that("the d+1-per-cohort gate still fires for etwfe() (its full design needs it) (#395)", {
	# A 1-unit cohort (d = 1): etwfe()'s full ETWFE design genuinely needs d + 1
	# units per cohort, so the gate legitimately stops -- unchanged by #395 (only
	# twfeCovs() skips it).
	df <- .mk_panel_395(c(rep(0L, 4), rep(2L, 6), 3L), T = 4L, seed = 2)
	expect_error(
		etwfe(
			pdata = df,
			time_var = "year",
			unit_var = "unit",
			treatment = "treat",
			response = "y",
			covs = "x1",
			verbose = FALSE
		),
		"fewer than d \\+ 1 units"
	)
})
