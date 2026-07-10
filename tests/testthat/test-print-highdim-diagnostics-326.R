# Tests for #326: print.simultaneous_cis / print.debiased_att surface the
# experimental high-dimensional (p >= NT) desparsified-band diagnostics, and
# debiasedATT() returns a classed "debiased_att" object.
#
# Both are additive UX changes: no behavior change to returned values. The class
# is added WITHOUT a new list element (the regime is inferred from `lambda_node`),
# so debiasedATT()'s `names()` / `$`-accessor contract is unchanged (pinned by the
# back-compat assertions below + the exact-names tests in test-debiased-att*.R).
# Output is checked by substring (not snapshot) so the seed-dependent numbers in
# the diagnostics aren't a source of fragility.

# Shared fixtures (built once): a fixed-p fit and a high-dim scattered fit.
.fp_fit_326 <- local({
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 1)
	sim <- simulateData(
		cf,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	fetwfe(
		pdata = sim$pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		q = 0.5,
		verbose = FALSE
	)
})

.hd_fit_326 <- local({
	set.seed(5)
	d <- 10L
	N <- 24L
	T <- 8L
	adopt <- c(2L, 4L, 7L)
	covs <- matrix(stats::rnorm(N * d), N, d)
	cou <- c(
		rep(0L, N - 18L),
		rep(adopt[1], 6),
		rep(adopt[2], 6),
		rep(adopt[3], 6)
	)
	eff <- stats::setNames(c(0.5, 2, 3.5), as.character(adopt))
	panel <- do.call(
		rbind,
		lapply(seq_len(N), function(i) {
			g <- cou[i]
			df <- data.frame(
				unit = sprintf("u%02d", i),
				year = 1:T,
				treat = as.integer(g > 0 & (1:T) >= g)
			)
			for (j in 1:d) {
				df[[paste0("x", j)]] <- covs[i, j]
			}
			df$y <- 0.3 *
				(1:T) /
				T +
				(if (g > 0) eff[[as.character(g)]] else 0) *
					df$treat +
				0.2 * covs[i, 1] +
				stats::rnorm(T, 0, 0.4)
			df
		})
	)
	fetwfe(
		pdata = panel,
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:10),
		response = "y",
		q = 0.5,
		verbose = FALSE,
		gls = FALSE
	)
})

test_that("debiasedATT() returns a classed, back-compat 'debiased_att' object (#326)", {
	db <- debiasedATT(.fp_fit_326)
	expect_s3_class(db, "debiased_att")
	# back-compat: the fixed-p list contents are UNCHANGED (no new field; same
	# names; all length-1 numeric). This is the contract test-debiased-att.R pins.
	expect_identical(
		names(db),
		c("att", "se", "ci_low", "ci_high", "var_reg", "var_weight")
	)
	expect_true(all(vapply(
		db,
		function(z) is.numeric(z) && length(z) == 1L,
		NA
	)))
	# tidy surface matches the tidy.simultaneous_cis column convention.
	td <- tidy(db)
	expect_s3_class(td, "data.frame")
	expect_identical(
		names(td),
		c("term", "estimate", "std.error", "conf.low", "conf.high")
	)
	expect_equal(td$estimate, db$att)
	expect_equal(td$conf.low, db$ci_low)
})

test_that("print.debiased_att shows the high-dim diagnostics only in the high-dim regime (#326)", {
	out_fp <- capture.output(print(debiasedATT(.fp_fit_326)))
	expect_true(any(grepl("Debiased overall ATT \\(fixed-p regime\\)", out_fp)))
	expect_true(any(grepl("Estimate:", out_fp)))
	# fixed-p: NO experimental diagnostics block.
	expect_false(any(grepl("EXPERIMENTAL|nodewise", out_fp)))

	db_hd <- suppressWarnings(debiasedATT(.hd_fit_326, lambda_c = "cv"))
	expect_s3_class(db_hd, "debiased_att")
	out_hd <- capture.output(print(db_hd))
	expect_true(any(grepl("high-dimensional regime", out_hd)))
	# The scalar overall-ATT interval was de-experimented in #387 (Theorem
	# `debiased.highdim.thm` validated in simulation): the interval banner no
	# longer says EXPERIMENTAL -- it cites the theorem + validated coverage. (The
	# family-wise *band* banner keeps EXPERIMENTAL; see the simultaneous_cis test.)
	expect_false(any(grepl("EXPERIMENTAL", out_hd)))
	expect_true(any(grepl("validated near-nominally", out_hd)))
	expect_true(any(grepl("nodewise penalty: lambda_c", out_hd)))
	expect_true(any(grepl("converged.*KKT-feasible", out_hd)))

	# #343: tidy.debiased_att in the HIGH-DIM regime (the fixed-p tidy schema is
	# pinned above; confirm the high-dim object stays broom-compatible too).
	td_hd <- tidy(db_hd)
	expect_s3_class(td_hd, "data.frame")
	expect_equal(nrow(td_hd), 1L)
	expect_named(
		td_hd,
		c("term", "estimate", "std.error", "conf.low", "conf.high")
	)
})

test_that("print.simultaneous_cis shows the high-dim diagnostics block only in the high-dim regime (#326)", {
	out_fp <- capture.output(print(simultaneousCIs(
		.fp_fit_326,
		family = "cohort"
	)))
	# fixed-p analytic band: unchanged print, no diagnostics block.
	expect_false(any(grepl("EXPERIMENTAL|nodewise penalty", out_fp)))

	sc_hd <- suppressWarnings(simultaneousCIs(
		.hd_fit_326,
		family = "cohort",
		method = "bootstrap",
		B = 200,
		seed = 1,
		lambda_c = "cv"
	))
	expect_identical(sc_hd$regime, "high-dimensional")
	out_hd <- capture.output(print(sc_hd))
	expect_true(any(grepl(
		"High-dimensional \\(p >= NT\\) desparsified band",
		out_hd
	)))
	expect_true(any(grepl(
		"nodewise directions:.*converged.*KKT-feasible",
		out_hd
	)))
	# the CI table is still printed (the diagnostics are additive).
	expect_true(any(grepl("simultaneous_ci_low", out_hd)))
})

test_that("the high-dim event_study propensity-channel feasibility count uses the gate's slack (#326)", {
	# The propensity (#309) channel's KKT-feasible count must use the same
	# feasibility tolerance (via `.riesz_feasible()`) as the per-effect line and
	# the band's warning gate. The
	# nodewise solver binds the constraint to ~1e-10, so a strict `<=` reports
	# 0/K feasible exactly when no feasibility warning fired -- the most alarming
	# possible output, shown when everything is fine.
	sc_es <- suppressWarnings(simultaneousCIs(
		.hd_fit_326,
		family = "event_study",
		method = "bootstrap",
		B = 200,
		seed = 1,
		lambda_c = "cv"
	))
	expect_identical(sc_es$regime, "high-dimensional")
	expect_false(is.null(sc_es$propensity_feasibility)) # the channel exists here
	out_es <- capture.output(print(sc_es))
	prop_line <- grep("propensity channel", out_es, value = TRUE)
	expect_length(prop_line, 1L)
	# parse "... X/Y KKT-feasible": with the slack X is the gate's feasible count
	# (most directions); a strict `<=` (the bug) collapses it to 0.
	feas <- as.integer(sub(
		".*converged, ([0-9]+)/[0-9]+ KKT-feasible.*",
		"\\1",
		prop_line
	))
	expect_gt(feas, 0L)
})
