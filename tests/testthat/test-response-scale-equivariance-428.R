# Response-scale equivariance (#428). Multiplying the response by k > 0, and any
# supplied variance by k^2, must multiply every estimate, standard error,
# interval bound and coefficient by k, every lambda by k^(2 - q) and every
# variance by k^2, and leave model sizes, critical values and p-values unchanged.
# Each block compares the fits at k against the fit at k = 1, after dividing out
# that power of k, and first asserts that the k = 1 reference is not degenerate:
# a null model is equivariant trivially.
#
# Two exceptions, each cited where it applies: simultaneous-band quantities are
# compared only at scales where the band's absolute variance floor cannot bind
# (#489), and `gls = FALSE` fits keep the unstandardized lambda grid (#490).
#
# The fits call the estimators directly rather than the `*WithSimulatedData()`
# wrappers, which pass the simulator's true variances through and would skip
# the variance estimation that is part of the pipeline under test.

.RSE428_TOL <- 1e-6
# Scales at which pointwise quantities are compared with k = 1.
.RSE428_SCALES <- c(1e-9, 1e-6, 1e-3, 1e3, 1e6)
# Scales at which simultaneous-band quantities are compared with k = 1 (#489).
.RSE428_BAND_SCALES <- c(1e-2, 1e3, 1e6)

.rse428_sim <- simulateData(
	genCoefs(G = 3, T = 4, d = 2, density = 0.5, eff_size = 2, seed = 123),
	N = 120,
	sig_eps_sq = 0.5,
	sig_eps_c_sq = 0.5,
	seed = 456
)

# `estimator` fitted on `sim` with its response multiplied by `k`; `...` is
# passed to the estimator. Warnings are not suppressed: an unexpected one at an
# extreme scale is one of the things this file exists to notice.
.rse428_scaled_fit <- function(estimator, k, ..., sim = .rse428_sim) {
	pdata <- sim$pdata
	pdata[[sim$response]] <- pdata[[sim$response]] * k
	estimator(
		pdata = pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		...
	)
}

# `make(k)` at k = 1 and at each of `scales`, built now; returns the accessor
# `at(k)`, which errors on a scale it was not built at.
.rse428_fits <- function(make, scales) {
	ks <- c(1, scales)
	fits <- stats::setNames(lapply(ks, make), as.character(ks))
	function(k) {
		key <- as.character(k)
		if (!key %in% names(fits)) {
			stop("no fit at k = ", key)
		}
		fits[[key]]
	}
}

# The default-route `fetwfe()` fits, with default and with pointwise intervals,
# that several blocks share; each is built on first use.
.rse428_cache <- new.env(parent = emptyenv())
.rse428_shared <- function(config, k) {
	id <- paste0(config, "@", as.character(k))
	if (!exists(id, envir = .rse428_cache, inherits = FALSE)) {
		fit <- switch(
			config,
			default = .rse428_scaled_fit(fetwfe, k),
			pointwise = .rse428_scaled_fit(fetwfe, k, ci_type = "pointwise"),
			stop("unknown configuration: ", config)
		)
		assign(id, fit, envir = .rse428_cache)
	}
	get(id, envir = .rse428_cache, inherits = FALSE)
}

.rse428_q <- function(power, get) {
	list(power = power, get = get)
}

# Expects each quantity of `at(k)`, divided by k^power, to equal the same
# quantity of `at(1)`, at every k in `scales`.
.rse428_expect_equivariant <- function(at, quantities, scales) {
	ref <- at(1)
	for (k in scales) {
		obj <- at(k)
		for (nm in names(quantities)) {
			spec <- quantities[[nm]]
			expect_equal(
				spec$get(obj) / k^spec$power,
				spec$get(ref),
				tolerance = .RSE428_TOL,
				info = sprintf("%s at k = %g", nm, k)
			)
		}
	}
}

# The pointwise quantities of a fit. `bounds` reads `catt_df`'s interval, so it
# belongs to `ci_type = "pointwise"` fits only. A fit returns a supplied
# variance as given, so `variances = FALSE` is for fits that supply them.
.rse428_fit_quantities <- function(
	q = 0.5,
	penalized = TRUE,
	ses = TRUE,
	bounds = TRUE,
	variances = TRUE
) {
	out <- list(
		att_hat = .rse428_q(1, function(f) f$att_hat),
		beta_hat = .rse428_q(1, function(f) f$beta_hat),
		catt_hats = .rse428_q(1, function(f) f$catt_hats)
	)
	if (variances) {
		out$sig_eps_sq <- .rse428_q(2, function(f) f$sig_eps_sq)
		out$sig_eps_c_sq <- .rse428_q(2, function(f) f$sig_eps_c_sq)
	}
	if (ses) {
		out$att_se <- .rse428_q(1, function(f) f$att_se)
		out$catt_ses <- .rse428_q(1, function(f) f$catt_ses)
	}
	if (bounds) {
		out$ci_low <- .rse428_q(1, function(f) f$catt_df$ci_low)
		out$ci_high <- .rse428_q(1, function(f) f$catt_df$ci_high)
	}
	if (penalized) {
		out$lambda_star <- .rse428_q(2 - q, function(f) f$lambda_star)
		out$lambda.max <- .rse428_q(2 - q, function(f) f$lambda.max)
		out$lambda.min <- .rse428_q(2 - q, function(f) f$lambda.min)
		out$lambda_star_model_size <- .rse428_q(
			0,
			function(f) f$lambda_star_model_size
		)
	}
	out
}

# The simultaneous-band quantities of a default fit's `catt_df`.
.rse428_band_quantities <- list(
	ci_low = .rse428_q(1, function(f) f$catt_df$ci_low),
	ci_high = .rse428_q(1, function(f) f$catt_df$ci_high),
	p_value = .rse428_q(0, function(f) f$catt_df$p_value)
)

# The reference is not degenerate: it selected something and has a nonzero
# CATT. An unpenalized fit selects nothing out, so it gets the CATT check only.
.rse428_expect_live <- function(ref, penalized = TRUE) {
	if (penalized) {
		expect_gt(ref$lambda_star_model_size, 0L)
	}
	expect_true(any(ref$catt_hats != 0))
}

# The precondition of a band comparison (#489): at every scale it uses, the
# family's smallest variance exceeds the floor `var_tol` takes when every
# variance is below 1, so the floor cannot class an effect degenerate.
.rse428_expect_band_ready <- function(at, scales, ses) {
	for (k in c(1, scales)) {
		expect_gt(
			min(ses(at(k))^2),
			sqrt(.Machine$double.eps),
			label = sprintf("smallest squared SE at k = %g", k)
		)
	}
}

# The q = 1 nuisance that `debiasedATT()` and the bootstrap band refit when
# p >= NT is not null, so their comparisons cannot pass trivially.
.rse428_expect_nuisance_live <- function(ref) {
	X <- ref$internal$X_final
	y <- as.numeric(ref$internal$y_final)[seq_len(ref$N * ref$T)]
	expect_true(any(fetwfe:::.fit_q1_nuisance(X, y, ref$N, ref$T)[-1] != 0))
}

# A penalized estimator on one selection route. `pw` and `df` are accessors from
# `.rse428_fits()`, of `ci_type = "pointwise"` and of default fits.
.rse428_expect_route <- function(pw, df, band = .rse428_band_quantities) {
	.rse428_expect_live(pw(1))
	.rse428_expect_equivariant(pw, .rse428_fit_quantities(), .RSE428_SCALES)
	# Band quantities only where the variance floor cannot bind (#489).
	.rse428_expect_live(df(1))
	.rse428_expect_band_ready(df, .RSE428_BAND_SCALES, function(f) f$catt_ses)
	.rse428_expect_equivariant(df, band, .RSE428_BAND_SCALES)
}

test_that("fetwfe() on the default route (CV, q = 0.5) is equivariant (#428)", {
	.rse428_expect_route(
		pw = .rse428_fits(
			function(k) .rse428_shared("pointwise", k),
			.RSE428_SCALES
		),
		df = .rse428_fits(
			function(k) .rse428_shared("default", k),
			.RSE428_BAND_SCALES
		)
	)
})

test_that("fetwfe() with lambda_selection = \"bic\" is equivariant (#428)", {
	.rse428_expect_route(
		pw = .rse428_fits(
			function(k) {
				.rse428_scaled_fit(
					fetwfe,
					k,
					lambda_selection = "bic",
					ci_type = "pointwise"
				)
			},
			.RSE428_SCALES
		),
		df = .rse428_fits(
			function(k) .rse428_scaled_fit(fetwfe, k, lambda_selection = "bic"),
			.RSE428_BAND_SCALES
		),
		band = .rse428_band_quantities[c("ci_low", "ci_high")]
	)
})

test_that("betwfe() is equivariant on both selection routes (#428)", {
	for (route in c("cv", "bic")) {
		.rse428_expect_route(
			pw = .rse428_fits(
				function(k) {
					.rse428_scaled_fit(
						betwfe,
						k,
						lambda_selection = route,
						ci_type = "pointwise"
					)
				},
				.RSE428_SCALES
			),
			df = .rse428_fits(
				function(k) {
					.rse428_scaled_fit(betwfe, k, lambda_selection = route)
				},
				.RSE428_BAND_SCALES
			)
		)
	}
})

# add_ridge = TRUE on all four estimators, one block each, so a failure names
# its estimator.
.rse428_ridge_rows <- list(
	fetwfe = list(estimator = fetwfe, penalized = TRUE),
	betwfe = list(estimator = betwfe, penalized = TRUE),
	etwfe = list(estimator = etwfe, penalized = FALSE),
	twfeCovs = list(estimator = twfeCovs, penalized = FALSE)
)
for (.rse428_row in names(.rse428_ridge_rows)) {
	test_that(
		sprintf("%s(add_ridge = TRUE) is equivariant (#428)", .rse428_row),
		{
			row <- .rse428_ridge_rows[[.rse428_row]]
			at <- .rse428_fits(
				function(k) {
					.rse428_scaled_fit(
						row$estimator,
						k,
						add_ridge = TRUE,
						ci_type = "pointwise"
					)
				},
				.RSE428_SCALES
			)
			.rse428_expect_live(at(1), penalized = row$penalized)
			.rse428_expect_equivariant(
				at,
				.rse428_fit_quantities(penalized = row$penalized),
				.RSE428_SCALES
			)
		}
	)
}

test_that("fetwfe() with supplied variances rescaled with the response is equivariant (#428)", {
	at <- .rse428_fits(
		function(k) {
			.rse428_scaled_fit(
				fetwfe,
				k,
				sig_eps_sq = 0.5 * k^2,
				sig_eps_c_sq = 0.5 * k^2,
				ci_type = "pointwise"
			)
		},
		.RSE428_SCALES
	)
	.rse428_expect_live(at(1))
	.rse428_expect_equivariant(
		at,
		.rse428_fit_quantities(variances = FALSE),
		.RSE428_SCALES
	)
})

test_that("fetwfe(gls = FALSE) keeps the unstandardized grid: non-regression pin (#490)", {
	# A gls = FALSE fit has no noise variance to standardize by, so it keeps the
	# unstandardized grid (#490). The values origin/main returns at 7179af7.
	fit <- .rse428_scaled_fit(fetwfe, 1, gls = FALSE)
	expect_equal(fit$att_hat, -1.40732975154555029, tolerance = .RSE428_TOL)
	expect_equal(fit$lambda_star, 0.00918281458106647, tolerance = .RSE428_TOL)
	expect_identical(fit$lambda_star_model_size, 26L)
	expect_equal(
		fit$beta_hat,
		c(
			2.19615299720982,
			2.19615299720982,
			0,
			3.8267095458097,
			1.89762581065613,
			0,
			1.99669787032021,
			1.92575218665297,
			6.07871321915834,
			2.71025519933139,
			3.96541945686076,
			1.72683597050094,
			2.26139313749105,
			1.72683597050094,
			0,
			4.20975053605285,
			0,
			2.18775908186378,
			0,
			1.89061238365982,
			-1.70156921162222,
			0.214200113177695,
			2.03652964817023,
			-3.8944685390395,
			-2.24506491738206,
			-2.05358550672409,
			2.35069506263859,
			0.537917712025784,
			2.35069506263859,
			0.537917712025784,
			0.221154547650776,
			-0.49262956958145,
			4.28042141324742,
			0.537917712025784,
			3.89810003947083,
			0.537917712025784,
			5.95655388270437,
			0.537917712025784
		),
		tolerance = .RSE428_TOL
	)
})

# The other two penalty exponents. Standard errors and bounds are `NA` at
# q >= 1, so neither is compared.
for (.rse428_q_val in c(1, 2)) {
	test_that(
		sprintf("fetwfe(q = %g) is equivariant (#428)", .rse428_q_val),
		{
			q <- .rse428_q_val
			at <- .rse428_fits(
				function(k) .rse428_scaled_fit(fetwfe, k, q = q),
				.RSE428_SCALES
			)
			.rse428_expect_live(at(1))
			.rse428_expect_equivariant(
				at,
				.rse428_fit_quantities(q = q, ses = FALSE, bounds = FALSE),
				.RSE428_SCALES
			)
		}
	)
}

test_that("eventStudy(), cohortStudy() and cohortTimeATTs() are equivariant (#428)", {
	df <- function(k) .rse428_shared("default", k)
	pw <- function(k) .rse428_shared("pointwise", k)
	.rse428_expect_live(df(1))
	.rse428_expect_live(pw(1))
	es <- function(k) eventStudy(df(k))
	.rse428_expect_equivariant(
		es,
		list(
			estimate = .rse428_q(1, function(x) x$estimate),
			se = .rse428_q(1, function(x) x$se)
		),
		.RSE428_SCALES
	)
	.rse428_expect_equivariant(
		function(k) eventStudy(df(k), ci_type = "pointwise"),
		list(
			ci_low = .rse428_q(1, function(x) x$ci_low),
			ci_high = .rse428_q(1, function(x) x$ci_high)
		),
		.RSE428_SCALES
	)
	# The event-study band only where the variance floor cannot bind (#489).
	.rse428_expect_band_ready(es, .RSE428_BAND_SCALES, function(x) x$se)
	.rse428_expect_equivariant(
		es,
		list(
			ci_low = .rse428_q(1, function(x) x$ci_low),
			ci_high = .rse428_q(1, function(x) x$ci_high)
		),
		.RSE428_BAND_SCALES
	)
	.rse428_expect_equivariant(
		function(k) cohortStudy(pw(k)),
		list(
			estimate = .rse428_q(1, function(x) x$estimate),
			se = .rse428_q(1, function(x) x$se),
			ci_low = .rse428_q(1, function(x) x$ci_low),
			ci_high = .rse428_q(1, function(x) x$ci_high)
		),
		.RSE428_SCALES
	)
	.rse428_expect_equivariant(
		function(k) cohortTimeATTs(df(k)),
		list(
			estimate = .rse428_q(1, function(x) x$estimate),
			se = .rse428_q(1, function(x) x$se)
		),
		.RSE428_SCALES
	)
})

test_that("simultaneousCIs() is equivariant under both methods (#428)", {
	df <- function(k) .rse428_shared("default", k)
	.rse428_expect_live(df(1))
	# The band only where the variance floor cannot bind (#489).
	.rse428_expect_band_ready(
		df,
		.RSE428_BAND_SCALES,
		function(f) eventStudy(f)$se
	)
	bands <- list(
		analytic = function(k) simultaneousCIs(df(k), method = "analytic"),
		bootstrap = function(k) {
			simultaneousCIs(df(k), method = "bootstrap", seed = 1)
		}
	)
	for (band in bands) {
		.rse428_expect_equivariant(
			band,
			list(
				estimate = .rse428_q(1, function(x) x$ci$estimate),
				pointwise_ci_low = .rse428_q(1, function(x) {
					x$ci$pointwise_ci_low
				}),
				pointwise_ci_high = .rse428_q(1, function(x) {
					x$ci$pointwise_ci_high
				})
			),
			.RSE428_SCALES
		)
		.rse428_expect_equivariant(
			band,
			list(
				simultaneous_ci_low = .rse428_q(
					1,
					function(x) x$ci$simultaneous_ci_low
				),
				simultaneous_ci_high = .rse428_q(
					1,
					function(x) x$ci$simultaneous_ci_high
				),
				critical_value = .rse428_q(0, function(x) x$critical_value),
				adjusted_p_values = .rse428_q(0, function(x) {
					x$adjusted_p_values
				})
			),
			.RSE428_BAND_SCALES
		)
	}
})

test_that("debiasedATT() is equivariant under both methods (#428)", {
	df <- function(k) .rse428_shared("default", k)
	.rse428_expect_live(df(1))
	.rse428_expect_equivariant(
		function(k) debiasedATT(df(k)),
		list(
			att = .rse428_q(1, function(x) x$att),
			se = .rse428_q(1, function(x) x$se),
			var_reg = .rse428_q(2, function(x) x$var_reg),
			var_weight = .rse428_q(2, function(x) x$var_weight)
		),
		.RSE428_SCALES
	)
	.rse428_expect_equivariant(
		function(k) debiasedATT(df(k), method = "bootstrap", seed = 1),
		list(
			att = .rse428_q(1, function(x) x$att),
			se = .rse428_q(1, function(x) x$se),
			ci_low = .rse428_q(1, function(x) x$ci_low),
			ci_high = .rse428_q(1, function(x) x$ci_high),
			crit_value = .rse428_q(0, function(x) x$crit_value)
		),
		.RSE428_SCALES
	)
})

test_that("a supplied lambda.max keeps its meaning on the BIC route: non-regression pin (#428)", {
	# The values origin/main returns at 7179af7.
	fit <- .rse428_scaled_fit(
		fetwfe,
		1,
		lambda_selection = "bic",
		lambda.max = 2
	)
	expect_equal(fit$att_hat, -0.9828402446, tolerance = .RSE428_TOL)
	expect_equal(fit$att_se, 0.2132159268, tolerance = .RSE428_TOL)
	expect_identical(fit$lambda_star_model_size, 13L)
	expect_equal(fit$lambda.max, 2, tolerance = 1e-12)
	# gBridge's default lambda.min is the fraction 0.001 of lambda.max when n > p.
	expect_equal(fit$lambda.min, 0.002, tolerance = 1e-12)
})

test_that("a supplied lambda.min keeps its meaning on the BIC route: non-regression pin (#428)", {
	fit <- .rse428_scaled_fit(
		fetwfe,
		1,
		lambda_selection = "bic",
		lambda.max = 2,
		lambda.min = 0.01
	)
	expect_equal(fit$lambda.min, 0.02, tolerance = 1e-12)
})

test_that(".bridge_response_scale() returns the noise SD, else the response SD, else 1 (#428)", {
	scale_of <- fetwfe:::.bridge_response_scale
	y <- c(1, 3)
	# sqrt(sig_eps_sq) for a finite positive variance, whatever y holds.
	expect_identical(scale_of(4), 2)
	expect_identical(scale_of(0.25, y = y), 0.5)
	# sd(y) = sqrt(2) when sig_eps_sq is not a finite positive number.
	for (v in list(NA, NA_real_, 0, -1, Inf)) {
		expect_equal(
			scale_of(v, y = y),
			sqrt(2),
			info = sprintf("sig_eps_sq = %s", format(v))
		)
	}
	# 1 when there is no y, or y is constant.
	expect_identical(scale_of(NA_real_), 1)
	expect_identical(scale_of(), 1)
	expect_identical(scale_of(NA_real_, y = rep(3, 5)), 1)
})

test_that("fetwfe() rejects add_ridge = TRUE with gls = FALSE without a variance rationale (#428)", {
	msg <- msg_of(tryCatch(
		.rse428_scaled_fit(fetwfe, 1, add_ridge = TRUE, gls = FALSE),
		error = identity
	))
	expect_identical(
		msg,
		"fetwfe(): `add_ridge = TRUE` is not supported with `gls = FALSE`. Re-fit with `add_ridge = FALSE`."
	)
})

test_that("fetwfe() and betwfe() standardize by the sig_eps_sq their GLS step used (#428)", {
	# An equivariance assertion cannot tell which variance a core hands the
	# dispatcher, so record it.
	real <- fetwfe:::.bridge_response_scale
	for (est in c("fetwfe", "betwfe")) {
		rec <- new.env(parent = emptyenv())
		rec$sig_eps_sq <- list()
		fit <- testthat::with_mocked_bindings(
			.rse428_scaled_fit(get(est), 1),
			.bridge_response_scale = function(sig_eps_sq = NA_real_, y = NULL) {
				rec$sig_eps_sq[[length(rec$sig_eps_sq) + 1L]] <- sig_eps_sq
				real(sig_eps_sq, y)
			},
			.package = "fetwfe"
		)
		# The estimated variances differ, so handing over sig_eps_c_sq would show.
		expect_false(
			isTRUE(all.equal(fit$sig_eps_sq, fit$sig_eps_c_sq)),
			info = est
		)
		expect_length(rec$sig_eps_sq, 1L)
		expect_identical(rec$sig_eps_sq, list(fit$sig_eps_sq), info = est)
	}
})

# The p >= NT regime, where `debiasedATT()` and the bootstrap band refit the
# q = 1 nuisance. Each block builds a high-dimensional fit per scale, so it
# skips on CRAN, as the file that owns this fixture does.
.rse428_hd_sim <- function() {
	simulateData(
		genCoefs(
			G = 3,
			T = 5,
			d = 20,
			density = 0.08,
			eff_size = 6,
			seed = 11
		),
		N = 60,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = 1001
	)
}

.rse428_hd_band_quantities <- list(
	estimate = .rse428_q(1, function(x) x$ci$estimate),
	pointwise_ci_low = .rse428_q(1, function(x) x$ci$pointwise_ci_low),
	pointwise_ci_high = .rse428_q(1, function(x) x$ci$pointwise_ci_high)
)

test_that("high dimensions, gls = TRUE with supplied variances: the fit, debiasedATT() and the bootstrap band are equivariant (#428)", {
	skip_on_cran()
	sim <- .rse428_hd_sim()
	at <- .rse428_fits(
		function(k) {
			.rse428_scaled_fit(
				fetwfe,
				k,
				sim = sim,
				sig_eps_sq = 0.5 * k^2,
				sig_eps_c_sq = 0.5 * k^2,
				ci_type = "pointwise"
			)
		},
		.RSE428_SCALES
	)
	ref <- at(1)
	expect_gte(ref$p, ref$N * ref$T)
	.rse428_expect_live(ref)
	.rse428_expect_nuisance_live(ref)
	.rse428_expect_equivariant(
		at,
		.rse428_fit_quantities(variances = FALSE),
		.RSE428_SCALES
	)
	.rse428_expect_equivariant(
		function(k) debiasedATT(at(k)),
		list(
			att = .rse428_q(1, function(x) x$att),
			se = .rse428_q(1, function(x) x$se),
			var_reg = .rse428_q(2, function(x) x$var_reg),
			var_weight = .rse428_q(2, function(x) x$var_weight)
		),
		.RSE428_SCALES
	)
	.rse428_expect_equivariant(
		function(k) simultaneousCIs(at(k), method = "bootstrap", seed = 1),
		.rse428_hd_band_quantities,
		.RSE428_SCALES
	)
})

test_that("high dimensions, gls = FALSE: debiasedATT() and the bootstrap band are equivariant (#428)", {
	skip_on_cran()
	sim <- .rse428_hd_sim()
	at <- .rse428_fits(
		function(k) .rse428_scaled_fit(fetwfe, k, sim = sim, gls = FALSE),
		.RSE428_SCALES
	)
	ref <- at(1)
	expect_gte(ref$p, ref$N * ref$T)
	.rse428_expect_nuisance_live(ref)
	# Not debiasedATT()'s se: its V2 term reads the fused fit, which keeps the
	# unstandardized grid on this route (#490).
	.rse428_expect_equivariant(
		function(k) debiasedATT(at(k)),
		list(
			att = .rse428_q(1, function(x) x$att),
			var_reg = .rse428_q(2, function(x) x$var_reg)
		),
		.RSE428_SCALES
	)
	.rse428_expect_equivariant(
		function(k) simultaneousCIs(at(k), method = "bootstrap", seed = 1),
		.rse428_hd_band_quantities,
		.RSE428_SCALES
	)
})

# The other fit options, one block each.
.rse428_no_never_treated <- function() {
	sim <- .rse428_sim
	pd <- sim$pdata
	ever <- unique(pd[[sim$unit_var]][pd[[sim$treatment]] == 1])
	sim$pdata <- pd[pd[[sim$unit_var]] %in% ever, ]
	sim
}

# The fits warn at every scale that the panel was truncated. Capture rather
# than expect: at edition 2, expect_warning() would also swallow a second,
# unexpected warning.
.rse428_fit_no_never_treated <- function(k, sim) {
	ws <- testthat::capture_warnings(
		fit <- .rse428_scaled_fit(
			fetwfe,
			k,
			sim = sim,
			allow_no_never_treated = TRUE,
			ci_type = "pointwise"
		)
	)
	expect_length(ws, 1L)
	expect_true(grepl("No never-treated units in input data", ws, fixed = TRUE))
	fit
}

.rse428_option_rows <- list(
	`fetwfe(se_type = "cluster")` = function(k) {
		.rse428_scaled_fit(
			fetwfe,
			k,
			se_type = "cluster",
			ci_type = "pointwise"
		)
	},
	`fetwfe(se_type = "conservative")` = function(k) {
		.rse428_scaled_fit(
			fetwfe,
			k,
			se_type = "conservative",
			ci_type = "pointwise"
		)
	},
	`fetwfe(fusion_structure = "event_study")` = function(k) {
		.rse428_scaled_fit(
			fetwfe,
			k,
			fusion_structure = "event_study",
			ci_type = "pointwise"
		)
	},
	`fetwfe(indep_counts)` = function(k) {
		.rse428_scaled_fit(
			fetwfe,
			k,
			indep_counts = .rse428_sim$indep_counts,
			ci_type = "pointwise"
		)
	},
	`fetwfe(allow_no_never_treated = TRUE)` = function(k) {
		.rse428_fit_no_never_treated(k, .rse428_no_never_treated())
	},
	`fetwfe(fusion_matrix = <identity>)` = function(k) {
		# The identity of the default fit's treatment block.
		n_treat <- length(.rse428_shared("pointwise", 1)$treat_inds)
		.rse428_scaled_fit(
			fetwfe,
			k,
			fusion_matrix = diag(n_treat),
			ci_type = "pointwise"
		)
	},
	`betwfe(se_type = "cluster")` = function(k) {
		.rse428_scaled_fit(
			betwfe,
			k,
			se_type = "cluster",
			ci_type = "pointwise"
		)
	}
)
for (.rse428_row in names(.rse428_option_rows)) {
	test_that(sprintf("%s is equivariant (#428)", .rse428_row), {
		at <- .rse428_fits(
			.rse428_option_rows[[.rse428_row]],
			.RSE428_SCALES
		)
		.rse428_expect_live(at(1))
		.rse428_expect_equivariant(
			at,
			.rse428_fit_quantities(),
			.RSE428_SCALES
		)
	})
}
