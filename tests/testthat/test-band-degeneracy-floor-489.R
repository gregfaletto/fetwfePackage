# The simultaneous band's degeneracy rule, `.band_nondegenerate()`, and the
# routes that apply it (#489, #501).

test_that(".band_nondegenerate()'s tolerance is eps times the larger of the family's largest variance per unit of squared weight and the reference (#489)", {
	nondeg <- fetwfe:::.band_nondegenerate
	tol <- .Machine$double.eps
	# Exact zeros are degenerate.
	expect_identical(nondeg(c(0, 0, 0), 1, 1), c(FALSE, FALSE, FALSE))
	expect_identical(nondeg(c(1, 0), 1, 1), c(TRUE, FALSE))
	# Residue-sized variances are degenerate against a real reference variance.
	expect_identical(nondeg(c(1e-35, 2e-35), 1e-3, 1), c(FALSE, FALSE))
	# With the reference below the family's largest variance, the tolerance is
	# relative to that variance.
	big <- 1e-4
	expect_identical(
		nondeg(c(big, 0.5 * tol * big, 2 * tol * big), 1e-3 * big, 1),
		c(TRUE, FALSE, TRUE)
	)
	# Multiplying both arguments by k^2 leaves the classification of a real
	# variance, a residue and a zero unchanged.
	for (k in c(1e-9, 1e-6, 1e-3, 1, 1e3, 1e6)) {
		expect_identical(
			nondeg(c(0.3, 1e-20, 0) * k^2, 0.01 * k^2, 1),
			c(TRUE, FALSE, FALSE),
			info = sprintf("k = %g", k)
		)
	}
	# Scaling one effect's contrast row by s multiplies its variance and its
	# `w2` by s^2, and leaves every effect's classification unchanged.
	v <- c(0.3, 0.2, 1e-20, 0)
	for (j in seq_along(v)) {
		for (s in c(1e-9, 1e9)) {
			w2 <- replace(rep(1, 4), j, s^2)
			expect_identical(
				nondeg(v * w2, 0.01, w2),
				c(TRUE, TRUE, FALSE, FALSE),
				info = sprintf("effect %d scaled by %g", j, s)
			)
		}
	}
	# Two rows whose weights differ by 1e8 are both non-degenerate.
	expect_identical(
		nondeg(c(0.3, 0.2 * 1e-16), 0.01, c(1, 1e-16)),
		c(TRUE, TRUE)
	)
	# A zero row is degenerate.
	expect_identical(nondeg(c(0.3, 0), 0.01, c(1, 0)), c(TRUE, FALSE))
	# A scalar `w2` applies to every effect.
	expect_identical(nondeg(c(0.03, 0.035 * tol), 0.01, 4), c(TRUE, FALSE))
})

test_that(".simultaneous_bootstrap_crit() applies var_ref in variance units (#489)", {
	N <- 60L
	n <- N * 4L
	var_ref <- 0.5 / n
	tol <- .Machine$double.eps
	set.seed(489)
	u <- stats::rnorm(N)
	u <- u / sqrt(sum(u^2))
	# A one-column family whose variance colSums(F^2) / n^2 is v.
	nondeg_at <- function(v) {
		fetwfe:::.simultaneous_bootstrap_crit(
			matrix(u * n * sqrt(v), N, 1),
			n = n,
			alpha = 0.05,
			B = 10,
			seed = 1,
			var_ref = var_ref,
			w2 = 1
		)$nondeg
	}
	expect_false(nondeg_at(0.5 * tol * var_ref))
	expect_true(nondeg_at(2 * tol * var_ref))
})

# A fixture on which the bridge fuses cells to one value.
.bdf489_sim <- simulateData(
	genCoefs(G = 3, T = 5, d = 2, density = 0.15, eff_size = 2, seed = 3),
	N = 150,
	sig_eps_sq = 0.5,
	sig_eps_c_sq = 0.5,
	seed = 1003
)

# `fetwfe()` on `sim` with its response multiplied by `k`; `...` is passed to
# `fetwfe()`.
.bdf489_fit <- function(k, ..., sim = .bdf489_sim) {
	pdata <- sim$pdata
	pdata[[sim$response]] <- pdata[[sim$response]] * k
	fetwfe(
		pdata = pdata,
		time_var = sim$time_var,
		unit_var = sim$unit_var,
		treatment = sim$treatment,
		response = sim$response,
		covs = sim$covs,
		...
	)
}

test_that("a custom contrast whose variance is rounding error gets an NA adjusted p-value at every scale (#489)", {
	ref <- .bdf489_fit(1)
	est <- cohortTimeATTs(ref)$estimate
	# The first three cells that share one nonzero estimate.
	cells <- integer(0)
	for (i in which(est != 0)) {
		tied <- which(abs(est - est[i]) <= 1e-12 * abs(est[i]))
		if (length(tied) >= 3L) {
			cells <- tied[1:3]
			break
		}
	}
	expect_length(cells, 3L)
	contrast <- matrix(0, 1, length(ref$treat_inds))
	contrast[1, cells] <- c(0.1, 0.2, -0.3)
	bands <- list(
		analytic = list(
			se_type = "default",
			band = function(fit) {
				simultaneousCIs(
					fit,
					family = "custom",
					contrasts = contrast,
					method = "analytic"
				)
			}
		),
		bootstrap = list(
			se_type = "default",
			band = function(fit) {
				simultaneousCIs(
					fit,
					family = "custom",
					contrasts = contrast,
					method = "bootstrap",
					seed = 1
				)
			}
		),
		conservative = list(
			se_type = "conservative",
			band = function(fit) {
				suppressMessages(simultaneousCIs(
					fit,
					family = "custom",
					contrasts = contrast,
					method = "analytic"
				))
			}
		)
	)
	for (k in c(1, 1e-9, 1e-4, 1e6)) {
		fits <- list(
			default = .bdf489_fit(k),
			conservative = .bdf489_fit(k, se_type = "conservative")
		)
		for (nm in names(bands)) {
			info <- sprintf("%s at k = %g", nm, k)
			fit <- fits[[bands[[nm]]$se_type]]
			# The three cells still share one nonzero estimate.
			cell_est <- cohortTimeATTs(fit)$estimate[cells]
			expect_true(cell_est[1] != 0, info = info)
			expect_equal(cell_est, rep(cell_est[1], 3L), info = info)
			sc <- bands[[nm]]$band(fit)
			# The contrast's standard error is positive, so the classification
			# rests on the tolerance rather than on an exact zero.
			expect_gt(
				sc$ci$pointwise_ci_high - sc$ci$pointwise_ci_low,
				0,
				label = paste("pointwise interval width,", info)
			)
			expect_true(is.na(sc$adjusted_p_values), info = info)
		}
	}
})

# The default fixture of `test-response-scale-equivariance-428.R`.
.bdf489_sim428 <- simulateData(
	genCoefs(G = 3, T = 4, d = 2, density = 0.5, eff_size = 2, seed = 123),
	N = 120,
	sig_eps_sq = 0.5,
	sig_eps_c_sq = 0.5,
	seed = 456
)

test_that("a custom family's band does not depend on the scale of its rows' weights (#501)", {
	fits <- list(
		default = .bdf489_fit(1, sim = .bdf489_sim428),
		conservative = .bdf489_fit(
			1,
			se_type = "conservative",
			sim = .bdf489_sim428
		)
	)
	routes <- list(
		analytic = list(fit = "default", method = "analytic"),
		bootstrap = list(fit = "default", method = "bootstrap", seed = 1),
		conservative = list(fit = "conservative", method = "analytic")
	)
	band <- function(route, contrasts) {
		r <- routes[[route]]
		suppressMessages(simultaneousCIs(
			fits[[r$fit]],
			family = "custom",
			contrasts = contrasts,
			method = r$method,
			seed = r$seed
		))
	}
	# `C0` with `m` rows: a 1 on cell i in row i.
	num_treats <- length(fits$default$treat_inds)
	c0 <- function(m) {
		out <- matrix(0, m, num_treats)
		out[cbind(seq_len(m), seq_len(m))] <- 1
		out
	}
	ws <- list(
		rep(1e-9, 3),
		rep(1e-5, 3),
		rep(1e6, 3),
		rep(1e9, 3),
		c(1e-3, 1e-5),
		c(1, 1e-5),
		c(1, 1e-9),
		c(1, 1e-18),
		c(1e-3, 1e-3, 1e-5)
	)
	bounds <- c(
		"estimate",
		"simultaneous_ci_low",
		"simultaneous_ci_high",
		"pointwise_ci_low",
		"pointwise_ci_high"
	)
	for (route in names(routes)) {
		refs <- list()
		for (m in 2:3) {
			refs[[m]] <- band(route, c0(m))
			# `C0`'s band is live: more than one effect has an adjusted p-value.
			expect_gt(
				sum(is.finite(refs[[m]]$adjusted_p_values)),
				1L,
				label = sprintf(
					"finite adjusted p-values of C0 (%s, %d rows)",
					route,
					m
				)
			)
		}
		for (w in ws) {
			info <- sprintf(
				"%s, w = (%s)",
				route,
				paste(format(w), collapse = ", ")
			)
			ref <- refs[[length(w)]]
			sc <- band(route, diag(w, nrow = length(w)) %*% c0(length(w)))
			expect_equal(
				sc$critical_value,
				ref$critical_value,
				tolerance = 1e-6,
				info = info
			)
			expect_equal(
				sc$adjusted_p_values,
				ref$adjusted_p_values,
				tolerance = 1e-6,
				info = info
			)
			# Each bound, divided by its row's weight, is `C0`'s bound.
			for (b in bounds) {
				expect_equal(
					sc$ci[[b]] / w,
					ref$ci[[b]],
					tolerance = 1e-6,
					info = paste(info, b)
				)
			}
		}
	}
})

test_that("the rule receives the reference sig_eps_sq / (N * T), or var(y) / (N * T) from a fit without sig_eps_sq, the conservative branch's squared standard errors, and a custom family's squared row norms on the analytic, bootstrap, conservative and desparsified routes (#501)", {
	fit <- .bdf489_fit(1, sim = .bdf489_sim428)
	fit_cons <- .bdf489_fit(1, se_type = "conservative", sim = .bdf489_sim428)
	real <- fetwfe:::.band_nondegenerate
	rec <- new.env(parent = emptyenv())
	recorder <- function(v, v_ref, w2) {
		rec$calls[[length(rec$calls) + 1L]] <- list(
			v = v,
			v_ref = v_ref,
			w2 = w2
		)
		real(v, v_ref, w2)
	}
	nt <- fit$N * fit$T
	v_ref <- fit$sig_eps_sq / nt

	rec$calls <- list()
	testthat::with_mocked_bindings(
		simultaneousCIs(fit),
		.band_nondegenerate = recorder,
		.package = "fetwfe"
	)
	expect_length(rec$calls, 1L)
	expect_equal(rec$calls[[1]]$v_ref, v_ref)

	# On `.simultaneous_bootstrap_crit()`'s scale; see its `var_ref`.
	rec$calls <- list()
	testthat::with_mocked_bindings(
		simultaneousCIs(fit, method = "bootstrap", seed = 1),
		.band_nondegenerate = recorder,
		.package = "fetwfe"
	)
	expect_length(rec$calls, 1L)
	expect_equal(rec$calls[[1]]$v_ref, nt^2 * v_ref)

	rec$calls <- list()
	sc <- testthat::with_mocked_bindings(
		suppressMessages(simultaneousCIs(fit_cons)),
		.band_nondegenerate = recorder,
		.package = "fetwfe"
	)
	# The conservative branch built this band.
	expect_equal(sc$critical_value, sc$bonferroni_critical_value)
	expect_length(rec$calls, 1L)
	expect_equal(
		rec$calls[[1]]$v_ref,
		fit_cons$sig_eps_sq / (fit_cons$N * fit_cons$T)
	)
	expect_equal(
		rec$calls[[1]]$v,
		((sc$ci$pointwise_ci_high - sc$ci$pointwise_ci_low) /
			(2 * sc$pointwise_critical_value))^2
	)

	# A custom family whose rows' squared norms are 4 and 1.25, on each route.
	C <- matrix(0, 2, length(fit$treat_inds))
	C[1, 1] <- 2
	C[2, 2:3] <- c(1, 0.5)
	routes <- list(
		analytic = function() {
			simultaneousCIs(fit, family = "custom", contrasts = C)
		},
		bootstrap = function() {
			simultaneousCIs(
				fit,
				family = "custom",
				contrasts = C,
				method = "bootstrap",
				seed = 1
			)
		},
		conservative = function() {
			suppressMessages(
				simultaneousCIs(fit_cons, family = "custom", contrasts = C)
			)
		}
	)
	for (route in names(routes)) {
		rec$calls <- list()
		testthat::with_mocked_bindings(
			routes[[route]](),
			.band_nondegenerate = recorder,
			.package = "fetwfe"
		)
		expect_equal(
			lapply(rec$calls, "[[", "w2"),
			list(c(4, 1.25)),
			info = route
		)
	}

	# The rest builds a `gls = FALSE` fit on #428's high-dimensional fixture, so it
	# skips on CRAN, as that file's high-dimensional blocks do.
	skip_on_cran()
	sim_hd <- simulateData(
		genCoefs(G = 3, T = 5, d = 20, density = 0.08, eff_size = 6, seed = 11),
		N = 60,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = 1001
	)
	fit_hd <- .bdf489_fit(1, gls = FALSE, sim = sim_hd)
	nt_hd <- fit_hd$N * fit_hd$T
	y <- fit_hd$internal$y_final[seq_len(nt_hd)]
	rec$calls <- list()
	testthat::with_mocked_bindings(
		simultaneousCIs(fit_hd, method = "bootstrap", seed = 1),
		.band_nondegenerate = recorder,
		.package = "fetwfe"
	)
	expect_length(rec$calls, 1L)
	expect_equal(rec$calls[[1]]$v_ref, nt_hd^2 * stats::var(y) / nt_hd)

	# The custom family above, on the desparsified route.
	C_hd <- matrix(0, 2, length(fit_hd$treat_inds))
	C_hd[1, 1] <- 2
	C_hd[2, 2:3] <- c(1, 0.5)
	rec$calls <- list()
	testthat::with_mocked_bindings(
		simultaneousCIs(
			fit_hd,
			family = "custom",
			contrasts = C_hd,
			method = "bootstrap",
			seed = 1
		),
		.band_nondegenerate = recorder,
		.package = "fetwfe"
	)
	expect_equal(lapply(rec$calls, "[[", "w2"), list(c(4, 1.25)))
})
