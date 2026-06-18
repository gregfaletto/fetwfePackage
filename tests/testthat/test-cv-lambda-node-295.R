# Tests for CV-selecting the high-dimensional debiasing constant lambda_c (#295).
#
# The high-dim (p >= NT) nodewise penalty is
#   lambda_node = lambda_c * max(|a|) * sqrt(log(p) / N).
# Historically lambda_c was a fixed constant (default the theory scale 1.0). #295
# adds `lambda_c = "cv"`, which selects the constant per fit by cross-validating
# the unit-level Riesz loss 0.5 v'Sig v - a'v over the KKT-feasible region (the
# desparsified-lasso / auto-DML standard), with a fall back to 1.0 when no grid
# penalty is feasible. The same CV (on the overall-ATT direction) serves both the
# debiasedATT() point estimate and every simultaneousCIs() band effect, so one
# constant is shared and the selection is deterministic.

# ---- A fast synthetic high-dim case for direct .cv_lambda_node() tests --------
# p > n so Sig = X'X/n is singular (the regime the nodewise solve exists for).
.cv295_synth <- function() {
	set.seed(42)
	N_units <- 12L
	T <- 5L
	n <- N_units * T
	p <- 80L
	X <- matrix(stats::rnorm(n * p), n, p)
	list(
		Sig = crossprod(X) / n,
		a = c(stats::rnorm(5), rep(0, p - 5)) * 0.3,
		X = X,
		N_units = N_units,
		T = T
	)
}

# ============================ helper-level tests =============================

test_that(".cv_lambda_node selects a feasible grid constant", {
	s <- .cv295_synth()
	r <- .cv_lambda_node(s$Sig, s$a, s$X, s$N_units, s$T)
	expect_false(r$fallback)
	# The selection is one of the grid multipliers ...
	expect_true(r$lambda_c %in% r$mult_grid)
	# ... and it is one that actually passed the KKT-feasibility gate.
	expect_true(r$feasible[match(r$lambda_c, r$mult_grid)])
})

test_that(".cv_lambda_node only scores feasible grid points", {
	s <- .cv295_synth()
	r <- .cv_lambda_node(s$Sig, s$a, s$X, s$N_units, s$T)
	# A CV loss exists for exactly the feasible grid points (NA elsewhere): the
	# minimiser is taken only over feasible candidates.
	expect_true(all(is.na(r$cv_loss[!r$feasible])))
	expect_true(all(!is.na(r$cv_loss[r$feasible])))
	# The picked multiplier minimises the loss over the feasible set.
	feas <- which(r$feasible)
	expect_identical(
		r$lambda_c,
		r$mult_grid[feas[which.min(r$cv_loss[feas])]]
	)
})

test_that(".cv_lambda_node is deterministic and preserves the ambient RNG", {
	s <- .cv295_synth()
	set.seed(7L)
	state0 <- .Random.seed
	r1 <- .cv_lambda_node(s$Sig, s$a, s$X, s$N_units, s$T)
	# The data-seeded fold draw is restored: ambient .Random.seed is untouched.
	expect_identical(.Random.seed, state0)
	# Re-running (even from a different ambient seed) gives the same selection.
	set.seed(123L)
	r2 <- .cv_lambda_node(s$Sig, s$a, s$X, s$N_units, s$T)
	expect_identical(r1$lambda_c, r2$lambda_c)
	expect_identical(r1$cv_loss, r2$cv_loss)
})

test_that(".cv_lambda_node falls back to the theory scale when nothing is feasible", {
	s <- .cv295_synth()
	# A grid of vanishingly small penalties: the solver cannot reach KKT
	# feasibility (||Sig v - a||_inf <= lambda) for any of them in the gate budget.
	r <- .cv_lambda_node(
		s$Sig,
		s$a,
		s$X,
		s$N_units,
		s$T,
		mult_grid = c(1e-8, 1e-7)
	)
	expect_true(r$fallback)
	expect_identical(r$lambda_c, 1.0)
	expect_true(all(!r$feasible))
	expect_true(all(is.na(r$cv_loss)))
})

# ===================== end-to-end tests on a real p >= NT fit =================
# A genuine high-dim fetwfe fit (N=20, T=8, d=12 => p=376 >> NT=160). The
# simulator enforces N >= (G+1)(d+1), so build the raw panel directly and supply
# the noise variances to bypass REML (mirrors test-debiased-att-highdim.R).
.cv295_panel <- function() {
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

.cv295_fit <- function(gls = TRUE) {
	args <- list(
		pdata = .cv295_panel(),
		time_var = "year",
		unit_var = "unit",
		treatment = "treat",
		covs = paste0("x", 1:12),
		response = "y",
		q = 0.5,
		verbose = FALSE
	)
	if (gls) {
		args$sig_eps_sq <- 0.16
		args$sig_eps_c_sq <- 0.1
	} else {
		args$gls <- FALSE
	}
	do.call(fetwfe, args)
}

hd295 <- .cv295_fit(gls = TRUE)
hd295_nogls <- .cv295_fit(gls = FALSE)

test_that("debiasedATT(lambda_c = 'cv') selects and reports the CV constant", {
	db <- debiasedATT(hd295, lambda_c = "cv")
	expect_identical(db$lambda_c_selection, "cv")
	expect_true(is.numeric(db$lambda_c) && db$lambda_c > 0)
	expect_true(is.finite(db$att))
	expect_true(is.finite(db$se) && db$se > 0)
	# The CV diagnostics ride along only on the "cv" path.
	expect_false(is.null(db$lambda_cv))
	expect_false(db$lambda_cv$fallback)
	expect_true(db$lambda_c %in% db$lambda_cv$mult_grid)
})

test_that("debiasedATT defaults to the fixed theory scale (cv is opt-in)", {
	db <- debiasedATT(hd295)
	expect_identical(db$lambda_c_selection, "fixed")
	expect_identical(db$lambda_c, 1.0)
	# No CV diagnostics when the constant is fixed.
	expect_null(db$lambda_cv)
})

test_that("debiasedATT(lambda_c = 'cv') works on a gls = FALSE (real-data) fit", {
	db <- debiasedATT(hd295_nogls, lambda_c = "cv")
	expect_identical(db$lambda_c_selection, "cv")
	expect_true(is.finite(db$att))
	expect_true(is.finite(db$se) && db$se > 0)
})

test_that("debiasedATT is deterministic at lambda_c = 'cv'", {
	a <- debiasedATT(hd295, lambda_c = "cv")
	b <- debiasedATT(hd295, lambda_c = "cv")
	expect_identical(a$lambda_c, b$lambda_c)
	expect_identical(a$att, b$att)
	expect_identical(a$se, b$se)
})

test_that("point estimate and bootstrap band share one CV constant (#295)", {
	# The fixture's cohorts adopt in consecutive periods, so the band's overall-ATT
	# direction equals debiasedATT()'s exactly; CVing it with the same data seed
	# resolves to the identical constant -- one constant for the point estimate and
	# every band effect (the auto-DML convention). (For scattered adoption times
	# high-dim debiasedATT() is unsupported; the band alone still CVs deterministically.)
	db <- debiasedATT(hd295, lambda_c = "cv")
	sc <- simultaneousCIs(
		hd295,
		family = "cohort",
		method = "bootstrap",
		B = 100L,
		seed = 1L,
		lambda_c = "cv"
	)
	expect_identical(sc$lambda_c_selection, "cv")
	expect_identical(sc$lambda_c, db$lambda_c)
})

test_that("the fixed lambda_c escape hatch is unchanged by the cv feature", {
	# Selecting the (formerly only) fixed constant explicitly still reports
	# "fixed" and reproduces the default path byte-for-byte.
	d_default <- debiasedATT(hd295)
	d_fixed1 <- debiasedATT(hd295, lambda_c = 1.0)
	expect_identical(d_fixed1$lambda_c_selection, "fixed")
	expect_identical(d_default$att, d_fixed1$att)
	expect_identical(d_default$se, d_fixed1$se)
})

# ============================== validation ===================================

test_that("lambda_c rejects anything but a positive number or 'cv'", {
	expect_error(debiasedATT(hd295, lambda_c = "nope"))
	expect_error(debiasedATT(hd295, lambda_c = -1))
	expect_error(debiasedATT(hd295, lambda_c = 0))
	expect_error(debiasedATT(hd295, lambda_c = c(1, 2)))
	expect_error(
		simultaneousCIs(
			hd295,
			method = "bootstrap",
			lambda_c = "nope",
			B = 50L
		)
	)
	expect_error(
		simultaneousCIs(
			hd295,
			method = "bootstrap",
			lambda_c = -1,
			B = 50L
		)
	)
	# Validated regardless of method (symmetric with debiasedATT), even though
	# lambda_c is only USED on the high-dim bootstrap path: a garbage value under
	# method = "analytic" should error rather than be silently ignored.
	expect_error(
		simultaneousCIs(hd295, method = "analytic", lambda_c = "garbage")
	)
})
