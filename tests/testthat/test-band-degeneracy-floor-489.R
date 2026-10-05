# The simultaneous band's degeneracy rule, `.band_nondegenerate()` (#489): the
# rule itself, the units in which the bootstrap applies its reference, and a
# custom contrast whose variance is rounding error.

test_that(".band_nondegenerate() floors its relative tolerance at the reference (#489)", {
	nondeg <- fetwfe:::.band_nondegenerate
	tol <- sqrt(.Machine$double.eps)
	# Exact zeros are degenerate.
	expect_identical(nondeg(c(0, 0, 0), 1), c(FALSE, FALSE, FALSE))
	expect_identical(nondeg(c(1, 0), 1), c(TRUE, FALSE))
	# Residue-sized variances are degenerate against a real reference variance.
	expect_identical(nondeg(c(1e-35, 2e-35), 1e-3), c(FALSE, FALSE))
	# With the reference below the family's largest variance, the tolerance is
	# relative to that variance.
	big <- 1e-4
	expect_identical(
		nondeg(c(big, 0.5 * tol * big, 2 * tol * big), 1e-3 * big),
		c(TRUE, FALSE, TRUE)
	)
	# Multiplying both arguments by k^2 leaves the classification of a real
	# variance, a residue and a zero unchanged.
	for (k in c(1e-9, 1e-6, 1e-3, 1, 1e3, 1e6)) {
		expect_identical(
			nondeg(c(0.3, 1e-20, 0) * k^2, 0.01 * k^2),
			c(TRUE, FALSE, FALSE),
			info = sprintf("k = %g", k)
		)
	}
})

test_that(".simultaneous_bootstrap_crit() floors degeneracy at var_ref in variance units (#489)", {
	N <- 60L
	n <- N * 4L
	var_ref <- 0.5 / n
	tol <- sqrt(.Machine$double.eps)
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
			var_ref = var_ref
		)$nondeg
	}
	expect_false(nondeg_at(0.5 * tol * var_ref))
	expect_true(nondeg_at(2 * tol * var_ref))
})

# A fixture on which the bridge fuses cells to one value. A custom contrast
# whose weights sum to zero across fused cells then has a variance that is
# rounding error rather than exactly 0.
.bdf489_sim <- simulateData(
	genCoefs(G = 3, T = 5, d = 2, density = 0.15, eff_size = 2, seed = 3),
	N = 150,
	sig_eps_sq = 0.5,
	sig_eps_c_sq = 0.5,
	seed = 1003
)

# `fetwfe()` on `.bdf489_sim` with its response multiplied by `k`; `...` is
# passed to `fetwfe()`.
.bdf489_fit <- function(k, ...) {
	sim <- .bdf489_sim
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
