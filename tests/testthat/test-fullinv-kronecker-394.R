# Tests for #394: genFullInvFusionTransformMat() must agree with the
# authoritative main-path untransformCoefImproved(). The builder is a standalone
# block-diagonal D_N^{-1}; its three covariate-interaction blocks (cohort x X,
# time x X, treatment x X) must use the package's LEVEL-major Kronecker order
# `block %x% I_d`, not the feature-major `I_d %x% block`. The bug (feature-major)
# made `fetwfe(add_ridge = TRUE)` with d >= 2 penalize coordinate-permuted
# coefficients. This transform-agreement invariant is exactly what would have
# caught it (it fails by O(1) under the wrong order).

.check_fullinv_agrees <- function(G, T, d, fs, seed) {
	nt <- fetwfe:::getNumTreats(G, T)
	fi <- fetwfe:::getFirstInds(G, T)
	p <- G + (T - 1L) + d + d * G + d * (T - 1L) + nt + d * nt
	set.seed(seed)
	theta <- stats::rnorm(p)
	Dfull <- as.matrix(fetwfe:::genFullInvFusionTransformMat(
		fi,
		T,
		G,
		d,
		nt,
		fusion_structure = fs
	))
	via_D <- as.numeric(Dfull %*% theta)
	via_fn <- as.numeric(fetwfe:::untransformCoefImproved(
		theta,
		T,
		G,
		p,
		d,
		nt,
		first_inds = fi,
		fusion_structure = fs
	))
	# Byte-exact on the dev platform (0/1 arithmetic); the 1e-12 tolerance guards
	# cross-platform summation-order ULP. The bug is O(1), so it bites regardless.
	expect_equal(via_D, via_fn, tolerance = 1e-12)
}

test_that("genFullInvFusionTransformMat() == untransformCoefImproved() for d>=2, both fusion structures (#394)", {
	for (cfg in list(c(3L, 6L, 2L), c(2L, 4L, 3L), c(4L, 7L, 4L))) {
		for (fs in c("cohort", "event_study")) {
			.check_fullinv_agrees(cfg[1], cfg[2], cfg[3], fs, seed = 394)
		}
	}
})

test_that("genFullInvFusionTransformMat() agrees with untransform at d in {0,1} too (low-d guard) (#394)", {
	# At d in {0,1} the two Kronecker orders coincide, so this passed pre-fix as
	# well; it is an invariant guard against a future edit breaking the low-d path.
	.check_fullinv_agrees(3L, 6L, 1L, "cohort", seed = 1)
	.check_fullinv_agrees(3L, 6L, 0L, "cohort", seed = 1)
})
