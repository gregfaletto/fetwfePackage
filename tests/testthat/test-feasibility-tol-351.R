library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #351: the high-dim nodewise KKT feasibility check `.riesz_feasible()` floors
# its relative slack ABOVE coordinate-descent round-off, decoupled from the
# convergence control `riesz_tol`. The desparsified L1 direction binds its
# constraint (||Sig v - a||_inf = lambda_node) by construction, so a converged
# solver overshoots only by ~1e-9..1e-8 (relative); the old slack (riesz_tol /
# hardcoded 1e-9) sat below that and spuriously flagged binding directions as
# infeasible (gregfaletto/fetwfe#88). The fixed 1e-6 tolerance separates
# round-off from a genuine constraint violation (which overshoots by >> 1e-6).
# ------------------------------------------------------------------------------
test_that(".riesz_feasible() accepts round-off but rejects genuine infeasibility (#351)", {
	ln <- 0.19 # a representative lambda_node (the #88 anchor scale)

	# Binding direction, overshooting only by coordinate-descent round-off:
	# must read FEASIBLE under the 1e-6 tolerance.
	expect_true(fetwfe:::.riesz_feasible(ln, ln)) # exactly binding
	expect_true(fetwfe:::.riesz_feasible(ln * (1 + 1e-9), ln))
	expect_true(fetwfe:::.riesz_feasible(ln * (1 + 1e-8), ln))
	expect_true(fetwfe:::.riesz_feasible(ln * (1 + 1e-7), ln))

	# Genuine infeasibility (lambda_c too small / non-convergence) overshoots by
	# >> 1e-6: must read INFEASIBLE.
	expect_false(fetwfe:::.riesz_feasible(ln * (1 + 1e-3), ln))
	expect_false(fetwfe:::.riesz_feasible(ln * 17, ln)) # the test-142 magnitude

	# Vectorized over per-effect / per-cell directions.
	expect_identical(
		fetwfe:::.riesz_feasible(
			c(ln * (1 + 1e-8), ln * (1 + 1e-3)),
			c(ln, ln)
		),
		c(TRUE, FALSE)
	)

	# The bug, pinned: the old slack (riesz_tol = 1e-9) wrongly rejects a binding
	# direction overshooting by 1e-8 round-off; the 1e-6 default accepts it.
	expect_false(fetwfe:::.riesz_feasible(ln * (1 + 1e-8), ln, rtol = 1e-9))
	expect_true(fetwfe:::.riesz_feasible(ln * (1 + 1e-8), ln, rtol = 1e-6))
})
