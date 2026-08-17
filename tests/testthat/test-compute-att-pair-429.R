# Issue #429 item 5: a direct unit test of .compute_att_pair()
# (R/variance_machinery.R:1036-1082), the helper the #400/#401 consolidation
# created by folding the in-sample / independent ATT pair out of every
# estimator. It has no direct unit test of its own contract: the three call
# sites at tests/testthat/test-small-fixes-431.R:108-221 are all ERROR-path
# fixtures that exercise .call_te()'s naming, never the helper's field mapping,
# its assert_att_se guard, or its precheck ordering.
#
# Four named mutations of R/variance_machinery.R:1036-1082 must redden this
# file:
#   (1) map `in_sample_att_se_no_prob` to `att_te_se` instead of
#       `att_te_se_no_prob`                                    -> block (a)
#   (2) drop the `assert_att_se` guard                          -> block (b)
#   (3) drop the `indep_precheck()` call                        -> block (c)
#   (4) move that call AFTER the independent .call_te()         -> block (c)
#
# Every fixture is a hand-built argument list driving the real
# getTeResultsOLS(). No simulation, no fitting, nothing random.
#
# .compute_att_pair()'s first parameter is `te_fn_name`, a CHARACTER SCALAR and
# not a function object (renamed in #431 item 2; .call_te() resolves the name in
# the package namespace). Calling shape copied from
# tests/testthat/test-small-fixes-431.R:130-141.

# ------------------------------------------------------------------------------
# Shared fixture. Deliberately degenerate in the arithmetic, not in the code
# path: psi_mat and gram_inv are both the 2x2 identity, which makes att_var_1
# hand-derivable (see block (a)) without changing which branches run.
# ------------------------------------------------------------------------------
.cap429_base_args <- list(
	sig_eps_sq = 1,
	N = 50L,
	T = 3L,
	G = 2L,
	num_treats = 2L,
	cohort_tes = c(1, 2),
	psi_mat = matrix(c(1, 0, 0, 1), nrow = 2L, ncol = 2L),
	gram_inv = diag(2),
	tes = c(0.5, 1.5),
	calc_ses = TRUE
)

.cap429_in_sample_probs <- list(
	cohort_probs = c(0.5, 0.5),
	cohort_probs_overall = c(0.2, 0.2)
)

# ------------------------------------------------------------------------------
# (a) The five in-sample fields land in the right result slots.
# ------------------------------------------------------------------------------
test_that(".compute_att_pair() maps the in-sample fields correctly (#429 item 5)", {
	out <- fetwfe:::.compute_att_pair(
		"getTeResultsOLS",
		.cap429_base_args,
		.cap429_in_sample_probs,
		list(),
		FALSE
	)

	expect_identical(
		names(out),
		c(
			"in_sample_att_hat",
			"in_sample_att_se",
			"in_sample_att_se_no_prob",
			"in_sample_att_var_1",
			"in_sample_att_var_2",
			"indep_att_hat",
			"indep_att_se",
			"indep_att_var_1",
			"indep_att_var_2"
		)
	)

	# THE assertion mutation (1) has to walk past: in_sample_att_se_no_prob must
	# be the no-probability SE, i.e. sqrt(att_var_1). Both sides here are the
	# same double produced by the same code path (R/variance_machinery.R
	# computes att_te_se_no_prob <- sqrt(att_var_1) from the same local), so
	# exactness is legitimate under .workflow/PROFILE.md section 12 gotcha 8.
	expect_identical(
		out$in_sample_att_se_no_prob,
		sqrt(out$in_sample_att_var_1)
	)

	# NON-VACUITY CONTROL, and it is load-bearing rather than decorative: the
	# two SEs COINCIDE when att_var_2 is zero, and mutation (1) would then be
	# invisible. Measured on this fixture: att_var_1 = 0.003333, att_var_2 =
	# 0.0125, se_no_prob = 0.057735, att_se = 0.125831.
	expect_true(out$in_sample_att_var_2 > 0)
	expect_false(isTRUE(all.equal(
		out$in_sample_att_se,
		out$in_sample_att_se_no_prob
	)))

	# Absolute, hand-derived reference. With gram_inv and psi_mat both the
	# identity, att_var_1 collapses to
	#     sig_eps_sq * sum(cohort_probs^2) / (N * T) = 1 * 0.5 / 150 = 1/300,
	# so this literal is derived from the formula rather than recovered from the
	# function's own output. expect_equal() at the default tolerance rather than
	# expect_identical(): the two sides are NOT the same arithmetic -- one is a
	# BLAS t(psi) %*% gram_inv %*% psi chain, the other a scalar literal
	# division -- and .workflow/PROFILE.md section 12 gotcha 8 records that
	# exactly this shape passed macOS/Windows and failed all four Linux jobs on
	# #427's first cross-platform run. (It happens to be bit-identical on
	# macOS/Accelerate, which is precisely the evidence that misleads.)
	expect_equal(out$in_sample_att_var_1, 1 / 300)

	# With no independent count data, all four independent fields are NA.
	expect_true(is.na(out$indep_att_hat))
	expect_true(is.na(out$indep_att_se))
	expect_true(is.na(out$indep_att_var_1))
	expect_true(is.na(out$indep_att_var_2))
})

# ------------------------------------------------------------------------------
# (b) assert_att_se is what throws.
# ------------------------------------------------------------------------------
test_that(".compute_att_pair()'s assert_att_se guard fires on an NA SE (#429 item 5)", {
	no_se_args <- .cap429_base_args
	no_se_args$calc_ses <- FALSE

	# The anchor is the FIELD NAME, not base R's full stopifnot() deparse
	# ("!is.na(out$in_sample_att_se) is not TRUE"). That format belongs to base
	# R, not to this package, and CI runs R-devel. The field name still
	# discriminates: without the guard there is no error at all, and any other
	# error would not name this field.
	expect_error(
		fetwfe:::.compute_att_pair(
			"getTeResultsOLS",
			no_se_args,
			.cap429_in_sample_probs,
			list(),
			FALSE,
			assert_att_se = TRUE
		),
		"in_sample_att_se",
		fixed = TRUE
	)

	# CONTROL that makes the assertion above non-vacuous: the identical call
	# with the guard off returns normally, with the NA the guard rejects. So the
	# error is attributable to assert_att_se and to nothing else in the fixture.
	out <- fetwfe:::.compute_att_pair(
		"getTeResultsOLS",
		no_se_args,
		.cap429_in_sample_probs,
		list(),
		FALSE,
		assert_att_se = FALSE
	)
	expect_true(is.na(out$in_sample_att_se))
})

# ------------------------------------------------------------------------------
# (c) indep_precheck runs exactly once, and BEFORE the independent .call_te().
# ------------------------------------------------------------------------------
test_that(".compute_att_pair() runs indep_precheck once, before the independent call (#429 item 5)", {
	indep_probs_args <- list(
		cohort_probs = c(0.4, 0.6),
		cohort_probs_overall = c(0.2, 0.2)
	)

	# --- called exactly once when independent count data is available --------
	n_calls <- 0L
	out <- fetwfe:::.compute_att_pair(
		"getTeResultsOLS",
		.cap429_base_args,
		.cap429_in_sample_probs,
		indep_probs_args,
		TRUE,
		indep_precheck = function() n_calls <<- n_calls + 1L
	)
	expect_identical(n_calls, 1L)
	expect_true(is.finite(out$indep_att_hat))
	# The independent branch really ran on the independent probabilities:
	# c(1, 2) %*% c(0.4, 0.6) = 1.6, against the in-sample c(0.5, 0.5) -> 1.5.
	expect_equal(out$indep_att_hat, 1.6)
	expect_equal(out$in_sample_att_hat, 1.5)

	# --- CONTROL: not called at all when there is no independent count data --
	n_calls_off <- 0L
	out_off <- fetwfe:::.compute_att_pair(
		"getTeResultsOLS",
		.cap429_base_args,
		.cap429_in_sample_probs,
		indep_probs_args,
		FALSE,
		indep_precheck = function() n_calls_off <<- n_calls_off + 1L
	)
	expect_identical(n_calls_off, 0L)
	expect_true(is.na(out_off$indep_att_hat))
	expect_true(is.na(out_off$indep_att_se))
	expect_true(is.na(out_off$indep_att_var_1))
	expect_true(is.na(out_off$indep_att_var_2))

	# --- ORDERING, and this is the half that tests ordering AS ordering ------
	# Poison the independent probabilities: a length-3 cohort_probs against
	# length-2 cohort_tes. (NA-poisoning does NOT work here -- the
	# stopifnot(all(!is.na(cohort_probs))) lives in getTeResults2(), not in
	# getTeResultsOLS().)
	poisoned <- list(
		cohort_probs = c(0.3, 0.3, 0.4),
		cohort_probs_overall = c(0.2, 0.2)
	)

	# Half one: with NO precheck, the independent call itself raises.
	expect_error(
		fetwfe:::.compute_att_pair(
			"getTeResultsOLS",
			.cap429_base_args,
			.cap429_in_sample_probs,
			poisoned,
			TRUE
		),
		"non-conformable arguments",
		fixed = TRUE
	)

	# Half two: with a precheck that stops, the precheck's error is the one that
	# surfaces -- so the precheck DISPLACED an error the independent call would
	# otherwise have raised, which is what "runs before" means operationally.
	# Both halves are required; either alone shows only that something errored.
	n_calls_order <- 0L
	expect_error(
		fetwfe:::.compute_att_pair(
			"getTeResultsOLS",
			.cap429_base_args,
			.cap429_in_sample_probs,
			poisoned,
			TRUE,
			indep_precheck = function() {
				n_calls_order <<- n_calls_order + 1L
				stop("PRECHECK-RAN-FIRST")
			}
		),
		"PRECHECK-RAN-FIRST",
		fixed = TRUE
	)
	expect_identical(n_calls_order, 1L)
})
