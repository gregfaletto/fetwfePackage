# Call-site policy test for the shared singular-Gram branch (#429 item 4).
#
# .recompute_gram_and_sandwich() carries the package's one genuine inference
# policy split in its `on_singular` argument: when the recomputed Gram on the
# selected support is singular, eventStudy() and cohortTimeATTs() DEGRADE
# (all-NA standard errors) while simultaneousCIs() STOPs. Until now that split
# was asserted only by test-recompute-gram-sandwich-400.R, which calls the
# helper directly with a hand-built rank-deficient matrix -- nothing pinned it
# where a user actually observes it, at the three public accessors.
#
# The way in is the scenario the method-entry validators exist for.
# R/class_helpers.R says they re-run as "defense-in-depth; the object may have
# been hand-modified between construction and method call". So: take a real
# fit and overwrite one SELECTED column of X_final with another selected
# column. The selected support becomes rank-deficient while every slot the
# validators inspect is untouched, so the object stays valid and the branch
# becomes reachable. See the roxygen on .recompute_gram_and_sandwich() in
# R/variance_machinery.R for the canonical statement of what is reachable from
# where.
#
# cohortStudy() is deliberately NOT driven here: its `se` is a FIT-TIME value
# read straight out of the stored catt_df slot (R/cohort_study.R recomputes
# nothing), so modifying X_final after construction cannot change it.

.cs429_dat <- local({
	cf <- genCoefs(G = 3, T = 5, d = 2, density = 0.5, eff_size = 2, seed = 429)
	simulateData(cf, N = 80, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 429)
})

# q = 0.5 is load-bearing on the two fetwfe cases: at q = 1 (and q = 2) the
# UNMODIFIED fit already has calc_ses = FALSE, and then every degrade assertion
# below would pass while measuring nothing. That is what the PRECOND blocks
# exist to catch.
.cs429_cases <- list(
	list(
		label = "fetwfe/default",
		fit = function() {
			fetwfeWithSimulatedData(
				.cs429_dat,
				q = 0.5,
				se_type = "default",
				verbose = FALSE
			)
		}
	),
	list(
		label = "fetwfe/cluster",
		fit = function() {
			fetwfeWithSimulatedData(
				.cs429_dat,
				q = 0.5,
				se_type = "cluster",
				verbose = FALSE
			)
		}
	),
	# etwfe performs no selection, so it drives the scalar-NA "all features"
	# sentinel branch of the scaffold rather than an index vector.
	list(
		label = "etwfe",
		fit = function() etwfeWithSimulatedData(.cs429_dat, verbose = FALSE)
	)
)

# X_final lives under $internal for fetwfe and at top level for etwfe/betwfe.
.cs429_calc_ses <- function(x) {
	if (inherits(x, "fetwfe")) x$internal$calc_ses else x$calc_ses
}
.cs429_get_X <- function(x) {
	if (inherits(x, "fetwfe")) x$internal$X_final else x$X_final
}
.cs429_set_X <- function(x, M) {
	if (inherits(x, "fetwfe")) {
		x$internal$X_final <- M
	} else {
		x$X_final <- M
	}
	x
}
# The selected support, derived the way .resolve_selected_support() derives it
# for the two cases this file uses: the nonzero theta entries after the
# intercept for `fetwfe`, and every column for `etwfe` (whose support is the NA
# sentinel). A wrong answer here fails safe -- it would leave the Gram
# invertible and turn the degrade assertions red rather than green.
#
# The branch is `fetwfe` vs everything-else, NOT bridge vs OLS. `betwfe` is a
# bridge estimator but its class is "betwfe" alone -- it does not inherit
# "fetwfe" -- so it would take the every-column branch and get the wrong
# support. Adding a `betwfe` case means changing this predicate first; the
# natural form is `if (!is.null(x$internal$theta_hat))`.
.cs429_sel <- function(x) {
	if (inherits(x, "fetwfe")) {
		which(x$internal$theta_hat[2:(x$p + 1)] != 0)
	} else {
		seq_len(ncol(x$X_final))
	}
}

for (.cs429_case in .cs429_cases) {
	# Precondition, on the UNMODIFIED fit. Kept in its own test_that() so a
	# failure is attributable: without it the POLICY block below would pass
	# vacuously on any fit that reported no standard errors to begin with.
	test_that(
		paste0(
			"PRECOND ",
			.cs429_case$label,
			": the unmodified fit reports real SEs (#429)"
		),
		{
			fit <- .cs429_case$fit()
			expect_true(isTRUE(.cs429_calc_ses(fit)))
			expect_true(any(is.finite(eventStudy(fit)$se)))
		}
	)

	test_that(
		paste0(
			"POLICY ",
			.cs429_case$label,
			": a singular selected-support Gram degrades at eventStudy()/cohortTimeATTs() and stops at simultaneousCIs() (#429)"
		),
		{
			fit <- .cs429_case$fit()
			sel <- .cs429_sel(fit)
			X_mod <- .cs429_get_X(fit)
			X_mod[, sel[2]] <- X_mod[, sel[1]]
			fit_mod <- .cs429_set_X(fit, X_mod)

			# 1. The modified object is still a VALID fit -- no validator
			# inspects the rank of X_final -- so what follows is the accessors'
			# policy on a legitimate object, not on a broken one.
			expect_no_error(fetwfe:::.assert_estimator_object(fit_mod))

			# 2-3. Degrade. The warning assertion is NOT hygiene: it is what
			# attributes the all-NA SEs to the singular-Gram branch
			# specifically, rather than to the calc_ses = FALSE path that Part 3
			# of test-cross-accessor-scaffold-guardrail-400.R already covers. Do
			# not replace it with suppressWarnings(). expect_warning() returns
			# its expression's value, so each accessor is called exactly once.
			# The `length(x) > 0 &&` clause matters because
			# all(is.na(numeric(0))) is TRUE.
			es_se <- expect_warning(eventStudy(fit_mod), "not invertible")$se
			expect_true(length(es_se) > 0 && all(is.na(es_se)))
			cta_se <- expect_warning(
				cohortTimeATTs(fit_mod),
				"not invertible"
			)$se
			expect_true(length(cta_se) > 0 && all(is.na(cta_se)))

			# 4. Stop, with EXACTLY the shared constant. The nesting is
			# prescribed and the obvious simplifications fail on correct code:
			# expect_warning() OUTSIDE and tryCatch() INSIDE is the only order
			# that works (suppressWarnings() inside makes the warning
			# expectation fail; with no tryCatch the error escapes and aborts
			# the block).
			#
			# tryCatch rather than expect_error() because the caught condition
			# is reused by the two content anchors below -- one call, three
			# assertions. Note expect_error() would need fixed = TRUE if it were
			# used: its pattern defaults to a regex and this message contains
			# "(p >= NT)", so the default form does not match its own subject.
			# Measured: default FAILS, fixed = TRUE PASSES.
			e <- expect_warning(
				tryCatch(
					simultaneousCIs(fit_mod, family = "cohort"),
					error = function(e) e
				),
				"not invertible"
			)
			expect_s3_class(e, "error")
			# conditionMessage() errors on a non-condition, which would abort
			# the block before the two content anchors ran -- same shape and
			# same reason as the length guard above.
			emsg <- if (inherits(e, "condition")) {
				conditionMessage(e)
			} else {
				NA_character_
			}
			expect_identical(emsg, fetwfe:::.SINGULAR_GRAM_ANALYTIC_STOP_MSG)
			expect_true(grepl("not invertible", emsg, fixed = TRUE))
			expect_true(grepl("method = 'bootstrap'", emsg, fixed = TRUE))
		}
	)
}
