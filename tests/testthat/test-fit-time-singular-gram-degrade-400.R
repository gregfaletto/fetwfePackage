library(testthat)
library(fetwfe)

# The FIT-TIME singular-Gram degrade in getCohortATTsFinal() (via
# .recompute_gram_and_sandwich(), #400) is REACHABLE through an ordinary public
# fit. It is the backstop R/utility.R's #395 note relies on: "A genuinely
# singular collapsed design is caught downstream: getGramInv() degrades calc_ses
# to FALSE."
#
# The reachable window exists because the two rank tests disagree. The #395 gate
# runs anyNA(coef(lm(y ~ . + 0, df))) on the UNCENTERED, column-scaled design at
# lm.fit's 1e-7 QR tolerance; getGramInv() tests the CENTERED Gram at
# max(dim) * .Machine$double.eps relative. Covariates collinear at 1e-6 sit
# between them.
#
# This test exists because the roxygen on .recompute_gram_and_sandwich() twice
# asserted the opposite (that the branch was defensive-only / unreachable). An
# untested claim about reachability rots; this one pins it.

.singular_gram_panel <- function(near_dup) {
	set.seed(3)
	N <- 40L
	T <- 5L
	cohorts <- sample(c(0, 3, 4, 5), N, replace = TRUE)
	df <- expand.grid(t = seq_len(T), id = seq_len(N))
	df$cohort <- cohorts[df$id]
	df$id <- as.character(df$id)
	df$changed <- as.integer(df$cohort > 0 & df$t >= df$cohort)
	df$x1 <- stats::rnorm(nrow(df))
	# Near-collinear, NOT lm()-aliased at near_dup = 1e-6.
	df$x2 <- df$x1 + near_dup * stats::rnorm(nrow(df))
	df$y <- 0.5 * df$x1 + 0.3 * df$changed + stats::rnorm(nrow(df), sd = 0.4)
	df
}

test_that("a public etwfe() fit reaches the fit-time singular-Gram degrade (#400)", {
	df <- .singular_gram_panel(1e-6)

	expect_warning(
		fit <- etwfe(
			pdata = df,
			time_var = "t",
			unit_var = "id",
			treatment = "changed",
			covs = c("x1", "x2"),
			response = "y",
			se_type = "cluster"
		),
		"not invertible"
	)

	# The point estimate survives; only inference degrades.
	expect_false(fit$calc_ses)
	expect_true(is.na(fit$att_se))
	expect_true(is.finite(fit$att_hat))
})

test_that("tighter collinearity is intercepted by the #395 rank gate instead (#400)", {
	# The complement of the case above: at 1e-8 the design IS lm()-aliased, so
	# #395 errors before the Gram is ever formed. Pins the ordering of the two
	# tolerances -- a future change to either threshold moves this boundary, and
	# without this nothing would notice.
	df <- .singular_gram_panel(1e-8)

	expect_error(
		etwfe(
			pdata = df,
			time_var = "t",
			unit_var = "id",
			treatment = "changed",
			covs = c("x1", "x2"),
			response = "y",
			se_type = "cluster"
		)
	)
})
