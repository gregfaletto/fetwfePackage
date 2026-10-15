# Cross-validation folds drawn over a label-free unit order (#524).
#
# `.cv_unit_order()` orders the units by their data; `.cv_row_folds()` and
# `.cv_unit_folds()` assign a seeded draw over that order. The first part tests
# those helpers directly; the second checks that fits built on them do not
# change when the units are relabeled. `.cv_unit_order()`'s roxygen describes
# the key.

# ---- Helpers ------------------------------------------------------------------

# A constant one-column design: every unit's first design row ties, so the
# response keys alone decide the order.
.ord524_const_X <- function(N, T) {
	matrix(1, N * T, 1L)
}

# The rows of a unit-major panel, in the sequence of the units listed in `units`.
.ord524_rows <- function(units, T) {
	rep((units - 1L) * T, each = T) + rep.int(seq_len(T), length(units))
}

# The data at each rank of `ord`: each unit's response block and first design
# row, in rank order.
.ord524_at_rank <- function(ord, X, y, T) {
	list(
		y = if (is.null(y)) NULL else y[.ord524_rows(ord, T)],
		X = X[(ord - 1L) * T + 1L, , drop = FALSE]
	)
}

# ---- .cv_unit_order() ---------------------------------------------------------

test_that(".cv_unit_order() permutes with the units: the data at each rank is unchanged (#524)", {
	set.seed(5240)
	T <- 3L

	# A continuous response and design.
	N <- 30L
	X <- matrix(stats::rnorm(N * T * 2L), N * T, 2L)
	y <- stats::rnorm(N * T)
	o <- .cv_unit_order(X, N, T, y)
	o_x <- .cv_unit_order(X, N, T)
	for (k in 1:5) {
		rows <- .ord524_rows(sample(N), T)
		Xp <- X[rows, , drop = FALSE]
		yp <- y[rows]
		expect_identical(
			.ord524_at_rank(.cv_unit_order(Xp, N, T, yp), Xp, yp, T),
			.ord524_at_rank(o, X, y, T)
		)
		# The design key alone, as `.cv_lambda_node()` uses it.
		expect_identical(
			.ord524_at_rank(.cv_unit_order(Xp, N, T), Xp, NULL, T),
			.ord524_at_rank(o_x, X, NULL, T)
		)
	}

	# A count response on a two-valued, time-invariant design, with one unit
	# repeating another exactly: response ties that the design breaks, and a
	# complete tie. The two identical units can trade places under a
	# permutation, so the assertion is on the data at each rank.
	N <- 40L
	X <- matrix(rep(sample(0:1, N, replace = TRUE), each = T), N * T, 1L)
	y <- stats::rpois(N * T, 3)
	y[.ord524_rows(N, T)] <- y[.ord524_rows(1L, T)]
	X[.ord524_rows(N, T), ] <- X[.ord524_rows(1L, T), ]
	o <- .cv_unit_order(X, N, T, y)
	for (k in 1:5) {
		rows <- .ord524_rows(sample(N), T)
		Xp <- X[rows, , drop = FALSE]
		yp <- y[rows]
		expect_identical(
			.ord524_at_rank(.cv_unit_order(Xp, N, T, yp), Xp, yp, T),
			.ord524_at_rank(o, X, y, T)
		)
	}
})

test_that(".cv_unit_order() is unchanged by a positive rescaling or a shift of the response (#524)", {
	set.seed(5241)
	N <- 30L
	T <- 4L
	X <- matrix(stats::rnorm(N * T), N * T, 1L)
	y <- stats::rnorm(N * T)
	o <- .cv_unit_order(X, N, T, y)
	expect_identical(.cv_unit_order(X, N, T, 3 * y + 100), o)
	expect_identical(.cv_unit_order(X, N, T, 7.1 * y - 2), o)

	# Counts on a constant design, where the order rests on the response keys
	# and their tie-breaks alone.
	Xc <- .ord524_const_X(N, T)
	yc <- stats::rpois(N * T, 2)
	oc <- .cv_unit_order(Xc, N, T, yc)
	expect_identical(.cv_unit_order(Xc, N, T, 3 * yc + 100), oc)
	expect_identical(.cv_unit_order(Xc, N, T, 7.1 * yc - 2), oc)
})

test_that(".cv_unit_order() is unchanged by negating a response with no exact mirror twins (#524)", {
	# T = 3, a constant design. Units 1 and 2 agree in absolute deviation in
	# every period and differ in sign in the last one, without being exact
	# negatives of each other; unit 3 makes the mean exactly 0 without adding an
	# exact mirror pair.
	y <- c(0, 1, 2, 0, 1, -2, 0, -2, 0)
	X <- .ord524_const_X(3L, 3L)
	expect_identical(
		.cv_unit_order(X, 3L, 3L, -y),
		.cv_unit_order(X, 3L, 3L, y)
	)
})

test_that(".cv_unit_order() rounds each key to 10 digits of its scale (#524)", {
	X <- .ord524_const_X(4L, 2L)

	# Below the quantum: the units' period-1 values differ by 2^-50, which the
	# rounding absorbs, so period 2 decides: the `u1` pair sorts ahead of the
	# `u2` pair. Within each exact mirror pair the signed key puts the negative
	# twin first. The four units sum to exactly 0.
	u1 <- c(0.75 + 2^-50, -0.25)
	u2 <- c(0.75, -0.5)
	expect_identical(
		.cv_unit_order(X, 4L, 2L, c(u1, u2, -u1, -u2)),
		c(3L, 1L, 4L, 2L)
	)

	# Above the quantum: a period-1 gap of 2^-31 (about 4.7e-10) survives the
	# rounding, so period 1 decides and the `u2` pair sorts first.
	u1 <- c(1, -0.25)
	u2 <- c(1 - 2^-31, -0.5)
	expect_identical(
		.cv_unit_order(X, 4L, 2L, c(u1, u2, -u1, -u2)),
		c(4L, 2L, 3L, 1L)
	)
})

test_that(".cv_unit_order() breaks ties by the design, then the sign-normalized and signed deviations (#524)", {
	# The design key: the two units tie on the absolute key, and the unit with
	# the smaller first design row comes first, whichever way round they are
	# given (their signed keys would order them the other way).
	Xd <- matrix(c(2, 2, 1, 1), 4L, 1L)
	expect_identical(.cv_unit_order(Xd, 2L, 2L, c(-1, 1, 1, -1)), c(2L, 1L))
	expect_identical(
		.cv_unit_order(
			Xd[c(3, 4, 1, 2), , drop = FALSE],
			2L,
			2L,
			c(1, -1, -1, 1)
		),
		c(1L, 2L)
	)

	# The sign-normalized key: units 1 and 2 tie on the absolute and design
	# keys and differ in sign in the last period only; unit 3 balances the mean.
	X3 <- .ord524_const_X(3L, 3L)
	y <- c(0, -1, 2, 0, -1, -2, 0, 2, 0)
	expect_identical(.cv_unit_order(X3, 3L, 3L, y), c(1L, 2L, 3L))
	expect_identical(
		.cv_unit_order(X3, 3L, 3L, y[c(4:6, 1:3, 7:9)]),
		c(2L, 1L, 3L)
	)

	# The signed key: exact mirror twins with the same design; the one whose
	# first deviation is negative comes first, whichever way round they are given.
	X2 <- .ord524_const_X(2L, 2L)
	expect_identical(.cv_unit_order(X2, 2L, 2L, c(1, -2, -1, 2)), c(2L, 1L))
	expect_identical(.cv_unit_order(X2, 2L, 2L, c(-1, 2, 1, -2)), c(1L, 2L))

	# A complete tie keeps the input order: units 1 and 2 are identical, and
	# unit 3 sorts ahead of both.
	expect_identical(
		.cv_unit_order(.ord524_const_X(3L, 2L), 3L, 2L, c(1, -1, 1, -1, 0, 0)),
		c(3L, 1L, 2L)
	)
})

# ---- .cv_row_folds() and .cv_unit_folds() --------------------------------------

test_that(".cv_row_folds() and .cv_unit_folds() assign the seeded draw over the order (#524)", {
	N <- 6L
	T <- 3L
	K <- 4L
	n <- N * T + 2L # two rows past N * T, as ridge pseudo-rows would be
	draw_rows <- function(seed) {
		set.seed(seed)
		f <- sample(1:n %% K)
		f[f == 0] <- K
		f
	}
	draw_units <- function(seed) {
		set.seed(seed)
		sample(rep(seq_len(K), length.out = N))
	}

	# The identity order returns the draw itself.
	set.seed(5242)
	expect_identical(
		.cv_row_folds(n = n, N = N, T = T, K = K, ord = seq_len(N)),
		draw_rows(5242)
	)
	set.seed(5242)
	expect_identical(.cv_unit_folds(N, K, seq_len(N)), draw_units(5242))

	# A cycle, which is not its own inverse: the r-th block of T draws goes to
	# the rows of unit ord[r], and the draws past N * T to the rows past N * T.
	ord <- c(2:N, 1L)
	f <- draw_rows(5243)
	expected <- f
	for (r in seq_len(N)) {
		expected[(ord[r] - 1L) * T + seq_len(T)] <- f[(r - 1L) * T + seq_len(T)]
	}
	# The inverse mapping gives a different vector on this draw, so the
	# comparison below can tell the two apart.
	expect_false(identical(
		expected,
		f[c(.ord524_rows(ord, T), (N * T + 1L):n)]
	))
	set.seed(5243)
	folds <- .cv_row_folds(n = n, N = N, T = T, K = K, ord = ord)
	expect_identical(folds, expected)
	expect_setequal(folds, seq_len(K))
	expect_identical(folds[(N * T + 1L):n], f[(N * T + 1L):n])

	d <- draw_units(5244)
	expected <- d
	for (r in seq_len(N)) {
		expected[ord[r]] <- d[r]
	}
	expect_false(identical(expected, d[ord]))
	set.seed(5244)
	folds <- .cv_unit_folds(N, K, ord)
	expect_identical(folds, expected)
	expect_setequal(folds, seq_len(K))
})
