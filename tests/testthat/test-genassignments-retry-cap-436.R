# Issue #436: genAssignments() drew rmultinom(1, N, rep(1/(G+1), G+1)) until
# every group cleared its bar -- 1 unit normally, `d + 1` under
# guarantee_rank_condition = TRUE -- with no upper bound on attempts. Its
# callers assert `N >= G + 1` and `N >= (G + 1) * (d + 1)`, which guarantees
# only that the bar is ACHIEVABLE, not that it is achievable in practice: at
# N = G + 1 every group must get exactly one unit, so
# P = (G+1)! / (G+1)^(G+1), which is 4.8e-13 at G = 30. Both floors are
# reachable from a plain simulateData() call and both used to hang. The loop is
# now capped at 1,000,000 attempts.
#
# EVERY simulateData() CALL HERE PASSES `seed = NA` WITH A `set.seed(k)`
# IMMEDIATELY BEFORE IT, and both halves are mandatory. simulateData() warns
# only when `seed` is NULL; `seed = NA` is silent by design and means "use the
# ambient stream", so the preceding set.seed() is what supplies
# reproducibility. A blanket suppressWarnings() would instead make the suite's
# WARN-count watch unfalsifiable on exactly the blocks where a new warning is
# most likely, and `seed = NA` WITHOUT the set.seed() would leave the boundary
# block below drawing an unseeded geometric with mean 7,148 -- a varying
# runtime that someone would later "fix" by passing `seed = 1`, silently
# re-introducing the warning this rule exists to prevent.
#
# PER DECISION LOG D8 THIS FILE CARRIES NO DRAW-COUNT PIN, unlike its
# genCoefsCore() sibling. Replaying 1,000,000 rmultinom draws to pin the
# boundary costs 0.85-1.32 s on every run of a 9-minute suite across six CI
# jobs, and what it would buy is detection of a `>` vs `>=` off-by-one that
# moves the budget by one draw in a million and has no user-visible
# consequence whatsoever. An ABSENT cap is caught by the timeout guard, and a
# WRONG cap value by the exact literal `in 1000000 attempts`.

.gar436_coefs_plain <- genCoefs(
	G = 30,
	T = 31,
	d = 1,
	density = 0.2,
	eff_size = 1,
	seed = 1L
)

.gar436_coefs_grc <- genCoefs(
	G = 12,
	T = 13,
	d = 4,
	density = 0.2,
	eff_size = 1,
	seed = 2L
)

.gar436_coefs_boundary <- genCoefs(
	G = 10,
	T = 11,
	d = 1,
	density = 0.2,
	eff_size = 1,
	seed = 3L
)

test_that("an infeasible N = G + 1 panel errors instead of hanging (#436)", {
	setTimeLimit(elapsed = 60, transient = TRUE)
	on.exit(setTimeLimit(), add = TRUE)

	# N = G + 1 is exactly the boundary genCohortTimeFE()'s
	# stopifnot(N >= G + 1) admits, and P(success) = 31!/31^31 = 4.8e-13 here.
	# The literal is `in 1000000 attempts` with fixed = TRUE: a bare "1000000"
	# would be a substring of a 10x-mutated cap's message.
	set.seed(1)
	expect_error(
		simulateData(
			.gar436_coefs_plain,
			N = 31,
			sig_eps_sq = 0.5,
			sig_eps_c_sq = 0.5,
			seed = NA
		),
		"in 1000000 attempts",
		fixed = TRUE
	)
})

test_that("an infeasible guarantee_rank_condition panel errors instead of hanging (#436)", {
	setTimeLimit(elapsed = 60, transient = TRUE)
	on.exit(setTimeLimit(), add = TRUE)

	# N = (G + 1) * (d + 1) = 65 is the boundary simulateDataCore()'s own
	# stopifnot() admits, and it is a MUCH harsher regime than the bar-1 floor:
	# P = 2.7e-6 at G = 8, d = 4. This is also the first thing in the package to
	# exercise guarantee_rank_condition = TRUE inside genAssignments() at all.
	set.seed(2)
	e <- tryCatch(
		simulateData(
			.gar436_coefs_grc,
			N = 65,
			sig_eps_sq = 0.5,
			sig_eps_c_sq = 0.5,
			guarantee_rank_condition = TRUE,
			seed = NA
		),
		error = function(e) e
	)

	expect_true(grepl("in 1000000 attempts", msg_of(e), fixed = TRUE))

	# TWO literals, each asserted on the thing it names, because all three of
	# the obvious single choices are defective:
	#   * "d + 1" would blind the two assertions in
	#     test-simulate-data-highdim-293.R that match the bare pattern
	#     "d \\+ 1" to prove the up-front stopifnot() fired -- a mutation that
	#     today HANGS would become a pass;
	#   * a bare "25%" matches a `%%`-escaped message and the correct one
	#     alike, and observes `suggestion` rather than the `bar_desc` this
	#     block exists to exercise;
	#   * "units in 1000000 attempts" discriminates from the other branch only
	#     by the plural `s` (the plain branch says "at least one unit"), which
	#     a future grammar fix would silently delete -- note d = 0 under
	#     guarantee_rank_condition = TRUE renders "at least 1 units".
	expect_true(grepl("at least 5 units", msg_of(e), fixed = TRUE))
	expect_true(grepl("about 25% is", msg_of(e), fixed = TRUE))

	# The `call. = FALSE` contract from #431 item 1.
	expect_null(call_of(e))
})

test_that("a feasible boundary panel still succeeds with the cap in place (#436)", {
	setTimeLimit(elapsed = 60, transient = TRUE)
	on.exit(setTimeLimit(), add = TRUE)

	# G = 10, N = G + 1: exact P = 1.399e-4, E[draws] = 7,148 -- feasible, and
	# the only thing here exercising the SUCCESS path near the boundary. It
	# cannot flake (P(false error) = 1.7e-61 at cap 1e6). It is NOT evidence
	# about "a cap set too low": its power against a 100000L cap is 8.4e-7, and
	# the cap value is pinned deterministically by the exact literal above.
	set.seed(3)
	sim <- simulateData(
		.gar436_coefs_boundary,
		N = 11,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = NA
	)

	expect_s3_class(sim, "FETWFE_simulated")
	expect_equal(nrow(sim$pdata), 11 * 11)
})

test_that("a comfortable panel returns promptly in both bar modes (#436)", {
	setTimeLimit(elapsed = 60, transient = TRUE)
	on.exit(setTimeLimit(), add = TRUE)

	# The controls from the reproduction: the same coefficient objects at an
	# N the multinomial clears easily. Without these, a mutation that made
	# genAssignments() error unconditionally would still satisfy every
	# expect_error() above.
	set.seed(4)
	sim_plain <- simulateData(
		.gar436_coefs_plain,
		N = 200,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = NA
	)
	expect_s3_class(sim_plain, "FETWFE_simulated")
	expect_equal(nrow(sim_plain$pdata), 200 * 31)

	# Re-armed because this block makes a SECOND call that a hang-mutation
	# could stall: a fired elapsed limit disarms itself (measured; see the
	# header of test-gencoefs-retry-cap-436.R).
	setTimeLimit(elapsed = 60, transient = TRUE)
	set.seed(5)
	sim_grc <- simulateData(
		.gar436_coefs_grc,
		N = 300,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		guarantee_rank_condition = TRUE,
		seed = NA
	)
	expect_s3_class(sim_grc, "FETWFE_simulated")
	expect_equal(nrow(sim_grc$pdata), 300 * 13)
})
