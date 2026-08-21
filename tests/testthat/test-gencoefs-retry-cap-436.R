# Issue #436: genCoefsCore()'s rejection-sampling loop over `density` had no
# upper bound on its attempts. A legal-but-degenerate `density` (e.g. 1e-300,
# which .validate_gen_coefs_scalars() accepts) made the call spin forever with
# no message, no error, and an apparently frozen session. The loop is now
# capped at 10,000 attempts and errors with an actionable message instead.
#
# THREE THINGS CARRY THE MUTATION POWER HERE, and each is written the way it is
# because a weaker version was measured not to work:
#
#   (1) THE LITERAL IS `in 10000 attempts`, WITH fixed = TRUE. Not "10000":
#       grepl("10000", "... in 100000 attempts ...", fixed = TRUE) is TRUE
#       (measured), so a cap mutated to 100000L would survive it. Not
#       "density" either -- the validator's own message for the SAME argument
#       contains that word, so a mutation that made the validator reject
#       1e-300 before the loop is entered would leave this green with the cap
#       never exercised. That is issue #458's defect class.
#
#   (2) THE TIMEOUT GUARD. Without it, "cap removed", "counter never
#       incremented" and "stop() unreachable" all HANG the suite -- putting
#       back into the suite the exact pathology #436 exists to remove. 60 s is
#       generous against a correct path measured at ~0.02 s, so runner jitter
#       cannot false-positive it, and the on.exit() is required to stop the
#       limit leaking into later blocks (measured). ONE setTimeLimit() CALL
#       GUARDS EXACTLY ONE HANGING CALL: measured, a fired elapsed limit is
#       disarmed, and the next long call in the same block runs unguarded to
#       completion. Any block below with two calls that could hang therefore
#       arms it twice.
#
#   (3) THE DRAW-COUNT PIN. `.apply_seed(NULL)` is a verified no-op on the RNG
#       (it returns before touching it) and everything between it and the loop
#       is RNG-free, so the ambient stream advances by exactly the loop's
#       draws. That pins the cap VALUE, the off-by-one in `tries > max_tries`,
#       AND that each retry is still `rbinom(n = p, size = 1, prob = density)`
#       with its arguments unchanged -- 10,000 iterations of stream, where the
#       #332 byte-identity guardrail in test-gencoefs-targeted-sparsity-332.R
#       only ever exercises one, so it cannot see a retry-counter defect.
#
# THERE IS DELIBERATELY NO WALL-CLOCK ASSERTION. It cannot observe the
# condition it would name -- if the cap is absent the call never returns, so
# the timer never completes -- and a load-dependent assertion on a gating
# six-job matrix whose measured job times swing 80% run to run is a flake
# generator. The timing evidence lives in the ExecPlan instead.

test_that("a degenerate density errors instead of hanging (#436)", {
	setTimeLimit(elapsed = 60, transient = TRUE)
	on.exit(setTimeLimit(), add = TRUE)

	expect_error(
		genCoefsCore(G = 3, T = 5, d = 2, density = 1e-300, eff_size = 1),
		"in 10000 attempts",
		fixed = TRUE
	)
})

test_that("the retry-cap error is identical from both doors and sets no call (#436)", {
	# ONE setTimeLimit() GUARDS EXACTLY ONE HANGING CALL, so this block arms it
	# twice -- once per door. Measured: inside a single top-level task (which is
	# what a test_that() block is under devtools::test()), an elapsed limit
	# fires at its deadline and is then DISARMED; a second long call in the same
	# block runs unguarded to completion. Without the re-arm, a cap-removing
	# mutation hangs this block forever on the second door, which is the exact
	# failure mode the guard exists to prevent. The on.exit() still clears the
	# limit at block end so it cannot leak into later blocks.
	on.exit(setTimeLimit(), add = TRUE)

	# Degrade to a readable value rather than raising, so a mutation that makes
	# the call SUCCEED fails these assertions cleanly instead of aborting the
	# block (same reasoning as test-gen-coefs-validator-429.R).
	setTimeLimit(elapsed = 60, transient = TRUE)
	e_core <- tryCatch(
		genCoefsCore(G = 3, T = 5, d = 2, density = 1e-300, eff_size = 1),
		error = function(e) e
	)
	setTimeLimit(elapsed = 60, transient = TRUE)
	e_gc <- tryCatch(
		genCoefs(G = 3, T = 5, d = 2, density = 1e-300, eff_size = 1),
		error = function(e) e
	)

	expect_true(grepl("in 10000 attempts", msg_of(e_core), fixed = TRUE))

	# Door parity: one error reaches both public doors with IDENTICAL text.
	# That is only true because the message names no function (Decision Log
	# D5) -- `call. = FALSE` suppresses conditionCall(), so a genCoefsCore()
	# prefix would tell a genCoefs() caller about a function they never
	# called.
	expect_identical(msg_of(e_gc), msg_of(e_core))

	# The `call. = FALSE` contract from #431 item 1, which
	# test-gen-coefs-validator-429.R pins on every one of its cells.
	expect_null(call_of(e_core))
	expect_null(call_of(e_gc))
})

test_that("the capped loop makes exactly 10,000 draws, no more and no fewer (#436)", {
	setTimeLimit(elapsed = 60, transient = TRUE)
	on.exit(setTimeLimit(), add = TRUE)

	# p = getP(G = 3, T = 5, d = 2) = 50, and every retry draws
	# rbinom(n = p, size = 1, prob = density). Comparing `.Random.seed` is an
	# integer-state comparison, so none of the floating-point tolerance
	# hazards apply.
	seed_now <- function() get(".Random.seed", envir = globalenv())

	set.seed(1)
	invisible(tryCatch(
		genCoefsCore(G = 3, T = 5, d = 2, density = 1e-300, eff_size = 1),
		error = function(e) NULL
	))
	after <- seed_now()

	set.seed(1)
	for (i in seq_len(10000L)) {
		invisible(rbinom(n = 50, size = 1, prob = 1e-300))
	}
	# Measured TRUE against exactly 10,000 reference draws, and FALSE against
	# 9,999, 10,001 and the 10x mutant.
	expect_identical(after, seed_now())
})

test_that("a workable density still returns with the cap in place (#436)", {
	# density = 0.01 at p = 50 needs ~2.5 attempts, four orders of magnitude
	# from the 10,000 budget. This is a "the cap does not fire when it should
	# not" control and NOTHING more -- it proves nothing about an off-by-one at
	# the boundary, which is what the draw-count pin above is for. Seeded so it
	# is reproducible.
	res <- genCoefsCore(
		G = 3,
		T = 5,
		d = 2,
		density = 0.01,
		eff_size = 1,
		seed = 1
	)

	expect_identical(names(res), c("beta", "theta"))
	expect_true(sum(res$theta != 0) > 0)
})
