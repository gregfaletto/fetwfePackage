# Issue #429 item 9: all SIXTEEN atomic conditions of
# .validate_gen_coefs_scalars() (R/gen_coefs.R:1200-1246), through BOTH public
# doors.
#
# The validator is the single-sourced scalar gate that #401 folded out of
# genCoefs() and genCoefsCore(). Ten of its sixteen atomic conditions had no
# message-anchored test, and the whole `G` message literal was unpinned, so
# weakening most of them was invisible to the two files that claim to cover this
# surface.
#
# Each cell asserts three things, and all three are needed:
#   (1) the error message contains the condition's own literal, `fixed = TRUE`;
#   (2) the two doors produce IDENTICAL messages (this is what shows the cell
#       reached the single-sourced validator rather than some door-local check);
#   (3) conditionCall() is NULL -- the `call. = FALSE` contract from #431 item 1.
#
# WHY THE MESSAGE LITERAL IS WHAT CARRIES THE POWER, not expect_error() alone.
# Both doors carry THREE consecutive redundant backstops --
# `stopifnot(G >= 1); stopifnot(T >= 2); stopifnot(G <= T - 1)` at
# R/gen_coefs.R:474-476 and :918-920 -- so three of the sixteen cells (G < 1,
# T < 2, G > T - 1) error SOMEWHERE regardless of what the validator does. A
# bare expect_error() on those would be satisfied by the backstop and the
# mutant would walk.
#
# WHY THE FIXTURES CARRY NO OTHER ARGUMENTS. The two doors run their prologues
# in different orders: genCoefsCore() runs .resolve_cohort_count_arg() (:898),
# match.arg(fusion_structure) (:900) and .validate_targeted_sparsity() (:906)
# BEFORE the scalar validator (:910), while genCoefs() runs the scalar validator
# (:375) BEFORE match.arg() (:383). Keep every argument that genCoefsCore()
# validates before :910 out of these fixtures -- `fusion_structure`,
# `n_signal_cohorts`, `treat_base_levels` and the deprecated `R` alias -- or a
# fixture bad in two ways at once will produce different messages from the two
# doors and break assertion (2) for a reason unrelated to this item.
#
# This file SUPERSEDES nothing in test-genCoefs.R / test-genCoefsCore.R. Six of
# these cells overlap with the ten assertions there, but those use bare
# `regexp =` (and one, "density must be", is a truncated prefix), so the cells
# here are strictly stronger; the older ones stay as door-level smoke tests.

.gcv429_base <- list(G = 3, T = 5, d = 2, density = 0.5, eff_size = 1)

.gcv429_T_msg <- "T must be a numeric value greater than or equal to 2"
.gcv429_G_msg <-
	"G must be a numeric value greater than or equal to 1 (at least one treated cohort)"
.gcv429_GT_msg <- "G must be less than or equal to T - 1"
.gcv429_d_msg <- "d must be a non-negative numeric value"
.gcv429_dens_msg <- "density must be numeric, greater than 0 and at most 1"
.gcv429_eff_msg <- "eff_size must be a numeric value"

# One row per atomic condition, in source order.
.gcv429_cells <- list(
	list(id = "T not numeric", over = list(T = "5"), lit = .gcv429_T_msg),
	list(id = "T length != 1", over = list(T = c(5, 6)), lit = .gcv429_T_msg),
	list(id = "T < 2", over = list(T = 1), lit = .gcv429_T_msg),
	list(id = "G not numeric", over = list(G = "3"), lit = .gcv429_G_msg),
	list(id = "G length != 1", over = list(G = c(2, 3)), lit = .gcv429_G_msg),
	list(id = "G < 1", over = list(G = 0), lit = .gcv429_G_msg),
	list(id = "G > T - 1", over = list(G = 5, T = 5), lit = .gcv429_GT_msg),
	list(id = "d not numeric", over = list(d = "2"), lit = .gcv429_d_msg),
	list(id = "d length != 1", over = list(d = c(1, 2)), lit = .gcv429_d_msg),
	list(id = "d < 0", over = list(d = -1), lit = .gcv429_d_msg),
	list(
		id = "density not numeric",
		over = list(density = "0.5"),
		lit = .gcv429_dens_msg
	),
	list(
		id = "density length != 1",
		over = list(density = c(0.4, 0.5)),
		lit = .gcv429_dens_msg
	),
	list(id = "density <= 0", over = list(density = 0), lit = .gcv429_dens_msg),
	list(
		id = "density > 1",
		over = list(density = 1.1),
		lit = .gcv429_dens_msg
	),
	list(
		id = "eff_size not numeric",
		over = list(eff_size = "1"),
		lit = .gcv429_eff_msg
	),
	list(
		id = "eff_size length != 1",
		over = list(eff_size = c(1, 2)),
		lit = .gcv429_eff_msg
	)
)

# ------------------------------------------------------------------------------
# SANITY BLOCK, first in the file. Without it, a change that made EVERY input
# fail would leave all sixteen negative cells green.
# ------------------------------------------------------------------------------
test_that("the unmodified base arguments succeed through both doors (#429 item 9)", {
	gc <- do.call(genCoefs, .gcv429_base)
	expect_s3_class(gc, "FETWFE_coefs")

	core <- do.call(genCoefsCore, .gcv429_base)
	expect_true(is.list(core))
	expect_identical(names(core), c("beta", "theta"))
})

test_that(".validate_gen_coefs_scalars() pins all sixteen conditions on both doors (#429 item 9)", {
	expect_length(.gcv429_cells, 16L)

	for (cell in .gcv429_cells) {
		args <- utils::modifyList(.gcv429_base, cell$over)

		e_gc <- expect_error(
			do.call(genCoefs, args),
			cell$lit,
			fixed = TRUE,
			info = paste("genCoefs:", cell$id)
		)
		e_core <- expect_error(
			do.call(genCoefsCore, args),
			cell$lit,
			fixed = TRUE,
			info = paste("genCoefsCore:", cell$id)
		)

		expect_identical(
			conditionMessage(e_gc),
			conditionMessage(e_core),
			info = cell$id
		)
		expect_null(conditionCall(e_gc), info = paste("genCoefs:", cell$id))
		expect_null(
			conditionCall(e_core),
			info = paste("genCoefsCore:", cell$id)
		)
	}
})
