# Issue #429 item 3: the three class validators must ROUTE through
# .check_c6_dims_toplevel().
#
# tests/testthat/test-small-fixes-431.R:467-509 already tests the helper
# DIRECTLY, with hand-built one-slot fixtures, so "helper body ->
# invisible(NULL)" is dead as a mutant. What survives is one level up: deleting
# the `.check_c6_dims_toplevel(x, cls)` call from any ONE of the three class
# validators (R/etwfe_class.R:147, R/betwfe_class.R:165,
# R/twfeCovs_class.R:158) leaves that direct unit test green. Nothing else
# covers it -- tests/testthat/test-class-validators.R:181-192 is the only
# validator-level C6 test and it drives .validate_fetwfe(), which has its own
# NESTED C6 variant and never calls this helper at all.
#
# Shape: twelve cells = three classes x four contracts, each cell a single
# mutation applied to a COPY of a live fit. One mutation per copy is mandatory
# here: .assert_contract() stops on the FIRST violated contract, so a fixture
# breaking several at once only ever demonstrates the earliest and the anchored
# assertions on the others pass vacuously (PR #452 hit exactly this trap).
#
# Nine of the twelve cells are redundant TODAY -- any single cell of a class's
# row kills that class's routing mutant. They stop being redundant the moment
# the shared helper is split per class, which is the refactor shape this issue
# exists because of.
#
# HOUSE RULE: every expect_error() anchor passes `fixed = TRUE`. The contract
# names contain parentheses and `*` (`C6 nrow(X_ints) == N * T`), which a
# default `regexp =` would read as metacharacters.

# ------------------------------------------------------------------------------
# Shared fixture: one simulated panel, three live fits. Recipe copied from
# tests/testthat/test-class-validators.R:21-36 (which is file-local there and so
# cannot be called from here; see the plan's Decision Log for why it was not
# promoted to a helper-*.R). suppressWarnings() alone, matching that source.
# ------------------------------------------------------------------------------
.c6r429_sim <- simulateData(
	genCoefs(G = 2, T = 3, d = 2, density = 0.5, eff_size = 1, seed = 17052026),
	N = 80,
	sig_eps_sq = 1,
	sig_eps_c_sq = 0.5,
	seed = 17052026
)

.c6r429_fits <- list(
	etwfe = suppressWarnings(etwfeWithSimulatedData(.c6r429_sim)),
	betwfe = suppressWarnings(betwfeWithSimulatedData(.c6r429_sim, q = 0.5)),
	twfeCovs = suppressWarnings(twfeCovsWithSimulatedData(.c6r429_sim))
)

.c6r429_validators <- list(
	etwfe = fetwfe:::.validate_etwfe,
	betwfe = fetwfe:::.validate_betwfe,
	twfeCovs = fetwfe:::.validate_twfeCovs
)

test_that("etwfe/betwfe/twfeCovs validators route through .check_c6_dims_toplevel() (#429 item 3)", {
	for (cls in names(.c6r429_fits)) {
		fit <- .c6r429_fits[[cls]]
		validate <- .c6r429_validators[[cls]]

		# Baseline. Without this the four errors below could be attributable to
		# a fixture that was already off-contract rather than to the mutation.
		expect_silent(validate(fit))

		# --- C6, three dimension contracts -----------------------------------
		cells <- list(
			beta_hat = list(
				mutate = function(o) {
					o$beta_hat <- o$beta_hat[-1]
					o
				},
				contract = "C6 length(beta_hat) == p"
			),
			y = list(
				mutate = function(o) {
					o$y <- o$y[-1]
					o
				},
				contract = "C6 length(y) == N * T"
			),
			X_ints = list(
				mutate = function(o) {
					o$X_ints <- o$X_ints[-1, , drop = FALSE]
					o
				},
				contract = "C6 nrow(X_ints) == N * T"
			)
		)

		for (cell in names(cells)) {
			expect_error(
				validate(cells[[cell]]$mutate(fit)),
				paste0(
					".validate_",
					cls,
					"(): contract violation `",
					cells[[cell]]$contract,
					"`"
				),
				fixed = TRUE,
				info = paste(cls, cell)
			)
		}

		# --- C8, calc_ses ----------------------------------------------------
		# C8 cannot be reached by the naive mutation, and the reason is
		# structural: isTRUE(x) is TRUE only for a length-1 logical with no
		# attributes -- exactly what C8 asserts -- so any object violating C8
		# makes isTRUE(calc_ses) FALSE, which sends .check_se_consistency() (C1,
		# which runs first) down its "standard errors must be NA" branch.
		# Reaching C8 therefore requires NA-ing the SE-bearing slots first,
		# which is the genuine calc_ses = FALSE shape. (`att_p_value` is the
		# slot; `p_value` is a COLUMN of catt_df.)
		na_only <- fit
		na_only$att_se <- NA_real_
		na_only$catt_ses <- rep(NA_real_, length(na_only$catt_ses))
		na_only$att_p_value <- NA_real_

		# Control: the NA-ing ALONE is legal. Without this, a future change that
		# made the NA-ing itself illegal would turn the C8 cell green for the
		# wrong reason.
		expect_silent(validate(na_only))

		mutated <- na_only
		mutated$calc_ses <- "TRUE"
		expect_error(
			validate(mutated),
			paste0(
				".validate_",
				cls,
				"(): contract violation `C8 calc_ses is length-1 logical`"
			),
			fixed = TRUE,
			info = paste(cls, "calc_ses")
		)
	}
})
