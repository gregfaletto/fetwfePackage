# Doc-vs-live parity test for the public estimator @return blocks plus an
# @examples-presence check for every exported function. Together these guard
# the package against the documentation-drift class of bug that accumulated
# silently from v1.5.5 through v1.9.0 and required three separate fix-PRs
# (v1.9.5 soft bugs, v1.9.7 vestigial params, v1.9.8 docs sweep) to
# unwind. See issue #70 for the rationale.

# ------------------------------------------------------------------------------
# Dual-mode Rd lookup: prefer source-tree `man/` (always fresh in dev),
# fall back to the installed package's compiled Rd (correct under R CMD
# check, which installs from current source before running tests).
# ------------------------------------------------------------------------------

.get_rd_db <- function() {
	here <- normalizePath(getwd(), mustWork = FALSE)
	for (i in 1:6) {
		if (
			file.exists(file.path(here, "DESCRIPTION")) &&
				dir.exists(file.path(here, "man"))
		) {
			return(tools::Rd_db(dir = here))
		}
		parent <- dirname(here)
		if (parent == here) {
			break
		}
		here <- parent
	}
	tools::Rd_db(package = "fetwfe")
}

# ------------------------------------------------------------------------------
# Helper: walk a parsed Rd node tree and extract every \item{KEY} first-arg
# encountered, expanding comma-separated forms like
# \item{time_var, unit_var, treatment}{...} into three entries.
# ------------------------------------------------------------------------------

.extract_value_items <- function(rd) {
	value_idx <- which(sapply(
		rd,
		function(n) identical(attr(n, "Rd_tag"), "\\value")
	))
	if (length(value_idx) == 0L) {
		return(character(0))
	}
	value_node <- rd[[value_idx[1L]]]
	result <- character(0)
	walker <- function(node) {
		if (is.list(node)) {
			tag <- attr(node, "Rd_tag")
			if (identical(tag, "\\item") && length(node) >= 1L) {
				key_block <- node[[1L]]
				key <- paste(
					unlist(lapply(key_block, as.character)),
					collapse = ""
				)
				result <<- c(result, key)
			}
			for (child in node) {
				walker(child)
			}
		}
	}
	walker(value_node)
	unique(trimws(unlist(strsplit(result, ","))))
}

# Direct \item children of the \value{} node -- top-level @return slots only,
# NOT the nested $internal sub-slots. A slot documented twice at the top level
# is a bug (the #206 ci_type duplicate); the intentional top-level/$internal
# dual-documentation of X_ints / y / X_final / y_final / calc_ses (#154, #180)
# is not, so this deliberately does not recurse the way .extract_value_items().
.extract_top_level_value_items <- function(rd) {
	value_idx <- which(sapply(
		rd,
		function(n) identical(attr(n, "Rd_tag"), "\\value")
	))
	if (length(value_idx) == 0L) {
		return(character(0))
	}
	keys <- character(0)
	for (child in rd[[value_idx[1L]]]) {
		if (
			is.list(child) &&
				identical(attr(child, "Rd_tag"), "\\item") &&
				length(child) >= 1L
		) {
			keys <- c(
				keys,
				paste(unlist(lapply(child[[1L]], as.character)), collapse = "")
			)
		}
	}
	trimws(unlist(strsplit(keys, ",")))
}

# ------------------------------------------------------------------------------
# Test 1: @return slot parity for the 4 public estimators and their 3 wrappers.
#
# For each, fit on a small simulated panel, dump names(res) (flatten with
# names(res$internal) for fetwfe's nested sublist), and assert setequal()
# against the documented slot keys from the parsed man/<fn>.Rd.
#
# Mutation test (executor checklist, not in CI): delete an \item{} entry
# from one of these docstrings, re-document, re-run this test; it must
# fail with the missing slot name in the info= message.
# ------------------------------------------------------------------------------

test_that("public estimator @return blocks match live names()", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 0,
		density = 0.5,
		eff_size = 1,
		seed = 20260517
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	cases <- list(
		list("fetwfe", function() fetwfeWithSimulatedData(sim, q = 0.5)),
		list(
			"fetwfeWithSimulatedData",
			function() fetwfeWithSimulatedData(sim, q = 0.5)
		),
		list("etwfe", function() etwfeWithSimulatedData(sim)),
		list("etwfeWithSimulatedData", function() etwfeWithSimulatedData(sim)),
		list("betwfe", function() betwfeWithSimulatedData(sim, q = 0.5)),
		list(
			"betwfeWithSimulatedData",
			function() betwfeWithSimulatedData(sim, q = 0.5)
		),
		list("twfeCovs", function() twfeCovsWithSimulatedData(sim)),
		list(
			"twfeCovsWithSimulatedData",
			function() twfeCovsWithSimulatedData(sim)
		),
		list("genCoefs", function() coefs),
		list("getTes", function() getTes(coefs)),
		list("simulateData", function() sim)
	)

	db <- .get_rd_db()
	for (case in cases) {
		fn_name <- case[[1]]
		res <- suppressWarnings(case[[2]]())
		live <- names(res)
		if ("internal" %in% live) {
			live <- c(live, names(res$internal))
		}
		rd_key <- paste0(fn_name, ".Rd")
		expect_true(
			rd_key %in% names(db),
			info = paste0("Missing Rd entry in db: ", rd_key)
		)
		doc <- .extract_value_items(db[[rd_key]])
		missing <- setdiff(live, doc)
		extra <- setdiff(doc, live)
		expect_true(
			length(missing) == 0,
			info = paste0(
				fn_name,
				": ",
				length(missing),
				" slot(s) in live names() but missing from @return: ",
				paste(missing, collapse = ", ")
			)
		)
		expect_true(
			length(extra) == 0,
			info = paste0(
				fn_name,
				": ",
				length(extra),
				" slot(s) documented but not in live names(): ",
				paste(extra, collapse = ", ")
			)
		)
	}
})

# ------------------------------------------------------------------------------
# Test 1b (#206): no @return \item is documented twice at the top level.
#
# .extract_value_items()'s unique() (and the setequal() comparisons above) is
# duplicate-insensitive, so a doubled top-level \item{} collapsed to one entry
# and shipped to the rendered manual undetected (the ci_type duplicate in
# ?fetwfe / ?fetwfeWithSimulatedData that #206 fixed). Check only the top-level
# slots, excluding the intentional $internal dual-doc (#154, #180).
# ------------------------------------------------------------------------------

test_that("public estimator @return blocks contain no duplicate \\item names (#206)", {
	db <- .get_rd_db()
	rd_keys <- c(
		"fetwfe.Rd",
		"fetwfeWithSimulatedData.Rd",
		"etwfe.Rd",
		"etwfeWithSimulatedData.Rd",
		"betwfe.Rd",
		"betwfeWithSimulatedData.Rd",
		"twfeCovs.Rd",
		"twfeCovsWithSimulatedData.Rd"
	)
	for (rd_key in rd_keys) {
		raw <- .extract_top_level_value_items(db[[rd_key]])
		dups <- unique(raw[duplicated(raw)])
		expect_true(
			length(dups) == 0L,
			info = paste0(
				rd_key,
				" documents these @return \\item(s) more than once: ",
				paste(dups, collapse = ", ")
			)
		)
	}
})

# ------------------------------------------------------------------------------
# Test 2: every exported function has an @examples block.
#
# Uses getNamespaceExports() for export enumeration (works in both
# devtools::test and R CMD check modes) and the parsed Rd db's
# \examples{} block for the presence check. The class Rd files
# (fetwfe-class.Rd etc.) are not aliased to any export, so they are
# correctly excluded.
# ------------------------------------------------------------------------------

test_that("every exported function has an @examples block", {
	exports <- sort(getNamespaceExports("fetwfe"))
	db <- .get_rd_db()
	# Map each exported function name to its Rd file (the convention is
	# `man/<fn>.Rd`).
	for (fn in exports) {
		rd_key <- paste0(fn, ".Rd")
		expect_true(
			rd_key %in% names(db),
			info = paste0("Missing Rd entry for exported function: ", fn)
		)
		if (rd_key %in% names(db)) {
			rd <- db[[rd_key]]
			has_examples <- any(sapply(
				rd,
				function(n) identical(attr(n, "Rd_tag"), "\\examples")
			))
			expect_true(
				has_examples,
				info = paste0(
					fn,
					": exported but man/",
					fn,
					".Rd has no \\examples block. ",
					"Public functions require examples (per AGENTS.md)."
				)
			)
		}
	}
})

# ------------------------------------------------------------------------------
# Test 3: cross-class slot parity (issue #93 Gap 1).
#
# For each slot name in the union of the four .EXPECTED_SLOTS_<CLASS> vectors
# (treating fetwfe's $internal sublist as flattened — its child slots are
# top-level in the other classes), assert the slot appears in EITHER all four
# classes (the default) OR in exactly the classes listed in a `DIVERGENT_SLOTS`
# table entry. The table records intentional divergences with a rationale.
#
# What this catches: a PR that adds a new slot to one estimator's output
# without adding it to the others (or that intentionally diverges without
# documenting the rationale). Test 1 above already catches per-class
# code-vs-doc drift; this Test 3 catches CROSS-class drift.
# ------------------------------------------------------------------------------

test_that("cross-class slot inventory matches the documented divergence table", {
	# Divergent slots: ones that DON'T appear in all 4 classes by design.
	# Each entry's `classes` is the exact set the slot should appear in;
	# `rationale` documents why.
	DIVERGENT_SLOTS <- list(
		# `alpha` was divergent until #204; twfeCovs now carries it too (it
		# threads the user's alpha into the default simultaneous catt_df band),
		# so `alpha` is a universal slot (no DIVERGENT_SLOTS entry).
		"att_selected" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only the selection-based (bridge-penalty) estimators perform selection."
		),
		"lambda.max" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda regularization path."
		),
		"lambda.max_model_size" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda regularization path."
		),
		"lambda.min" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda regularization path."
		),
		"lambda.min_model_size" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda regularization path."
		),
		"lambda_star" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda regularization path."
		),
		"lambda_star_model_size" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda regularization path."
		),
		"theta_hat" = list(
			classes = c("fetwfe"),
			rationale = "Only fetwfe applies the bridge fusion transform; theta_hat is the transformed coefficient vector."
		),
		# v1.13.0 (#164): the CV-tuned lambda-selection knob lives only on
		# the bridge-penalty estimators.
		"lambda_selection" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda penalty to tune."
		),
		"cv_folds" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda penalty to tune."
		),
		"cv_seed" = list(
			classes = c("fetwfe", "betwfe"),
			rationale = "Only bridge-penalty estimators have a lambda penalty to tune."
		)
	)

	# Effective slot set per class: the `$internal` sublist is flattened
	# for ALL four classes so its child slots count as top-level for
	# cross-class comparison. (User-access paths are different —
	# `fit$X_ints` vs `fit$internal$X_ints` — but the value, semantic,
	# and presence are equivalent across classes.) Post-#144 the
	# OLS-family also carries `$internal` (in addition to duplicating
	# the same slots at top level for backward compat); flattening
	# all four classes here keeps the comparison symmetric.
	effective_slots <- list(
		fetwfe = c(
			setdiff(fetwfe:::.EXPECTED_SLOTS_FETWFE, "internal"),
			fetwfe:::.EXPECTED_INTERNAL_SLOTS_FETWFE
		),
		etwfe = c(
			setdiff(fetwfe:::.EXPECTED_SLOTS_ETWFE, "internal"),
			fetwfe:::.EXPECTED_INTERNAL_SLOTS_ETWFE
		),
		betwfe = c(
			setdiff(fetwfe:::.EXPECTED_SLOTS_BETWFE, "internal"),
			fetwfe:::.EXPECTED_INTERNAL_SLOTS_BETWFE
		),
		twfeCovs = c(
			setdiff(fetwfe:::.EXPECTED_SLOTS_TWFECOVS, "internal"),
			fetwfe:::.EXPECTED_INTERNAL_SLOTS_TWFECOVS
		)
	)
	all_classes <- names(effective_slots)
	universe <- unique(unlist(effective_slots))

	for (slot in universe) {
		observed <- all_classes[sapply(
			effective_slots,
			function(s) slot %in% s
		)]
		if (slot %in% names(DIVERGENT_SLOTS)) {
			expected <- DIVERGENT_SLOTS[[slot]]$classes
			rationale <- DIVERGENT_SLOTS[[slot]]$rationale
		} else {
			expected <- all_classes
			rationale <- "(no DIVERGENT_SLOTS entry — slot expected in all 4 classes)"
		}
		expect_true(
			setequal(observed, expected),
			info = paste0(
				"Cross-class slot parity mismatch for `",
				slot,
				"`:\n",
				"  observed in classes: ",
				paste(observed, collapse = ", "),
				"\n",
				"  expected in classes: ",
				paste(expected, collapse = ", "),
				"\n",
				"  rationale: ",
				rationale,
				"\n",
				"  Fix: either add `",
				slot,
				"` to the missing class(es)' .EXPECTED_SLOTS_<CLASS> vector, ",
				"or add/update a DIVERGENT_SLOTS entry in this test with a rationale."
			)
		)
	}
})

# ------------------------------------------------------------------------------
# Test 4 (#84 item 4): eventStudy() output schema lock-in.
#
# Parallel to Test 1 but for the event-study aggregation surface. The
# eventStudy() @return uses \describe{\item{KEY}{...}} structure so we
# reuse the existing .extract_value_items() parser. Covers all three
# estimator classes that eventStudy() accepts (fetwfe, etwfe, betwfe;
# twfeCovs is rejected with a stop() per R/event_study.R:67-69).
# ------------------------------------------------------------------------------

test_that("eventStudy() @return matches live names() across estimator classes (#84 item 4)", {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 0,
		density = 0.5,
		eff_size = 1,
		seed = 20260517
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	fits <- list(
		fetwfe = suppressWarnings(fetwfeWithSimulatedData(sim, q = 0.5)),
		etwfe = suppressWarnings(etwfeWithSimulatedData(sim)),
		betwfe = suppressWarnings(betwfeWithSimulatedData(sim, q = 0.5))
	)

	db <- .get_rd_db()
	expect_true(
		"eventStudy.Rd" %in% names(db),
		info = "eventStudy.Rd missing from Rd db"
	)
	doc_cols <- .extract_value_items(db[["eventStudy.Rd"]])

	for (cls in names(fits)) {
		es <- suppressWarnings(eventStudy(fits[[cls]]))
		live <- names(es)
		missing <- setdiff(live, doc_cols)
		extra <- setdiff(doc_cols, live)
		expect_true(
			length(missing) == 0,
			info = paste0(
				cls,
				": live names() not in eventStudy.Rd @return: ",
				paste(missing, collapse = ", ")
			)
		)
		expect_true(
			length(extra) == 0,
			info = paste0(
				cls,
				": eventStudy.Rd @return slots not in live names(): ",
				paste(extra, collapse = ", ")
			)
		)
	}
})

# ------------------------------------------------------------------------------
# Test 5 (#84 item 4): tidy.eventStudy() output schema lock-in.
#
# The tidy method's @return is prose (not \describe{}-structured) and the
# schema is conditional on conf.int. Hard-coded expected column sets;
# future drift here forces a conscious update of either the test or the
# tidy implementation.
# ------------------------------------------------------------------------------

test_that("tidy.eventStudy() schema is locked across conf.int branches (#84 item 4)", {
	skip_if_not_installed("broom")
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 0,
		density = 0.5,
		eff_size = 1,
		seed = 20260517
	)
	sim <- simulateData(
		coefs,
		N = 120,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)

	fits <- list(
		fetwfe = suppressWarnings(fetwfeWithSimulatedData(sim, q = 0.5)),
		etwfe = suppressWarnings(etwfeWithSimulatedData(sim)),
		betwfe = suppressWarnings(betwfeWithSimulatedData(sim, q = 0.5))
	)

	# Documented schema (broom convention) — hard-coded because the
	# @return prose isn't \describe{}-structured. Source: roxygen at
	# R/broom_methods.R:423-425.
	expected_with_ci <- c(
		"term",
		"event_time",
		"n_cohorts",
		"estimate",
		"std.error",
		"statistic",
		"p.value",
		"conf.low",
		"conf.high"
	)
	expected_no_ci <- setdiff(expected_with_ci, c("conf.low", "conf.high"))

	for (cls in names(fits)) {
		es <- suppressWarnings(eventStudy(fits[[cls]]))

		td_ci <- broom::tidy(es)
		expect_identical(
			names(td_ci),
			expected_with_ci,
			info = paste0(
				cls,
				": tidy(eventStudy()) default columns drift"
			)
		)

		td_no_ci <- broom::tidy(es, conf.int = FALSE)
		expect_identical(
			names(td_no_ci),
			expected_no_ci,
			info = paste0(
				cls,
				": tidy(eventStudy(), conf.int = FALSE) columns drift"
			)
		)
	}
})
