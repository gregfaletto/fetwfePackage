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
		R = 3,
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
		)
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
