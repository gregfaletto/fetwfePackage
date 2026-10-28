# Checks #495's withholding of debiasedATT() and debiasedATTWithSimulatedData()
# until #527 re-exports them: the exports, the S3 registrations of their
# result's methods, and the text of the REML and calc_ses stops.

# The method registered for a generic, read from the S3 registry in the
# generic's namespace rather than looked up by name.
.registry_lookup <- function(generic_ns, method) {
	tbl <- get(".__S3MethodsTable__.", envir = asNamespace(generic_ns))
	get0(method, envir = tbl, inherits = FALSE)
}

# One simulated design; `N` sets the regime.
.sim_panel_495 <- function(N) {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 2,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	simulateData(coefs, N = N, sig_eps_sq = 1, sig_eps_c_sq = 0.5, seed = 1)
}

test_that("debiasedATT() and its wrapper are withheld but still present (#495)", {
	ex <- getNamespaceExports("fetwfe")
	expect_true("fetwfe" %in% ex)
	expect_false("debiasedATT" %in% ex)
	expect_false("debiasedATTWithSimulatedData" %in% ex)
	ns <- asNamespace("fetwfe")
	expect_true(is.function(get0("debiasedATT", envir = ns, inherits = FALSE)))
	expect_true(is.function(get0(
		"debiasedATTWithSimulatedData",
		envir = ns,
		inherits = FALSE
	)))
	expect_false(is.null(.registry_lookup("base", "print.fetwfe")))
	expect_false(is.null(.registry_lookup("base", "print.debiased_att")))
	expect_false(is.null(.registry_lookup("generics", "tidy.debiased_att")))
})

test_that("the REML and a fixed-p calc_ses stop name neither debiasedATT() nor #312, and the latter advises a q < 1 re-fit, not the bootstrap (#495)", {
	s12 <- .sim_panel_495(12)
	reml_msg <- msg_of(tryCatch(
		fetwfe(
			pdata = s12$pdata,
			time_var = s12$time_var,
			unit_var = s12$unit_var,
			treatment = s12$treatment,
			response = s12$response,
			covs = s12$covs,
			verbose = FALSE
		),
		error = identity
	))
	expect_match(
		reml_msg,
		"Not enough units available to estimate the noise variance",
		fixed = TRUE
	)
	expect_no_match(reml_msg, "debiasedATT", fixed = TRUE)
	expect_no_match(reml_msg, "#312", fixed = TRUE)

	fixed_p <- fetwfeWithSimulatedData(.sim_panel_495(30), q = 2)
	expect_false(isTRUE(fixed_p$internal$calc_ses))
	expect_lt(ncol(fixed_p$internal$X_final), nrow(fixed_p$internal$X_final))
	calc_ses_msg <- msg_of(tryCatch(
		simultaneousCIs(fixed_p, family = "cohort"),
		error = identity
	))
	expect_match(calc_ses_msg, "calc_ses = FALSE", fixed = TRUE)
	expect_no_match(calc_ses_msg, "debiasedATT", fixed = TRUE)
	expect_no_match(calc_ses_msg, "#312", fixed = TRUE)
	expect_match(calc_ses_msg, "re-fit with q < 1", fixed = TRUE)
	expect_no_match(
		calc_ses_msg,
		fetwfe:::.highdim_band_remedy(TRUE),
		fixed = TRUE
	)
})

test_that("a high-dimensional fit's calc_ses stop gives its estimator's band remedy and names neither debiasedATT() nor #312 (#495)", {
	s9 <- .sim_panel_495(9)

	hd_fetwfe <- fetwfeWithSimulatedData(s9, q = 2)
	expect_false(isTRUE(hd_fetwfe$internal$calc_ses))
	expect_gte(
		ncol(hd_fetwfe$internal$X_final),
		nrow(hd_fetwfe$internal$X_final)
	)
	fetwfe_msg <- msg_of(tryCatch(
		simultaneousCIs(hd_fetwfe, family = "cohort"),
		error = identity
	))
	expect_match(fetwfe_msg, "calc_ses = FALSE", fixed = TRUE)
	expect_no_match(fetwfe_msg, "debiasedATT", fixed = TRUE)
	expect_no_match(fetwfe_msg, "#312", fixed = TRUE)
	expect_match(fetwfe_msg, fetwfe:::.highdim_band_remedy(TRUE), fixed = TRUE)
	expect_no_match(
		fetwfe_msg,
		fetwfe:::.highdim_band_remedy(FALSE),
		fixed = TRUE
	)

	hd_betwfe <- betwfeWithSimulatedData(s9, q = 1)
	expect_false(isTRUE(hd_betwfe$calc_ses))
	expect_gte(ncol(hd_betwfe$X_final), nrow(hd_betwfe$X_final))
	betwfe_msg <- msg_of(tryCatch(
		simultaneousCIs(hd_betwfe, family = "cohort"),
		error = identity
	))
	expect_match(betwfe_msg, "calc_ses = FALSE", fixed = TRUE)
	expect_no_match(betwfe_msg, "debiasedATT", fixed = TRUE)
	expect_no_match(betwfe_msg, "#312", fixed = TRUE)
	expect_match(betwfe_msg, fetwfe:::.highdim_band_remedy(FALSE), fixed = TRUE)
	expect_no_match(
		betwfe_msg,
		fetwfe:::.highdim_band_remedy(TRUE),
		fixed = TRUE
	)
})
