.fake_fit_541 <- function(changed = NULL, dropped = NULL) {
	function(sim, se_type = NULL) {
		res <- list(
			att_se = 0.5,
			catt_ses = c(0.25, 0.75),
			att_hat = 2,
			catt_hats = c(1, 3)
		)
		if (identical(se_type, "default") && !is.null(changed)) {
			res[[changed]] <- res[[changed]] * (1 + 1e-6)
		}
		res[dropped] <- NULL
		res
	}
}

test_that("expect_se_type_omitted_matches_default() fails when one slot differs under se_type = 'default' (#541)", {
	for (s in c("att_se", "catt_ses", "att_hat", "catt_hats")) {
		expect_error(
			expect_se_type_omitted_matches_default(
				.fake_fit_541(changed = s),
				sim = list()
			),
			class = "expectation_failure",
			info = s
		)
	}
})

test_that("expect_se_type_omitted_matches_default() fails when a slot is missing from both fits (#541)", {
	expect_error(
		expect_se_type_omitted_matches_default(
			.fake_fit_541(dropped = "catt_ses"),
			sim = list()
		),
		class = "expectation_failure"
	)
})

test_that("expect_se_type_omitted_matches_default() fits with cohort probabilities from the same data (#541)", {
	fake <- function(sim, se_type = "conservative") {
		same_data <- all(is.na(sim$indep_counts))
		conservative <- identical(se_type, "conservative")
		list(
			att_se = if (conservative && same_data) 0.7 else 0.5,
			catt_ses = c(0.25, 0.75),
			att_hat = 2,
			catt_hats = c(1, 3)
		)
	}
	expect_error(
		expect_se_type_omitted_matches_default(
			fake,
			sim = list(indep_counts = c(10L, 20L))
		),
		class = "expectation_failure"
	)
})

test_that("estimator_fit_with_sim() forwards the simulation's fields and its extra arguments (#541)", {
	seen <- NULL
	fake <- function(...) {
		seen <<- list(...)
		list()
	}
	sim <- list(
		pdata = "P",
		time_var = "t",
		unit_var = "u",
		treatment = "w",
		response = "y",
		covs = "c",
		indep_counts = NA,
		sig_eps_sq = 1,
		sig_eps_c_sq = 2
	)
	estimator_fit_with_sim(fake)(sim)
	expect_identical(seen, sim)
	estimator_fit_with_sim(fake)(sim, se_type = "default")
	expect_identical(seen, c(sim, list(se_type = "default")))
})
