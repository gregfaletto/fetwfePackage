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
