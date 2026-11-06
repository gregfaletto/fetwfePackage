test_that(".compute_p_values() matches known answers, including P-values far below the tolerance", {
	expect_identical(fetwfe:::.compute_p_values(0, 1), 1)
	expect_equal(fetwfe:::.compute_p_values(stats::qnorm(0.975), 1), 0.05)
	expect_equal(
		log(fetwfe:::.compute_p_values(5, 0.5)),
		log(1.5239706048321054e-23),
		tolerance = 1e-12
	)
	expect_equal(
		log(fetwfe:::.compute_p_values(60, 2)),
		log(9.8134278542963744e-198),
		tolerance = 1e-12
	)
	expect_equal(
		fetwfe:::.compute_p_values(-2, 1),
		fetwfe:::.compute_p_values(2, 1)
	)
	expect_identical(fetwfe:::.compute_p_values(1, 0), NA_real_)
	expect_identical(fetwfe:::.compute_p_values(1, NA_real_), NA_real_)
	expect_identical(fetwfe:::.compute_p_values(1, -1), NA_real_)
})
