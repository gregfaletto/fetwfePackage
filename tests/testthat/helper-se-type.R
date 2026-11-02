# ------------------------------------------------------------------------------
# Helpers for se_type tests: small F1-conforming and AR(1)-corrupted panels.
# simulateData() produces F1-conforming shocks only; for the serially-correlated
# fixture we add an AR(1) shock per unit on top of the simulated response.
# ------------------------------------------------------------------------------
make_se_type_panel <- function(seed = 2026, N = 120) {
	set.seed(seed)
	sim_coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	# coefs built without a seed: continue the ambient RNG (seed = NA, #250).
	simulateData(
		sim_coefs,
		N = N,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = NA
	)
}

make_ar1_panel <- function(seed = 7, N = 150, rho = 0.85, sd_e = 1) {
	set.seed(seed)
	sim_coefs <- genCoefs(G = 3, T = 6, d = 2, density = 0.5, eff_size = 2)
	# coefs built without a seed: continue the ambient RNG (seed = NA, #250).
	sim <- simulateData(
		sim_coefs,
		N = N,
		sig_eps_sq = 0.5,
		sig_eps_c_sq = 0.5,
		seed = NA
	)
	pdata <- sim$pdata
	pdata <- pdata[order(pdata$unit, pdata$time), ]
	units <- unique(pdata$unit)
	T_ <- length(unique(pdata$time))
	set.seed(seed + 100)
	for (u in units) {
		idx <- which(pdata$unit == u)
		ar <- numeric(T_)
		ar[1] <- rnorm(1, sd = sd_e / sqrt(1 - rho^2))
		for (t in 2:T_) {
			ar[t] <- rho * ar[t - 1] + rnorm(1, sd = sd_e)
		}
		pdata[idx, sim$response] <- pdata[idx, sim$response] + ar
	}
	sim$pdata <- pdata
	sim
}

expect_se_type_omitted_matches_default <- function(
	fit_with_sim,
	sim = make_se_type_panel()
) {
	res_omit <- fit_with_sim(sim)
	res_def <- fit_with_sim(sim, se_type = "default")
	for (s in c("att_se", "catt_ses", "att_hat", "catt_hats")) {
		expect_false(is.null(res_omit[[s]]))
		expect_false(is.null(res_def[[s]]))
	}
	expect_equal(res_omit$att_se, res_def$att_se)
	expect_equal(res_omit$catt_ses, res_def$catt_ses)
	expect_equal(res_omit$att_hat, res_def$att_hat)
	expect_equal(res_omit$catt_hats, res_def$catt_hats)
}
