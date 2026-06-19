# Tests for #309: desparsify the high-dim (p >= NT) event_study propensity channel
# F_pi. The band CENTER is debiased (Theorem 6.6), but F_pi (the Sigma_2 cohort-
# weight variance) was built from the POST-SELECTION theta_sel, so a bridge-zeroed
# cohort contributed 0 to Sigma_2 while the debiased center counted its non-zero
# value -- a center/variance inconsistency that undercovers in the heterogeneous-
# cohort regime. The fix feeds the DEBIASED per-cohort-time effects (tau_db) into
# F_pi via the SAME q=1 nuisance + nodewise machinery the center uses.
#
# A[g,k] is a delta-method Jacobian contraction whose coefficients depend only on
# cohort_probs; only the per-cell tau values are post-selection. So the fix
# desparsifies the num_treats cells (tau_db) and re-pools them through the SAME
# Jacobian coefficients (j_cells, built on the identity per-cell map):
#   A_db[g, k] = (j_cells[[k]] %*% tau_db)[g].
#
# MUTATION ANCHOR (run manually to confirm the integration bites): in
# R/simultaneous_bootstrap.R's F_pi_mat dispatch, revert `A = A_db` to `A = NULL`
# (post-selection); the high-dim event_study band ses change and the end-to-end
# test below ("... matches the desparsified reconstruction") FAILS. Revert.

# ---- shared high-dim event_study fixture (G=3, T=5, d=22 => p=390 > NT=200) ----
.make_hd309_fit <- function() {
	coefs <- genCoefs(
		G = 3,
		T = 5,
		d = 22,
		density = 0.5,
		eff_size = 2,
		seed = 1
	)
	dat <- simulateData(
		coefs,
		N = 40,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5,
		seed = 1
	)
	dat$indep_counts <- NA
	fetwfe(
		pdata = dat$pdata,
		time_var = dat$time_var,
		unit_var = dat$unit_var,
		treatment = dat$treatment,
		response = dat$response,
		covs = dat$covs,
		q = 0.5,
		verbose = FALSE,
		sig_eps_sq = 1,
		sig_eps_c_sq = 0.5
	)
}

hd309 <- .make_hd309_fit()

# Reconstruct the dispatch quantities the desparsified F_pi is built from.
.hd309_parts <- function(hd) {
	X <- hd$internal$X_final
	y <- as.numeric(hd$internal$y_final)
	n <- nrow(X)
	p <- ncol(X)
	Tt <- hd$T
	G <- hd$G
	ti <- hd$treat_inds
	nt <- length(ti)
	fi <- getFirstInds(G = G, T = Tt)
	A <- genFullInvFusionTransformMat(
		first_inds = fi,
		T = Tt,
		G = G,
		d = hd$d,
		num_treats = nt,
		fusion_structure = hd$fusion_structure,
		d_inv_treat = hd$internal$d_inv_treat
	)
	theta_q1 <- fetwfe:::.fit_q1_nuisance(X, y, n / Tt, Tt)
	resid <- as.numeric(y - theta_q1[1] - X %*% theta_q1[-1])
	E <- matrix(0, p, nt)
	E[cbind(ti, seq_len(nt))] <- 1
	cell_targets <- crossprod(A, E) # p x num_treats
	offs <- fetwfe:::.resolve_event_study_offsets_and_first_inds(
		hd,
		G = G,
		T = Tt
	)
	j_cells <- fetwfe:::.build_j_list_for_family(
		family = "event_study",
		K = Tt - 1L,
		G = G,
		T = Tt,
		num_treats = nt,
		cohort_offsets_int = offs$cohort_offsets_int,
		first_inds = offs$first_inds,
		cohort_probs_overall = hd$cohort_probs_overall,
		d_inv_treat_sel = diag(nt)
	)
	list(
		X = X,
		n = n,
		p = p,
		T = Tt,
		G = G,
		ti = ti,
		nt = nt,
		theta_q1 = theta_q1,
		resid = resid,
		cell_targets = cell_targets,
		j_cells = j_cells,
		cohort_probs_overall = hd$cohort_probs_overall,
		fi = fi,
		cohort_offsets_int = offs$cohort_offsets_int,
		first_inds = offs$first_inds
	)
}

test_that("the #309 fixture is genuinely high-dim with valid SEs", {
	expect_gte(ncol(hd309$internal$X_final), nrow(hd309$internal$X_final))
	expect_true(isTRUE(hd309$internal$calc_ses))
})

test_that(".build_debiased_treat_cells_highdim matches an independent q=1 + nodewise reconstruction", {
	pr <- .hd309_parts(hd309)
	cells <- fetwfe:::.build_debiased_treat_cells_highdim(
		pr$X,
		pr$resid,
		pr$n / pr$T,
		pr$T,
		pr$theta_q1[-1],
		pr$cell_targets
	)
	# Independent reconstruction: per-cell q=1 plug-in + manual nodewise correction
	# (riesz_lasso + hand-rolled score), NOT the helper's internal loop.
	Sig <- crossprod(pr$X) / pr$n
	recon <- vapply(
		seq_len(pr$nt),
		function(j) {
			ath <- pr$cell_targets[, j]
			ln <- lambda_node_default(
				p = pr$p,
				N = pr$n / pr$T,
				c = 1.0,
				scale = max(abs(ath))
			)
			v <- riesz_lasso(Sig, ath, ln)
			sum(ath * pr$theta_q1[-1]) + mean((pr$X %*% v) * pr$resid)
		},
		numeric(1)
	)
	expect_equal(cells$tau_db, recon, tolerance = 1e-9)
	expect_length(cells$tau_db, pr$nt)
	expect_true(all(
		c("feasibility", "converged", "lambda_node") %in%
			names(cells$diagnostics)
	))
})

test_that(".build_propensity_if uses a supplied A; different A -> different Sigma_2", {
	# Pure assembler on synthetic G x K matrices (no fit needed).
	G <- 3L
	K <- 4L
	N <- 30L
	Tt <- 5L
	cpo <- c(0.3, 0.3, 0.2) # treated-cohort probs (never-treated mass implied)
	A1 <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1), nrow = G, ncol = K)
	A2 <- A1
	A2[1, 1] <- 5
	Fpi1 <- fetwfe:::.build_propensity_if(
		cohort_probs_overall = cpo,
		G = G,
		N = N,
		T = Tt,
		A = A1
	)
	Fpi2 <- fetwfe:::.build_propensity_if(
		cohort_probs_overall = cpo,
		G = G,
		N = N,
		T = Tt,
		A = A2
	)
	expect_equal(dim(Fpi1), c(N, K))
	s2_1 <- colSums(Fpi1^2) / (N * Tt)^2
	s2_2 <- colSums(Fpi2^2) / (N * Tt)^2
	expect_false(isTRUE(all.equal(s2_1, s2_2)))
})

test_that("zeroing a cohort's cells (the bug) changes Sigma_2 -- desparsification is the fix (#309)", {
	# The post-selection bug zeroes a cohort's cells -> 0 in F_pi. The fix uses the
	# non-zero DEBIASED tau_db. Show that the cohort's debiased effect is non-zero
	# AND that zeroing it (the bug's effect) materially changes Sigma_2 -- so using
	# the debiased values, not the zeroed post-selection ones, matters.
	pr <- .hd309_parts(hd309)
	cells <- fetwfe:::.build_debiased_treat_cells_highdim(
		pr$X,
		pr$resid,
		pr$n / pr$T,
		pr$T,
		pr$theta_q1[-1],
		pr$cell_targets
	)
	tau_db <- cells$tau_db
	G <- pr$G
	K <- pr$T - 1L
	A_db <- matrix(
		vapply(
			seq_len(K),
			function(k) as.numeric(pr$j_cells[[k]] %*% tau_db),
			numeric(G)
		),
		nrow = G,
		ncol = K
	)
	# Cohort 1's treatment cells; the debiased effect there is non-zero (so the
	# post-selection bug, which can zero it, would drop a real contribution).
	g1_cells <- pr$fi[1]:(pr$fi[2] - 1L)
	expect_gt(max(abs(tau_db[g1_cells])), 1e-6)
	tau_bug <- tau_db
	tau_bug[g1_cells] <- 0 # mimic post-selection zeroing of cohort 1
	A_bug <- matrix(
		vapply(
			seq_len(K),
			function(k) as.numeric(pr$j_cells[[k]] %*% tau_bug),
			numeric(G)
		),
		nrow = G,
		ncol = K
	)
	expect_false(isTRUE(all.equal(A_db, A_bug)))
	Fpi_db <- fetwfe:::.build_propensity_if(
		cohort_probs_overall = pr$cohort_probs_overall,
		G = G,
		N = pr$n / pr$T,
		T = pr$T,
		A = A_db
	)
	Fpi_bug <- fetwfe:::.build_propensity_if(
		cohort_probs_overall = pr$cohort_probs_overall,
		G = G,
		N = pr$n / pr$T,
		T = pr$T,
		A = A_bug
	)
	expect_false(isTRUE(all.equal(colSums(Fpi_db^2), colSums(Fpi_bug^2))))
})

test_that("high-dim event_study band matches the desparsified reconstruction (end-to-end, #309)", {
	# The load-bearing end-to-end check: independently reconstruct the per-effect
	# ses from BOTH channels using the DEBIASED A_db, and confirm the shipped band
	# matches. (Reverting A = A_db -> A = NULL in the dispatch -- the mutation
	# anchor -- changes css_pi and breaks this match.)
	pr <- .hd309_parts(hd309)
	N <- pr$n / pr$T
	K <- pr$T - 1L
	# Regression channel css_reg: the event-time family targets are the pooled cell
	# targets, targets[,k] = sum_g w_{g,k} cell_targets[,g-cell].
	offs <- fetwfe:::.resolve_event_study_offsets_and_first_inds(
		hd309,
		G = pr$G,
		T = pr$T
	)
	event_times <- 0:(pr$T - 2L)
	W <- matrix(0, pr$nt, K) # num_treats x K pooling weights
	for (kk in seq_len(K)) {
		e <- event_times[kk]
		V_e <- which(offs$cohort_offsets_int <= pr$T - e)
		if (!length(V_e)) {
			next
		}
		wv <- pr$cohort_probs_overall[V_e] / sum(pr$cohort_probs_overall[V_e])
		for (j in seq_along(V_e)) {
			W[offs$first_inds[V_e[j]] + e, kk] <- wv[j]
		}
	}
	targets <- pr$cell_targets %*% W # p x K
	F_reg <- fetwfe:::.build_regression_if_highdim(
		pr$X,
		pr$resid,
		N,
		pr$T,
		targets
	)
	css_reg <- colSums(F_reg^2)
	# Propensity channel css_pi from the DEBIASED A_db.
	cells <- fetwfe:::.build_debiased_treat_cells_highdim(
		pr$X,
		pr$resid,
		N,
		pr$T,
		pr$theta_q1[-1],
		pr$cell_targets
	)
	A_db <- matrix(
		vapply(
			seq_len(K),
			function(k) as.numeric(pr$j_cells[[k]] %*% cells$tau_db),
			numeric(pr$G)
		),
		nrow = pr$G,
		ncol = K
	)
	F_pi <- fetwfe:::.build_propensity_if(
		cohort_probs_overall = pr$cohort_probs_overall,
		G = pr$G,
		N = N,
		T = pr$T,
		A = A_db
	)
	css_pi <- colSums(F_pi^2)
	cadjust <- N / (N - 1)
	se_recon <- sqrt((cadjust * css_reg + css_pi) / pr$n^2)
	sc <- simultaneousCIs(
		hd309,
		family = "event_study",
		method = "bootstrap",
		B = 200,
		seed = 1
	)
	se_band <- (sc$ci$simultaneous_ci_high - sc$ci$estimate) /
		sc$critical_value
	expect_equal(se_band, se_recon, tolerance = 1e-7)
})

test_that("high-dim event_study band is finite and deterministic with the desparsified F_pi", {
	sc <- simultaneousCIs(
		hd309,
		family = "event_study",
		method = "bootstrap",
		B = 300,
		seed = 1
	)
	sc2 <- simultaneousCIs(
		hd309,
		family = "event_study",
		method = "bootstrap",
		B = 300,
		seed = 1
	)
	expect_identical(sc$regime, "high-dimensional")
	expect_true(all(is.finite(sc$ci$simultaneous_ci_low)))
	expect_true(all(is.finite(sc$ci$simultaneous_ci_high)))
	expect_identical(sc$critical_value, sc2$critical_value)
	expect_identical(sc$ci$estimate, sc2$ci$estimate)
})

test_that("the experimental warning COUNTS the per-cell propensity directions, and they are exposed (#309)", {
	# Count-aware (not vacuous): at a tiny lambda_c the warning's total-direction
	# count must include the num_treats per-cell propensity directions, not just the
	# K event-time effects. Reverting the cell fold-in (the mutation anchor) drops
	# the total to K and fails this. Also pins that the per-cell diagnostics are
	# now exposed on the returned object (so the warning's "inspect ..." is real).
	w <- NULL
	sc <- withCallingHandlers(
		simultaneousCIs(
			hd309,
			family = "event_study",
			method = "bootstrap",
			B = 100,
			seed = 1,
			lambda_c = 0.02
		),
		warning = function(cond) {
			w <<- conditionMessage(cond)
			invokeRestart("muffleWarning")
		}
	)
	expect_false(is.null(w)) # the experimental warning fired
	total <- as.integer(sub(".*of ([0-9]+) high-dimensional.*", "\\1", w))
	n_eff <- length(sc$converged) # K event-time effects
	n_cells <- length(hd309$treat_inds) # num_treats per-cell directions
	expect_equal(total, n_eff + n_cells) # cells are folded into the count
	# The per-cell diagnostics are exposed and align with num_treats.
	expect_length(sc$propensity_feasibility, n_cells)
	expect_length(sc$propensity_converged, n_cells)
	expect_length(sc$propensity_lambda_node, n_cells)
})

test_that("the default-lambda_c high-dim event_study fit exposes feasible per-cell diagnostics", {
	sc <- simultaneousCIs(
		hd309,
		family = "event_study",
		method = "bootstrap",
		B = 100,
		seed = 1
	)
	expect_length(sc$propensity_converged, length(hd309$treat_inds))
	expect_true(all(sc$propensity_converged))
	expect_true(all(
		sc$propensity_feasibility <= sc$propensity_lambda_node * (1 + 1e-9)
	))
})

test_that("the identity-map j_cells pooling matches the closed-form delta-method formula (#309 anchor)", {
	# Independent (non-circular) anchor for A_db[g,k] = (j_cells[[k]] %*% tau_db)[g]:
	# rebuild it from the closed-form per_effect_masked coefficients
	#   A[g,k] = cons_g * tau_{g,e} - sum_{g'!=g, g' in V_e} cons_off(g') * tau_{g',e},
	#   cons_g = (S_V - pi_g)/S_V^2,  cons_off(g') = pi_g'/S_V^2,  S_V = sum_{V_e} pi,
	# WITHOUT calling .build_j_list_for_family. A wrong Jacobian pooling fails this.
	pr <- .hd309_parts(hd309)
	cells <- fetwfe:::.build_debiased_treat_cells_highdim(
		pr$X,
		pr$resid,
		pr$n / pr$T,
		pr$T,
		pr$theta_q1[-1],
		pr$cell_targets
	)
	tau_db <- cells$tau_db
	G <- pr$G
	K <- pr$T - 1L
	pi <- pr$cohort_probs_overall
	A_jcells <- matrix(
		vapply(
			seq_len(K),
			function(k) as.numeric(pr$j_cells[[k]] %*% tau_db),
			numeric(G)
		),
		nrow = G,
		ncol = K
	)
	A_closed <- matrix(0, G, K)
	for (k in seq_len(K)) {
		e <- k - 1L
		V_e <- which(pr$cohort_offsets_int <= pr$T - e)
		if (length(V_e) <= 1L) {
			next # |V_e| <= 1 -> zero Jacobian (single/no cohort at this event time)
		}
		S_V <- sum(pi[V_e])
		for (g in V_e) {
			cons_g <- (S_V - pi[g]) / S_V^2
			val <- cons_g * tau_db[pr$first_inds[g] + e]
			for (gp in setdiff(V_e, g)) {
				val <- val - (pi[gp] / S_V^2) * tau_db[pr$first_inds[gp] + e]
			}
			A_closed[g, k] <- val
		}
	}
	expect_equal(A_jcells, A_closed, tolerance = 1e-10)
})
