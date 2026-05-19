# Shared panel-data test fixtures.
#
# testthat sources every `helper-*.R` in `tests/testthat/` once at the
# start of `devtools::test()` / `testthat::test_check()`, BEFORE any
# `test-*.R` runs. So the three helpers defined here are in scope for
# every test file in this directory without any `source()` call.
#
# This file consolidates what was previously 14 inline copies across 6
# test files (issue #91): 6 copies of `generate_panel_data`, 4 copies of
# `generate_bad_panel_data`, and 4 copies of `generate_minimal_panel_data`
# (with one logical-difference line — cov1 inclusion — now parameterized
# via the `include_cov1` argument).
#
# Behavior is byte-identical to the prior inline definitions:
# `generate_panel_data` and `generate_bad_panel_data` are direct copies
# of the canonical commented variants (the ones in `test-fetwfe.R` /
# `test-etwfe.R` / `test-betwfe.R` / `test-twfeCovs.R` pre-#91).
# `generate_minimal_panel_data(include_cov1 = TRUE)` matches the
# pre-#91 fetwfe/betwfe versions; `include_cov1 = FALSE` matches the
# pre-#91 etwfe/twfeCovs versions. The number and order of `rnorm()` /
# `runif()` / `sample()` calls is preserved in both branches so the
# fixed-seed output is unchanged.

# ------------------------------------------------------------------------------
# Helper function: generate a synthetic balanced panel data set.
# We now use N=30 units and T=10 periods so that the degrees‐of‐freedom condition
# (N * (T - 1) - p > 0) holds.
# For each unit, we randomly decide whether the unit is never treated (first_treat = Inf)
# or, if treated, assign a first treatment time (ensuring that no unit is treated in time 1).
# ------------------------------------------------------------------------------
generate_panel_data <- function(N = 30, T = 10, R = 9, seed = 123) {
	set.seed(seed)
	stopifnot(R <= T - 1)
	# Create a vector of unit IDs (as characters)
	unit_ids <- sprintf("unit%02d", 1:N)
	time_vals <- 1:T

	# For each unit, decide the first treatment time:
	# With probability 0.5 the unit is never treated (first_treat = Inf);
	# Otherwise, sample a first treatment time uniformly between 2 and T.
	first_treat <- sapply(unit_ids, function(u) {
		if (runif(1) < 0.5) {
			Inf
		} else {
			if (R > 1) {
				sample(2:(R + 1), 1)
			} else {
				2
			}
		}
	})

	# Build panel data (each unit appears for every time period)
	df <- do.call(
		rbind,
		lapply(seq_along(unit_ids), function(i) {
			unit <- unit_ids[i]
			ft <- first_treat[i]
			data.frame(
				time = as.integer(time_vals), # must be integer
				unit = as.character(unit), # must be character
				treatment = as.integer(ifelse(time_vals >= ft, 1, 0)), # must be integer 0/1
				cov1 = rnorm(T),
				cov2 = runif(T),
				y = rnorm(T) # outcome (numeric)
			)
		})
	)

	# The idCohorts function automatically removes any unit that was treated in period 1.
	# (This is simulated by removing any row with time == 1 and treatment == 1.)
	df <- df[!(df$time == 1 & df$treatment == 1), ]

	# Order rows by unit then time
	df <- df[order(df$unit, df$time), ]
	rownames(df) <- NULL
	return(df)
}

# ------------------------------------------------------------------------------
# Helper function: generate a panel where every unit is eventually treated
# (i.e., NO never-treated units). Used to exercise the `allow_no_never_treated`
# code path.
# ------------------------------------------------------------------------------
generate_bad_panel_data <- function(N = 30, T = 10, seed = 123) {
	set.seed(seed)
	# Create a vector of unit IDs (as characters)
	unit_ids <- sprintf("unit%02d", 1:N)
	time_vals <- 1:T

	# For each unit, decide the first treatment time:
	# Always sample a first treatment time uniformly between 2 and T
	# No never-treated units.
	first_treat <- sapply(unit_ids, function(u) {
		sample(2:T, 1)
	})

	# Build panel data (each unit appears for every time period)
	df <- do.call(
		rbind,
		lapply(seq_along(unit_ids), function(i) {
			unit <- unit_ids[i]
			ft <- first_treat[i]
			data.frame(
				time = as.integer(time_vals), # must be integer
				unit = as.character(unit), # must be character
				treatment = as.integer(ifelse(time_vals >= ft, 1, 0)), # must be integer 0/1
				cov1 = rnorm(T),
				cov2 = runif(T),
				y = rnorm(T) # outcome (numeric)
			)
		})
	)

	# Order rows by unit then time
	df <- df[order(df$unit, df$time), ]
	rownames(df) <- NULL
	return(df)
}

# ------------------------------------------------------------------------------
# Helper function: minimal hand-crafted panel (N = 3, T = 3) with one
# never-treated unit, one unit treated at t=2, and one unit treated at t=3.
# `include_cov1 = TRUE` adds a single `cov1 = rnorm(T)` column (matches the
# pre-#91 fetwfe/betwfe inline versions); `include_cov1 = FALSE` omits it
# (matches the pre-#91 etwfe/twfeCovs versions). The `cov1` rnorm draws
# happen BEFORE the `y` draws so the seed stream is preserved.
# ------------------------------------------------------------------------------
generate_minimal_panel_data <- function(seed = 123, include_cov1 = TRUE) {
	set.seed(seed)

	N <- 3
	T <- 3
	# Create a vector of unit IDs (as characters)
	unit_ids <- sprintf("unit%02d", 1:N)
	time_vals <- 1:T

	# For each unit, decide the first treatment time:
	# one treated, one untreated unit.
	first_treat <- sapply(unit_ids, function(u) {
		if (u == "unit01") {
			Inf
		} else if (u == "unit02") {
			2
		} else {
			3
		}
	})

	# Build panel data (each unit appears for every time period)
	df <- do.call(
		rbind,
		lapply(seq_along(unit_ids), function(i) {
			unit <- unit_ids[i]
			ft <- first_treat[i]
			base <- data.frame(
				time = as.integer(time_vals), # must be integer
				unit = as.character(unit), # must be character
				treatment = as.integer(ifelse(time_vals >= ft, 1, 0)) # must be integer 0/1
			)
			if (include_cov1) {
				base$cov1 <- rnorm(T)
			}
			base$y <- rnorm(T) # outcome (numeric)
			base
		})
	)

	# The idCohorts function automatically removes any unit that was treated in period 1.
	# (This is simulated by removing any row with time == 1 and treatment == 1.)
	df <- df[!(df$time == 1 & df$treatment == 1), ]

	# Order rows by unit then time
	df <- df[order(df$unit, df$time), ]
	rownames(df) <- NULL
	return(df)
}
