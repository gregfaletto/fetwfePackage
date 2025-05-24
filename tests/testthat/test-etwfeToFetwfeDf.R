# Tests for etwfeToFetwfeDf()
# Run with: testthat::test_file("tests/testthat/test-etwfeToFetwfeDf.R")

library(testthat)

## ---------------------------------------------------------------------------
## helper to create a minimal 3-unit × 3-year balanced panel
## ---------------------------------------------------------------------------
make_simple_df <- function() {
	grid <- expand.grid(
		county = 1:3,
		year = 2000:2002
	)

	transform(
		grid,
		first.treat = ifelse(county == 3, 2001, 0), # unit 3 treated from 2001
		lemp = rnorm(nrow(grid)),
		lpop = rnorm(nrow(grid))
	)
}

## ---------------------------------------------------------------------------
test_that("basic happy path works", {
	df <- make_simple_df()

	tidy <- etwfeToFetwfeDf(
		data = df,
		yvar = "lemp",
		tvar = "year",
		idvar = "county",
		gvar = "first.treat",
		covars = c("lpop")
	)

	expect_s3_class(tidy, "data.frame")
	expect_named(
		tidy,
		c("time_var", "unit_var", "treatment", "response", "lpop")
	)

	expect_true(all(tidy$time_var == as.integer(tidy$time_var)))
	expect_true(all(tidy$treatment %in% 0:1))

	## Unit 3 should switch on in 2001 and stay on
	unit3 <- subset(tidy, unit_var == "3")
	expect_equal(unit3$treatment, c(0L, 1L, 1L))
})

## ---------------------------------------------------------------------------
test_that("duplicate (id,time) rows error", {
	df <- make_simple_df()
	dup <- rbind(df, df[1, ]) # add a duplicate row

	expect_error(
		etwfeToFetwfeDf(dup, "lemp", "year", "county", "first.treat"),
		"pair must appear at most once"
	)
})

## ---------------------------------------------------------------------------
test_that("non-constant gvar within unit triggers error", {
	df <- make_simple_df()
	## make unit 2's gvar vary over time
	df$first.treat[df$county == 2 & df$year == 2002] <- 2002

	expect_error(
		etwfeToFetwfeDf(df, "lemp", "year", "county", "first.treat"),
		"irreversible treatment"
	)
})

## ---------------------------------------------------------------------------
test_that("dropping first-period treated units works", {
	df <- make_simple_df()
	df$first.treat[df$county == 1] <- 2000 # unit 1 treated from first period

	tidy <- etwfeToFetwfeDf(df, "lemp", "year", "county", "first.treat")
	expect_false(any(tidy$unit_var == "1"))
})

## ---------------------------------------------------------------------------
test_that("missing columns are caught", {
	df <- make_simple_df()[, -which(names(make_simple_df()) == "lemp")]

	expect_error(
		etwfeToFetwfeDf(df, "lemp", "year", "county", "first.treat"),
		"Column\\(s\\) not found"
	)
})

## ---------------------------------------------------------------------------
test_that("custom out_names are respected", {
	df <- make_simple_df()

	custom_names <- list(
		time = "yr",
		unit = "id",
		treatment = "D",
		response = "outcome"
	)

	tidy <- etwfeToFetwfeDf(
		data = df,
		yvar = "lemp",
		tvar = "year",
		idvar = "county",
		gvar = "first.treat",
		out_names = custom_names
	)

	expect_named(tidy, unlist(custom_names, use.names = FALSE))
	expect_identical(names(tidy)[1:4], unlist(custom_names, use.names = FALSE))
})
