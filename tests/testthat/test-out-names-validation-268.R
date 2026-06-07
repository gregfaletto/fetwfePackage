library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #268: malformed `out_names` (passed through attgtToFetwfeDf()/etwfeToFetwfeDf()
# to .fetwfe_df_core()) must error with the documented four-name contract, not a
# cryptic internal error ("attempt to select less than one element", "$ operator is
# invalid for atomic vectors"). The named-character-vector case is the subtle one:
# it carries the four names yet the assembly indexes with `$`, invalid for an atomic
# vector -- so the guard must require a *list*, not merely the right names.
# ------------------------------------------------------------------------------
test_that("malformed out_names errors with the four-name contract message (#268)", {
	df <- data.frame(
		county = rep(1:3, each = 3),
		year = rep(2000:2002, times = 3),
		first.treat = rep(c(0L, 0L, 2001L), each = 3),
		lemp = as.numeric(1:9)
	)
	call_with <- function(on) {
		etwfeToFetwfeDf(
			data = df,
			yvar = "lemp",
			tvar = "year",
			idvar = "county",
			gvar = "first.treat",
			out_names = on
		)
	}
	msg <- "must be a list containing the four names"

	# Missing a required key.
	expect_error(
		call_with(list(time = "t", unit = "u", treatment = "tr")),
		msg,
		fixed = TRUE
	)
	# Mistyped key (`treat` instead of `treatment`).
	expect_error(
		call_with(list(time = "t", unit = "u", treat = "tr", response = "r")),
		msg,
		fixed = TRUE
	)
	# Named character VECTOR: carries all four names but is not a list.
	expect_error(
		call_with(c(time = "t", unit = "u", treatment = "tr", response = "r")),
		msg,
		fixed = TRUE
	)

	# A well-formed list must NOT raise the contract error (any later error -- e.g.
	# from the tiny panel -- is unrelated to out_names validation).
	good_err <- tryCatch(
		{
			call_with(list(
				time = "time_var",
				unit = "unit_var",
				treatment = "treatment",
				response = "response"
			))
			NA_character_
		},
		error = function(e) conditionMessage(e)
	)
	expect_false(isTRUE(grepl(msg, good_err, fixed = TRUE)))
})
