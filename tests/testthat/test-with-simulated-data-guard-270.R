library(testthat)
library(fetwfe)

# ------------------------------------------------------------------------------
# #270: the FETWFE_simulated class guard moved from each *WithSimulatedData wrapper
# body into the shared .unpack_simulated_obj() helper. The raised error must still
# (a) carry the same message and (b) be reported at the CALLING wrapper, not the
# internal helper -- the helper uses stop(simpleError(..., call = sys.call(-1))) so
# the interactive `Error in <call>:` prefix matches the pre-#270 behavior and does
# not leak the helper name.
# ------------------------------------------------------------------------------
test_that("*WithSimulatedData reject non-FETWFE_simulated at the caller, not the helper (#270)", {
	bad <- list(a = 1)
	chk <- function(e, wrapper_name) {
		expect_s3_class(e, "error")
		expect_match(
			conditionMessage(e),
			"must be an object of class 'FETWFE_simulated'",
			fixed = TRUE
		)
		expect_identical(as.character(conditionCall(e)[[1]]), wrapper_name)
	}
	chk(
		tryCatch(etwfeWithSimulatedData(bad), error = identity),
		"etwfeWithSimulatedData"
	)
	chk(
		tryCatch(twfeCovsWithSimulatedData(bad), error = identity),
		"twfeCovsWithSimulatedData"
	)
	chk(
		tryCatch(fetwfeWithSimulatedData(bad), error = identity),
		"fetwfeWithSimulatedData"
	)
	chk(
		tryCatch(betwfeWithSimulatedData(bad), error = identity),
		"betwfeWithSimulatedData"
	)
})
