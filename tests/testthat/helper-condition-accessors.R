# Degrading accessors for condition objects captured with tryCatch().
#
# WHY THESE EXIST, rather than expect_error()-and-reuse. A cell whose mutation
# makes the call SUCCEED has to fail its loop cleanly rather than abort it.
# expect_error() on a non-erroring call records its failure and returns the
# expression's VALUE, and conditionMessage() on that value raises -- which ends
# the whole test_that() block and silently skips every later cell. The mutation
# battery would still show the row as detected, because the aborted block
# reports a failure either way, so the loss is invisible.
#
# These two degrade to a readable value instead of raising, so every assertion
# in a loop still runs and a red cell names itself.
#
# Shared by test-gen-coefs-validator-429.R, test-gencoefs-retry-cap-436.R and
# test-genassignments-retry-cap-436.R. Extracted here (#436) once the third
# byte-identical copy arrived; testthat sources helper-*.R before the tests.

msg_of <- function(x) {
	if (inherits(x, "condition")) {
		conditionMessage(x)
	} else {
		paste0("<no error; returned ", class(x)[1], ">")
	}
}

call_of <- function(x) {
	if (inherits(x, "condition")) {
		conditionCall(x)
	} else {
		quote(NO_ERROR_RAISED)
	}
}
