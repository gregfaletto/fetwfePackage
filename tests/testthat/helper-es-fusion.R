# Shared test helpers for the event-study fusion transform (issues #40, #237).
#
# `.es_difference()` is the forward event-study differencing operator
# theta = D^{(2)}_ES beta, implemented directly from the paper's *prose*
# (Faletto 2025, Lemma event.study.sing.val.lem) and therefore independent of
# the package's inverse-transform matrix code. Recovering theta from a generated
# beta with it is a genuine cross-check, not a tautological round-trip through
# the same matrix that built beta.
#
# testthat sources helper-*.R before any test-*.R, so these are available to
# every test file. tests/testthat/test-event-study-fusion-40.R predates this
# helper and keeps its own local copies; this file does not change that.

.es_difference <- function(beta, first_inds, n_k, G) {
	theta <- numeric(length(beta))
	c1 <- first_inds[1]:(first_inds[1] + n_k[1] - 1L)
	theta[c1[1]] <- beta[c1[1]]
	if (n_k[1] >= 2L) {
		for (i in 2:n_k[1]) {
			theta[c1[i]] <- beta[c1[i]] - beta[c1[i] - 1L]
		}
	}
	if (G >= 2L) {
		for (k in 2:G) {
			for (e in 0:(n_k[k] - 1L)) {
				theta[first_inds[k] + e] <-
					beta[first_inds[k] + e] - beta[first_inds[k - 1L] + e]
			}
		}
	}
	theta
}

.es_setup <- function(G, T) {
	num_treats <- fetwfe:::getNumTreats(G = G, T = T)
	first_inds <- fetwfe:::getFirstInds(G = G, T = T)
	n_k <- c(diff(first_inds), num_treats - first_inds[G] + 1L)
	list(num_treats = num_treats, first_inds = first_inds, n_k = n_k)
}
