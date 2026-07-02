# Warm up the namespaces used by the tests BEFORE the first test runs, so testthat's
# per-test reproducible-state check never observes a newly-loaded namespace mid-run.
# (In environments where some Suggested packages use delayed S3-method registration,
# describing such a state change during the comparison can error.)
suppressMessages(invisible(
	trestruct(ape::rcoal(20), fdr = 0.2, minCladeSize = 5, verbosity = 0)
))
