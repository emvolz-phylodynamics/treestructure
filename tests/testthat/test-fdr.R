# Tests for the FDR-calibrated partition (trestruct( ..., fdr = )), the global-null
# test, the heterochronous diagnostic, and composition with node support.

test_that("fdr mode returns a TreeStructure with fdr / global.test / hetero fields", {
	set.seed(1)
	tr <- ape::rcoal(80)
	s <- trestruct(tr, fdr = 0.05, minCladeSize = 15, verbosity = 0)
	expect_s3_class(s, "TreeStructure")
	expect_equal(s$fdr, 0.05)
	expect_false(is.null(s$global.test))
	expect_gte(s$global.test$p.value, 0)
	expect_lte(s$global.test$p.value, 1)
	expect_false(is.null(s$hetero))
})

test_that("fdr overrides level (with a message) and records the target", {
	set.seed(2)
	tr <- ape::rcoal(80)
	expect_message(s <- trestruct(tr, fdr = 0.05, level = 0.01, minCladeSize = 15, verbosity = 1),
		regexp = "false discovery rate", ignore.case = TRUE)
	expect_equal(s$fdr, 0.05)
})

test_that("invalid fdr values are rejected", {
	tr <- ape::rcoal(40)
	expect_error(trestruct(tr, fdr = 1.5, verbosity = 0))
	expect_error(trestruct(tr, fdr = 0, verbosity = 0))
	expect_error(trestruct(tr, fdr = -0.1, verbosity = 0))
})

test_that("the analytic fdr partition is deterministic", {
	set.seed(3)
	tr <- ape::rcoal(90)
	a <- trestruct(tr, fdr = 0.1, minCladeSize = 15, verbosity = 0)
	b <- trestruct(tr, fdr = 0.1, minCladeSize = 15, verbosity = 0)
	expect_identical(as.integer(a$clustering), as.integer(b$clustering))
})

test_that("fdr is the default threshold method", {
	set.seed(7)
	tr <- ape::rcoal(80)
	s <- trestruct(tr, minCladeSize = 15, verbosity = 0)   # neither fdr nor level supplied
	expect_equal(s$fdr, 0.2)
})

test_that("an explicit level (without fdr) selects level mode", {
	set.seed(4)
	tr <- ape::rcoal(80)
	s <- trestruct(tr, level = 0.01, minCladeSize = 15, verbosity = 0)
	expect_null(s$fdr)
})

test_that("isochronous trees are not flagged as heterochronous", {
	set.seed(5)
	expect_lt(treestructure:::.hetero_index(ape::rcoal(100)), 0.1)
})

test_that("node support composes with fdr: all-unsupported yields no structure", {
	set.seed(6)
	tr <- ape::rcoal(80)
	tr$node.label <- rep("0", tr$Nnode)                       # posterior 0 -> nothing clears support
	s <- trestruct(tr, fdr = 0.1, minCladeSize = 15,
		nodeSupportValues = TRUE, nodeSupportThreshold = 95, verbosity = 0)
	expect_equal(nlevels(s$clustering), 1L)
})

test_that("split selects the correction; BH is no less permissive than Bonferroni", {
	tr <- ape::read.tree(system.file("sim.nwk", package = "treestructure"))
	b <- trestruct(tr, fdr = 0.1, minCladeSize = 10, split = "bonferroni", verbosity = 0)
	h <- trestruct(tr, fdr = 0.1, minCladeSize = 10, split = "bh", verbosity = 0)
	expect_equal(b$split, "bonferroni")
	expect_equal(h$split, "bh")
	expect_gte(nlevels(h$clustering), nlevels(b$clustering))
	expect_equal(trestruct(tr, fdr = 0.1, minCladeSize = 10, verbosity = 0)$split, "bonferroni")  # default
	expect_error(trestruct(tr, fdr = 0.1, split = "nope", verbosity = 0))
})

test_that("split='bh' still controls the global-null FWER", {
	skip_on_cran()
	set.seed(43)
	anystruct <- vapply(seq_len(30), function(i)
		nlevels(trestruct(ape::rcoal(80), fdr = 0.1, minCladeSize = 15, split = "bh", verbosity = 0)$clustering) > 1,
		logical(1))
	expect_lt(mean(anystruct), 0.35)
})

test_that("global-null FWER is approximately controlled on isochronous trees", {
	skip_on_cran()
	set.seed(42)
	reps <- 30
	anystruct <- vapply(seq_len(reps), function(i) {
		tr <- ape::rcoal(80)
		nlevels(trestruct(tr, fdr = 0.1, minCladeSize = 15, verbosity = 0)$clustering) > 1
	}, logical(1))
	expect_lt(mean(anystruct), 0.35)                          # target 0.10; loose upper bound
})
