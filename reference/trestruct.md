# Detect cryptic population structure in time trees

Estimates a partition of a time-scaled tree by contrasting coalescent
patterns.

## Usage

``` r
trestruct(
  tre,
  fdr = 0.2,
  level = 0.01,
  minCladeSize = 10,
  nodeSupportValues = FALSE,
  nodeSupportThreshold = 95,
  minOverlap = -Inf,
  nsim = 10000,
  ncpu = 1,
  verbosity = 1,
  debugLevel = 0,
  levellb = 0.001,
  levelub = 0.1,
  res = 11,
  method = "analytic",
  split = c("bonferroni", "bh")
)
```

## Arguments

- tre:

  A tree of type ape::phylo. Must be rooted. If the tree has
  multifurcations, it will be converted to a binary tree before
  processing.

- fdr:

  Target false discovery rate for detected structure, a number in (0,1).
  This is the default way of choosing the split threshold (`fdr = 0.2`):
  the threshold at each scan is calibrated so that the whole-tree false
  discovery rate is controlled at this level (see details for the
  precise meaning). It is analytic and requires no simulation. An
  explicitly supplied `level` takes precedence unless `fdr` is also
  given; set `fdr = NULL` to use `level` instead.

- level:

  Significance level for finding a new split within a set of tips. Used
  when `fdr = NULL`, or when `level` is supplied without an explicit
  `fdr`. This is a subjective clustering threshold, not an error rate;
  prefer `fdr`. Can also be NULL, in which case the optimal level is
  found according to the CH index (see details).

- minCladeSize:

  All clusters within partition must have at least this many tips.

- nodeSupportValues:

  Node support values such as produced by bootstrap or Bayesian
  credibility scores. Must be logical or vector with length equal to
  number of internal nodes in the tree. If nodeSupportValues = TRUE,
  then the function will get the information on node support from the
  tree. If numeric vector, these values should be between 0 and 100.

- nodeSupportThreshold:

  Threshold node support value between 0 and 100. Nodes with support
  lower than this threshold will not be tested.

- minOverlap:

  Threshold time overlap required to find splits in a clade.

- nsim:

  Number of simulations for computing null distribution of test
  statistics. Only used when `method = 'sim'`.

- ncpu:

  If \> 1 will compute statistics in parallel using multiple CPUs.

- verbosity:

  If \> 0 will print information about progress of the algorithm.

- debugLevel:

  If \> 0 will produce additional data in return value.

- levellb:

  If optimizing the \`level\` parameter, this is the lower bound for the
  search.

- levelub:

  If optimizing the \`level\` parameter, this is the upper bound for the
  search.

- res:

  If optimizing the \`level\` parameter, this is the number of values to
  test.

- method:

  How the coalescent null of the rank-sum statistic is characterised at
  each test. `'analytic'` (the default) computes the exact mean and
  variance by a deterministic recursion (fast, and deterministic so
  results are reproducible); `'sim'` uses Monte-Carlo simulation with
  `nsim` replicates (the original behaviour).

- split:

  Multiple-testing correction applied at each scan in `fdr` mode.
  `'bonferroni'` (the default) uses a per-scan Bonferroni bound; `'bh'`
  uses a Benjamini-Hochberg step-up, which is less conservative and
  retains more power on large trees with abundant moderate structure
  (where the Bonferroni bound to make the first split can grow with the
  tree size). Ignored when a subjective `level` is used.

## Value

A TreeStructure object which includes cluster and partition assignment
for each tip of the tree.

## Details

Estimates a partition of a time-scaled tree by contrasting coalescent
patterns. The algorithm is premised on a Kingman coalescent null
hypothesis for the ordering of node heights when contrasting two clades,
and a test statistic is formulated based on the rank sum of node times
in the tree. If node support values are available (as computed by
bootstrap procedures), the method can optionally exclude designation of
structure on poorly supported nodes. The method will not designate
structure on nodes with zero branch length relative to their immediate
ancestor. The significance level for detecting significant partitions of
the tree can be provided, or a range of values can be examined. The [CH
index](https://en.wikipedia.org/wiki/Calinski-Harabasz_index) based on
within- and between-cluster variance in node heights can be used to
select a significance level if none is provided.

**Calibrating to a false discovery rate (the default).** By default, and
whenever `fdr` is supplied, the split threshold is calibrated to a
target false discovery rate rather than to a subjective `level`. At each
scan the algorithm splits at the most extreme eligible candidate clade
only if its standardised statistic clears the Bonferroni threshold
\\\Phi^{-1}(1 - fdr/(2k))\\, where \\k\\ is the number of *eligible*
candidates in that scan; candidates excluded by `minCladeSize`, node
support, or time overlap do not count towards \\k\\. The `fdr` refers to
the **whole tree**, not an individual scan or clade: under the global
null of one unstructured coalescent the probability of designating *any*
structure equals `fdr`, and when real structure is present `fdr` bounds
the expected fraction of spurious splits among all splits. This
whole-tree guarantee is obtained by controlling each scan; it is not a
per-clade p-value. Unlike `level`, the analytic default requires no
simulation and is deterministic.

The returned object also carries a global-null test in `$global.test`
(the root-scan \\\max\|z\|\\, the number of candidates \\k\\, and a
Bonferroni p-value for the presence of *any* structure), and, in `fdr`
mode, a heterochronous-sampling diagnostic in `$hetero`. Serially
sampled (heterochronous) trees can modestly inflate the realised FDR;
this is reported and discussed in the package vignette.

## References

Volz EM, Carsten W, Grad YH, Frost SDW, Dennis AM, Didelot X.
Identification of hidden population structure in time-scaled
phylogenies. Systematic Biology 2020; 69(5):884-896.

## Author

Erik M Volz

## Examples

``` r
tree <- ape::rcoal(50)
# subjective clustering threshold (default):
struct <-  trestruct( tree )
#> Finding splits under nodes: 51 
print(struct)
#> Call: 
#> trestruct(tre = tree)
#> 
#> Target FDR: 0.2 (bonferroni correction)
#> Global structure test: max|z| = 1.92 over 9 candidates, p = 0.499
#> Number of clusters: 1 
#> Number of partitions: 1 
#> Number of taxa in each cluster:
#> 
#>  1 
#> 50 
#> Number of taxa in each partition:
#> 
#>  1 
#> 50 
#> ...
#> For complete data, use `as.data.frame(...)` 
# calibrate the threshold to a target false discovery rate instead:
struct_fdr <- trestruct( tree, fdr = 0.05 )
#> Finding splits under nodes: 51 
print(struct_fdr)
#> Call: 
#> trestruct(tre = tree, fdr = 0.05)
#> 
#> Target FDR: 0.05 (bonferroni correction)
#> Global structure test: max|z| = 1.92 over 9 candidates, p = 0.499
#> Number of clusters: 1 
#> Number of partitions: 1 
#> Number of taxa in each cluster:
#> 
#>  1 
#> 50 
#> Number of taxa in each partition:
#> 
#>  1 
#> 50 
#> ...
#> For complete data, use `as.data.frame(...)` 
```
