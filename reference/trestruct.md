# Detect cryptic population structure in time trees

Estimates a partition of a time-scaled tree by contrasting coalescent
patterns.

## Usage

``` r
trestruct(
  tre,
  minCladeSize = 25,
  minOverlap = -Inf,
  nodeSupportValues = FALSE,
  nodeSupportThreshold = 95,
  nsim = 10000,
  level = 0.01,
  ncpu = 1,
  verbosity = 1,
  debugLevel = 0,
  levellb = 0.001,
  levelub = 0.1,
  res = 11
)
```

## Arguments

- tre:

  A tree of type ape::phylo. Must be rooted. If the tree has
  multifurcations, it will be converted to a binary tree before
  processing.

- minCladeSize:

  All clusters within partition must have at least this many tips.

- minOverlap:

  Threshold time overlap required to find splits in a clade.

- nodeSupportValues:

  Node support values such as produced by bootstrap or Bayesian
  credibility scores. Must be logical or vector with length equal to
  number of internal nodes in the tree. If nodeSupportValues = TRUE,
  then the function will get the information on node support from the
  tree. If numeric vector, these values should be between 0 and 100.

- nodeSupportThreshold:

  Threshold node support value between 0 and 100. Nodes with support
  lower than this threshold will not be tested.

- nsim:

  Number of simulations for computing null distribution of test
  statistics.

- level:

  Significance level for finding new split within a set of tips. Can
  also be NULL, in which case the optimal level is found according to
  the CH index (see details).

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

## References

Volz EM, Carsten W, Grad YH, Frost SDW, Dennis AM, Didelot X.
Identification of hidden population structure in time-scaled
phylogenies. Systematic Biology 2020; 69(5):884-896.

## Author

Erik M Volz

## Examples

``` r
tree <- ape::rcoal(50)
struct <-  trestruct( tree )
#> Finding splits under nodes: 51 
print(struct)
#> Call: 
#> .trestruct(tre = tre, minCladeSize = minCladeSize, minOverlap = minOverlap, 
#>     nodeSupportValues = nodeSupportValues, nodeSupportThreshold = nodeSupportThreshold, 
#>     nsim = nsim, level = level[1], ncpu = ncpu, verbosity = verbosity, 
#>     debugLevel = debugLevel, useNodeSupport = useNodeSupport, 
#>     tredat = tredat)
#> 
#> Significance level: 0.01 
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
