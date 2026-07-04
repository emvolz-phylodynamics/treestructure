# treestructure applied to structured coalescent simulation

## Structured coalescent simulation

This example shows the function `trestruct` applied to a simulated
structured coalescent tree that includes samples from a large constant
size population and samples from three small “outbreaks” which are
growing exponentially. These simulations were generated with the
[phydynR package](https://emvolz-phylodynamics.github.io/phydynR/).

``` r

library(ape)
library(treestructure)
library(ggtree)
```

Load the tree:

``` r

( tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') ) )
#> 
#> Phylogenetic tree with 275 tips and 274 internal nodes.
#> 
#> Tip labels:
#>   1, 1, 1, 1, 1, 1, ...
#> 
#> Rooted; includes branch length(s).
```

Note that the tip labels corresponds to the deme of each sample. ‘1’ is
the constant size reservoir, and ‘0’ is the exponentially growing deme.

This will run the `treestructure` algorithm under default settings:

``` r

s <- trestruct( tree )
#> Tree has NA or duplicated tip labels. Adding a unique id.
#> Finding splits under nodes: 276 
#> Finding splits under nodes: 276 526 
#> Finding splits under nodes: 276 421 526 530 
#> Finding splits under nodes: 276 421 422 473
```

The `treestructure` algorithm searches from root to tips of the
phylogenetic tree. The message printed on screen,
`Finding splits under nodes`, reports which nodes the algorithm is
currently testing.

By default `trestruct` calibrates the split threshold to a target false
discovery rate (`fdr = 0.2`), so that the expected proportion of
spuriously designated clusters is controlled across the whole tree. The
printed output reports this target and a global test for whether the
tree has any structure at all:

``` r

print(s)
#> Call: 
#> trestruct(tre = tree)
#> 
#> Target FDR: 0.2 (bonferroni correction)
#> Global structure test: max|z| = 5.60 over 52 candidates, p = 1.1e-06
#> Number of clusters: 6 
#> Number of partitions: 3 
#> Number of taxa in each cluster:
#> 
#>   1   2   3   4   5   6 
#>  12  13 200  10  15  25 
#> Number of taxa in each partition:
#> 
#>   1   2   3 
#>  12  53 210 
#> ...
#> For complete data, use `as.data.frame(...)`
```

### Plotting results

The default plotting behavior uses the `ggtree` package if available.

``` r

plot(s)  + ggtree::geom_tiplab()
```

![](treestructure_files/figure-html/unnamed-chunk-5-1.png)

If not, or if desired, `ape` plots are available

``` r

plot( s, use_ggtree = FALSE )
```

![](treestructure_files/figure-html/unnamed-chunk-6-1.png)

For subsequent analysis, you may want to turn the `treestructure` result
into a dataframe:

``` r

structureData <- as.data.frame( s )
head( structureData )
#>   taxon cluster partition
#> 1   1_1       3         3
#> 2   2_1       3         3
#> 3   3_1       3         3
#> 4   4_1       3         3
#> 5   5_1       3         3
#> 6   6_1       3         3
```

Each cluster and partition assignment is stored as a factor. You could
use `split` to get a data frame for each partition. Suppose we want a
tree corresponding to partition 1:

``` r

with ( structureData,
       ape::keep.tip(s$tree, taxon[ partition==1 ] )
       ) -> partition1
partition1
#> 
#> Phylogenetic tree with 12 tips and 11 internal nodes.
#> 
#> Tip labels:
#>   264_0, 265_0, 266_0, 267_0, 268_0, 269_0, ...
#> 
#> Rooted; includes branch length(s).
plot(partition1)
```

![](treestructure_files/figure-html/unnamed-chunk-8-1.png)

## Parameter choice and number of clusters

Two things have the largest influence on the number of clusters.

1.  `minCladeSize` controls the smallest allowed cluster in terms of the
    number of tips; the default is 10. With a smaller value, smaller
    clusters may be detected, but computation time increases. For
    example, structure that exists in clades smaller than `minCladeSize`
    will not be detected regardless of the threshold, so lowering it can
    reveal finer structure:

``` r

trestruct( tree, minCladeSize = 5 )
#> Tree has NA or duplicated tip labels. Adding a unique id.
#> Finding splits under nodes: 276 
#> Finding splits under nodes: 276 526 
#> Finding splits under nodes: 276 421 
#> Finding splits under nodes: 276 473
#> Call: 
#> trestruct(tre = tree, minCladeSize = 5)
#> 
#> Target FDR: 0.2 (bonferroni correction)
#> Global structure test: max|z| = 5.60 over 109 candidates, p = 2.3e-06
#> Number of clusters: 4 
#> Number of partitions: 2 
#> Number of taxa in each cluster:
#> 
#>   1   2   3   4 
#>  25  25 200  25 
#> Number of taxa in each partition:
#> 
#>   1   2 
#>  75 200 
#> ...
#> For complete data, use `as.data.frame(...)`
```

2.  The split threshold — how readily a clade is designated a new
    cluster. By default this is set by a false discovery rate (`fdr`);
    alternatively it can be set by a subjective significance `level`, or
    that `level` chosen automatically with the CH index.

#### False discovery rate (default)

By default the threshold is calibrated to a target false discovery rate,
so that `fdr` has an interpretable, whole-tree error-rate meaning:

``` r

trestruct( tree, fdr = 0.05 )
#> Tree has NA or duplicated tip labels. Adding a unique id.
#> Finding splits under nodes: 276 
#> Finding splits under nodes: 276 526 
#> Finding splits under nodes: 276 421 
#> Finding splits under nodes: 276 473
#> Call: 
#> trestruct(tre = tree, fdr = 0.05)
#> 
#> Target FDR: 0.05 (bonferroni correction)
#> Global structure test: max|z| = 5.60 over 52 candidates, p = 1.1e-06
#> Number of clusters: 4 
#> Number of partitions: 2 
#> Number of taxa in each cluster:
#> 
#>   1   2   3   4 
#>  25  25 200  25 
#> Number of taxa in each partition:
#> 
#>   1   2 
#>  75 200 
#> ...
#> For complete data, use `as.data.frame(...)`
```

The [false discovery
rate](https://emvolz-phylodynamics.github.io/treestructure/articles/false_discovery_rate.md)
vignette explains what the rate controls, the global test for any
structure, and the effect of serial (heterochronous) sampling.

#### Subjective significance level

Supplying a `level` (without `fdr`) uses it instead as a fixed per-test
significance. This is a subjective clustering threshold rather than an
error rate: raising it detects more clusters but also increases the
false positive rate.

``` r

trestruct( tree, level = 0.05 )
#> Tree has NA or duplicated tip labels. Adding a unique id.
#> Finding splits under nodes: 276 
#> Finding splits under nodes: 276 526 
#> Finding splits under nodes: 276 421 
#> Finding splits under nodes: 276 473
#> Call: 
#> trestruct(tre = tree, level = 0.05)
#> 
#> Significance level: 0.05 
#> Global structure test: max|z| = 5.60 over 52 candidates, p = 1.1e-06
#> Number of clusters: 4 
#> Number of partitions: 2 
#> Number of taxa in each cluster:
#> 
#>   1   2   3   4 
#>  25  25 200  25 
#> Number of taxa in each partition:
#> 
#>   1   2 
#>  75 200 
#> ...
#> For complete data, use `as.data.frame(...)`
```

The best value of `level` depends on your application. One way to choose
it is to use additional data associated with each sample, selecting the
`level` that gives clusters explaining the most variance in a variable
of interest (e.g. using the cluster as a factor in an ANOVA).

#### CH index

Alternatively, in the absence of any additional data, the
`treestructure` package supports using the [CH
index](https://en.wikipedia.org/wiki/Calinski%E2%80%93Harabasz_index) to
compare different `level`s. This statistic is based on the ratio of the
between-cluster and within-cluster variance of the time of each node
(distance from the root) and returns the `level` such that this ratio is
maximized. If you wish to use the CH index, pass `level = NULL` to
`trestruct`, and read documentation for the `levellb`, `levelub`, and
`res` parameters.
