# Compare and add tips into new treestructure object

Compares a new input tree to an old treestructure fit and merges tips
into a new treestructure object. Tips in the new tree that are not in
the new treestructure will be merged. Merging is carried out based on a
phylogenetic criterion. The new tips are added to the cluster which
shares its MRCA (most recent common ancestor).

## Usage

``` r
addtips(trst, tre)
```

## Arguments

- trst:

  Original treestructure fit that that will be updated.

- tre:

  A new tree (ape::phylo) which may contain samples not in trst. This
  tree must be rooted, but does not need to be time-scaled or binary.

## Value

A new treestructure fit.

## Author

Erik Volz

## Examples

``` r
set.seed(072023)
# simulate two trees and bind them to simulate structure
tr1 <- ape::rcoal( 50 )
tr2 <- ape::rcoal( 100 )
tr1$tip.label <- gsub(tr1$tip.label, patt = 't', rep = 's')
tr1$edge.length <- tr1$edge.length*.5
tr1$root.edge <- 1
tr2$root.edge <- 1
tr <- ape::bind.tree(tr1, tr2, position = .5 ) |> ape::multi2di()

# subsample the tree to simulating missing tips and estimate structure
ex <- sample( tr$tip.label, size = 30, replace = FALSE)
tr0 <- ape::drop.tip( tr, ex )
(s0 <- treestructure::trestruct( tr0 ))
#> Finding splits under nodes: 121 
#> Finding splits under nodes: 121 163 
#> Call: 
#> .trestruct(tre = tre, minCladeSize = minCladeSize, minOverlap = minOverlap, 
#>     nodeSupportValues = nodeSupportValues, nodeSupportThreshold = nodeSupportThreshold, 
#>     nsim = nsim, level = level[1], ncpu = ncpu, verbosity = verbosity, 
#>     debugLevel = debugLevel, useNodeSupport = useNodeSupport, 
#>     tredat = tredat)
#> 
#> Significance level: 0.01 
#> Number of clusters: 2 
#> Number of partitions: 2 
#> Number of taxa in each cluster:
#> 
#>  1  2 
#> 42 78 
#> Number of taxa in each partition:
#> 
#>  1  2 
#> 42 78 
#> ...
#> For complete data, use `as.data.frame(...)` 

# assign structure to the previously missing tips
(s <- treestructure::addtips( s0, tr ))
#> Call: 
#> .trestruct(tre = tre, minCladeSize = minCladeSize, minOverlap = minOverlap, 
#>     nodeSupportValues = nodeSupportValues, nodeSupportThreshold = nodeSupportThreshold, 
#>     nsim = nsim, level = level[1], ncpu = ncpu, verbosity = verbosity, 
#>     debugLevel = debugLevel, useNodeSupport = useNodeSupport, 
#>     tredat = tredat)
#> 
#> Significance level: 0.01 
#> Number of clusters: 2 
#> Number of partitions: 2 
#> Number of taxa in each cluster:
#> 
#>   1   2 
#>  50 100 
#> Number of taxa in each partition:
#> 
#>   1   2 
#>  50 100 
#> ...
#> For complete data, use `as.data.frame(...)` 
```
