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
#> Finding splits under nodes: 121 142 
#> Call: 
#> treestructure::trestruct(tre = tr0)
#> 
#> Target FDR: 0.2 (bonferroni correction)
#> Global structure test: max|z| = 3.34 over 24 candidates, p = 0.02
#> Number of clusters: 3 
#> Number of partitions: 2 
#> NOTE: heterochronous sampling detected (overlap index 0.86). Under the coalescent null this
#>       can modestly inflate the realized FDR (up to ~2x at target 5% in our real-data study);
#>       see vignette("treestructure"). Deep clades whose most recent sample is old are most affected.
#> Number of taxa in each cluster:
#> 
#>  1  2  3 
#> 78 28 14 
#> Number of taxa in each partition:
#> 
#>  1  2 
#> 78 42 
#> ...
#> For complete data, use `as.data.frame(...)` 

# assign structure to the previously missing tips
(s <- treestructure::addtips( s0, tr ))
#> Call: 
#> treestructure::trestruct(tre = tr0)
#> 
#> Target FDR: 0.2 (bonferroni correction)
#> Global structure test: max|z| = 3.34 over 24 candidates, p = 0.02
#> Number of clusters: 3 
#> Number of partitions: 2 
#> NOTE: heterochronous sampling detected (overlap index 0.86). Under the coalescent null this
#>       can modestly inflate the realized FDR (up to ~2x at target 5% in our real-data study);
#>       see vignette("treestructure"). Deep clades whose most recent sample is old are most affected.
#> Number of taxa in each cluster:
#> 
#>   1   2   3 
#> 100  36  14 
#> Number of taxa in each partition:
#> 
#>   1   2 
#> 100  50 
#> ...
#> For complete data, use `as.data.frame(...)` 
```
