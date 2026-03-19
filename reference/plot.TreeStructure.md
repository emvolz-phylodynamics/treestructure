# Plot TreeStructure tree with cluster and partition variables

Plot TreeStructure tree with cluster and partition variables

## Usage

``` r
# S3 method for class 'TreeStructure'
plot(x, use_ggtree = TRUE, ...)
```

## Arguments

- x:

  A TreeStructure object

- use_ggtree:

  Toggle ggtree or ape plotting behavior

- ...:

  Additional arguments passed to ggtree or ape::plot.phylo

## Examples

``` r
#tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') )
# you can run the example below before plotting
#struc <- trestruct( tree )

#because it can take a minute or so to run treestructure, we will load it here
struc <- readRDS( system.file('struc_plot_example.rds', package='treestructure') )
#plot treestructure object

suppressWarnings(plot(struc))
```
