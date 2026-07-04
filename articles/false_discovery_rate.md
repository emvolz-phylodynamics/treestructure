# Calibrating cluster detection to a false discovery rate

## Why a false discovery rate?

`treestructure` partitions a tree by repeatedly testing, at each node,
whether a clade differs from the rest of the tree in its coalescent
pattern. The `level` argument sets the significance of each individual
test. Because the algorithm examines many candidate clades and retains
the most extreme one at every step, `level` behaves like a subjective
sensitivity dial: lowering it detects more clusters, but the number of
*spurious* clusters grows as well. On its own it is not an error rate.

Specifying a target **false discovery rate** with the `fdr` argument
makes the threshold interpretable, and is the default behaviour of
`trestruct`. We illustrate with the simulated tree from the
[introductory
vignette](https://emvolz-phylodynamics.github.io/treestructure/articles/treestructure.md):

``` r

library(ape)
library(treestructure)

tree <- read.tree(system.file("sim.nwk", package = "treestructure"))
s <- trestruct(tree, fdr = 0.05)
#> Tree has NA or duplicated tip labels. Adding a unique id.
#> Finding splits under nodes: 276 
#> Finding splits under nodes: 276 526 
#> Finding splits under nodes: 276 421 
#> Finding splits under nodes: 276 473
print(s)
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

The output reports the target FDR, a global structure test (below), and
the resulting clusters.

## What `fdr` controls

`fdr` is used by default; supplying `level` without `fdr` selects a
subjective threshold instead. The false discovery rate is a property of
the **whole tree**, not of any single test:

- If the tree has no structure at all — a single, unstructured
  coalescent — the probability that `treestructure` designates *any*
  cluster is `fdr`.
- If the tree does contain real structure, `fdr` bounds the expected
  fraction of spurious clusters among all the clusters designated.

Internally this is achieved by requiring the most extreme clade at each
node to clear a threshold that is adjusted for the number of candidate
clades examined there, so that the error is controlled across the whole
tree rather than one clade at a time.

The calibration is analytic and deterministic — identical inputs give
identical results, with no simulation.

### A test for any structure at all

Alongside the partition, the result carries a single test of the null
hypothesis that the tree has *no* structure, reported in the printed
output and available as:

``` r

s$global.test
#> $statistic
#> [1] 5.602629
#> 
#> $k
#> [1] 52
#> 
#> $p.value
#> [1] 1.097851e-06
```

`statistic` is the most extreme standardized clade contrast at the root,
`k` the number of candidate clades, and `p.value` a p-value for the
presence of any structure.

## The false discovery rate is controlled

Unstructured coalescent trees (produced by
[`ape::rcoal`](https://rdrr.io/pkg/ape/man/rtree.html)) contain no
structure, so the fraction of such trees in which `treestructure`
designates any cluster should not exceed the target `fdr`. The following
pre-computed check confirms this:

``` r

readRDS(system.file("vignette_fwer.rds", package = "treestructure"))
#>   target_fdr realized reps ntip
#> 1       0.05    0.020  200  100
#> 2       0.10    0.070  200  100
#> 3       0.20    0.085  200  100
```

The realized rate stays at or below the target. It is somewhat
conservative: the threshold treats every candidate clade as a separate
test, whereas nearby clades are statistically correlated, so the true
error is smaller than the correction assumes.

## Serially sampled (heterochronous) trees

The coalescent null underlying the test assumes that all tips were
sampled at the same time. When samples are collected over a period
comparable to the timescale of the tree — as with most genomic
epidemiology data — the test statistic acquires a small upward bias,
which can push the realized false discovery rate somewhat above the
target.

In practice the inflation is modest. The table below reports the
realized false discovery rate at a target of 0.05 for unstructured trees
simulated to match the estimated population dynamics and sampling times
of several real, densely time-sampled pathogen datasets. The
`heterochronicity_index` (0 for contemporaneous sampling, approaching 1
as sampling spans the tree) measures how serially sampled each dataset
is:

``` r

readRDS(system.file("vignette_hetero.rds", package = "treestructure"))
#>          dataset heterochronicity_index realized_FDR
#> 1     SARS-CoV-2                   0.87         0.02
#> 2 influenza H3N2                   0.95         0.10
#> 3   B. pertussis                   0.76         0.03
#> 4          Ebola                   0.86         0.04
```

Even the most heterochronous datasets stay within roughly twice the
target, and most are at or below it. `treestructure` detects
heterochronous sampling automatically and prints a note when it is
present; the clade whose most recent sample is oldest is the most
affected. On such trees the target `fdr` should be read as approximate,
particularly for deep clusters.

## Large trees: Bonferroni vs Benjamini-Hochberg

At each step the algorithm corrects for the number of candidate clades
it examines. By default this is a Bonferroni correction, which is simple
and conservative. On very large trees the Bonferroni bar needed to make
the *first* split grows slowly with the number of candidates, so
genuinely structured but individually-moderate clades can be missed.
Setting `split = "bh"` uses a Benjamini–Hochberg step-up instead, which
is less conservative — it recovers more such structure while still
controlling the false discovery rate:

``` r

trestruct(tree, fdr = 0.05, split = "bh")
```

Bonferroni remains the default; `"bh"` is worth trying on large trees,
or whenever you suspect moderate structure is being missed. (Since BH
always includes the Bonferroni decision as a special case, it never
designates fewer clusters.)

## Combining with node support

Whether a clade’s coalescent pattern differs significantly is a separate
question from whether the clade is well supported topologically. The two
criteria are complementary and can be combined: supply node support
values (for example bootstrap proportions or Bayesian posterior
probabilities) and a minimum threshold, and only sufficiently supported
clades are eligible to be designated as clusters.

``` r

# tree$node.label holds support values in [0, 100] (or [0, 1])
trestruct(tree, fdr = 0.05,
          nodeSupportValues = TRUE, nodeSupportThreshold = 95)
```

Requiring both a controlled false discovery rate and adequate support is
the most conservative option, and removes clusters designated at poorly
resolved nodes. The [node support
values](https://emvolz-phylodynamics.github.io/treestructure/articles/supportValues.md)
vignette develops this with a worked example.
