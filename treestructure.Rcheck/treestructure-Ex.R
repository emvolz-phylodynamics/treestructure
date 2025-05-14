pkgname <- "treestructure"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('treestructure')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("addtips")
### * addtips

flush(stderr()); flush(stdout())

### Name: addtips
### Title: Compare and add tips into new treestructure object
### Aliases: addtips

### ** Examples

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

# assign structure to the previously missing tips
(s <- treestructure::addtips( s0, tr ))



cleanEx()
nameEx("plot.TreeStructure")
### * plot.TreeStructure

flush(stderr()); flush(stdout())

### Name: plot.TreeStructure
### Title: Plot TreeStructure tree with cluster and partition variables
### Aliases: plot.TreeStructure

### ** Examples


tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') )

(struc <- trestruct( tree ))

#plot treestructure object

plot(struc)



cleanEx()
nameEx("treestructure.test")
### * treestructure.test

flush(stderr()); flush(stdout())

### Name: treestructure.test
### Title: Test treestructure hypothesis
### Aliases: treestructure.test

### ** Examples

tree <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') )

(struc <- trestruct( tree ))

#run the test

results <- treestructure.test(tree, x = struc$clusterSets[[1]],
                              y = struc$clusterSets[[2]])

print(results)



cleanEx()
nameEx("trestruct")
### * trestruct

flush(stderr()); flush(stdout())

### Name: trestruct
### Title: Detect cryptic population structure in time trees
### Aliases: trestruct

### ** Examples

tree <- ape::rcoal(50)
struct <-  trestruct( tree )
print(struct)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
