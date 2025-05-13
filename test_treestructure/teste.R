#script to test that after fixing the global variables I still get the same
#answer using treestructure
library(ape)
library(treestructure)

( tre <- ape::read.tree( system.file('sim.nwk', package = 'treestructure') ) )

tredat <- .tredat(tre)

s <- trestruct( tre )

saveRDS(s, "treestr_code_as_it_is.RDS")

test_s <- readRDS("test_treestructure/treestr_code_as_it_is.RDS")
