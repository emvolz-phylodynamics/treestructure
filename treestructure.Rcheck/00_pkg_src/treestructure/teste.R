
library(treestructure)
library(ape)

tr1 <- read.tree("~/Downloads/seqs2_adj1_up_2020_03_15.fa.tree")
tr2 <- read.tree("~/Downloads/seqs1_up_2020_02_29.fa.tree")



mltr2_outl_rm <- readRDS( system.file('mltr2_outl_rm_sc2_feb2020.rds',
                                      package='treestructure') )

na_nodes <- which(mltr2_outl_rm$node.label == "")

na_tips <- numeric()
for (node in na_nodes) {
  na_tips <- c(na_tips, tips(mltr2_outl_rm, node))
}
na_tips <- unique(na_tips)




timetr2_phylo <- readRDS( system.file('timetr2_phylo_sc2_feb2020.rds',
                                      package='treestructure') )

tail(mltr2_outl_rm$edge)
mltr2_outl_rm$tip.label[899]
